#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define STENCIL_SIZE 500
//The STENCIL_SIZE value can be modified to 20, 200, 300 or 500 for example. When doing this, it is recommended to comment the display_matrix function, and remove the call of this function in stencil_init_mpi()


typedef float stencil_t;

/** conduction coeff used in computation */
static const stencil_t alpha = 0.02;

/** threshold for convergence */
static const stencil_t epsilon = 0.0001;

/** max number of steps */
static const int stencil_max_steps = 100000;

static stencil_t*values = NULL;
static stencil_t*prev_values = NULL;
static stencil_t* global_values = NULL;

static int* sendcounts;
static int* displs;

static int size_x = STENCIL_SIZE;
static int size_y = STENCIL_SIZE;


static int rank;
static int nprocs;

static int local_ny;     // number of real rows of the process. size_x keeps being global. size_y no (local_ny)




/** display the global matrix from a pointer, after using gatherv*/
static void stencil_display_global(stencil_t *buf)
{
    int x, y;
    for (y = 0; y < size_y; y++) {
        for (x = 0; x < size_x; x++) {
            printf("%8.5g ", buf[x + size_x * y]);
        }
        printf("\n");
    }
}


/** init stencil values to 0, borders to non-zero */
static void stencil_init_mpi(void)
{
    int base = size_y / nprocs;
    int rest = size_y % nprocs;

    local_ny = base + (rank < rest ? 1 : 0);   //the first 'rest' processes receives 1 more row

    //Step 2. ok, now EVERY process reserves local memory (with halos, that's why +2)
    values = malloc((local_ny + 2) * size_x * sizeof(stencil_t));
    prev_values = malloc((local_ny + 2) * size_x * sizeof(stencil_t));
        

    if(rank==0)
    {
        global_values = malloc(size_y * size_x * sizeof(stencil_t));

        int x, y;
        for(x = 0; x < size_x; x++)   //"inside" the table
        {
          for(y = 0; y < size_y; y++)
            {
              global_values[x + size_x * y] = 0.0;
            }
        }
        for(x = 0; x < size_x; x++) //top and bottom borders
        {
          global_values[x + size_x * 0] = x;
          global_values[x + size_x * (size_y - 1)] = size_x - x - 1;
        }
        for(y = 0; y < size_y; y++)   //right and legt borders
        {
          global_values[0 + size_x * y] = y;
          global_values[size_x - 1 + size_x * y] = size_y - y - 1;
        }

        //2 vectors to make the distribuiton of rows with scatterv later.
        sendcounts = malloc(nprocs * sizeof(int));
        displs     = malloc(nprocs * sizeof(int));

        int offset = 0;
        for (int p = 0; p < nprocs; p++) {
            int rows_p = base + (p < rest ? 1 : 0);
            sendcounts[p] = rows_p * size_x;
            displs[p] = offset * size_x;
            offset += rows_p;
        }

    }
    else    //the rank is != 0
    {
        sendcounts = NULL;
        displs = NULL;
    }


    MPI_Scatterv(
        global_values,  //only root uses it
        sendcounts,     //root
        displs,         //root  //the r0 process will insert in global_values+displs[p] the data.
        MPI_FLOAT,
        &values[size_x],         // it starts at row y=1. th y=0 will be already filled in the halo exchange phase, for each iteration, with MPI_Sendrecv
        local_ny * size_x,       //Each process wants to receive that number of elements.
        MPI_FLOAT,
        0,  //root process, in this case rank0
        MPI_COMM_WORLD
    );

    memcpy(prev_values, values, (local_ny + 2) * size_x * sizeof(stencil_t)); //the previous state is the current state when initializing...

    if (rank == 0) {

        free(global_values);    //free, it's no longer needed
        free(sendcounts);
        free(displs);
    }
}

static void stencil_free(void)
{
  free(values);
  free(prev_values);
}




/*each process sends their 1st and last real row, and each one receives its superior halo, and its inferior halo (row).*/
/*two calls. the first to send the superior halo And receive the */
void exchange_halos(void)
{
    int up   = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int down = (rank == nprocs - 1) ? MPI_PROC_NULL : rank + 1;

    MPI_Sendrecv(
        &prev_values[1 * size_x],        // envío fila real superior sendbuf
        size_x,
        MPI_FLOAT,
        up,     //dest, mando mi primera fila real al proceso de arriba, el anterior proceso (que la recibirá como su halo inferior, aha.)
        0,      //sendtag
        &prev_values[(local_ny + 1) * size_x], //recibo halo inferior recvbuf, Soy elprde arriba, gracias. Ya +2 no tendría sentido, se va. reservé +2 sí, pero 1 parriba y otro pa abajo. 
        size_x,
        MPI_FLOAT,
        down,   //source, el source es down, sería el siguiente tio, por qué? porque es ese quien te ha enviado el halo inf loco
        0,      //recvtag
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE
    );

    MPI_Sendrecv(
        &prev_values[local_ny * size_x], // envío fila real inferior
        size_x,
        MPI_FLOAT,
        down,   //envio mi ultima fila al siguiente bobi, que lo tomará como su halo superior claro el siguiente tio claro
        1,
        &prev_values[0 * size_x],         // recibo halo superior
        size_x,
        MPI_FLOAT,
        up,     //recibo mi halo sup de quién, del tio de antes que me envió su última fila real. gracias
        1,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE
    );
}

/** compute the next stencil step, return 1 if computation has converged */
static int stencil_step_mpi(void)
{
  int convergence = 1;

  /* switch buffers */
  stencil_t*tmp = prev_values;
  prev_values = values;
  values = tmp;
    /*despues de este swap:  prev_values → t+1    values      → buffer libre para t+2.
    Cada paso necesita:   leer estado anterior completo    escribir estado nuevo completo

    Por tanto, en cada iteración:    usar el resultado anterior como entrada   producir uno nuevo
    
    prev_values = foto congelada del sistema
    values = lienzo vacío
    exchange_halos() siempre debe trabajar sobre el estado que se va a leer.  Es decir: sobre prev_values. Nunca sobre values. Porque los halos representan datos vecinos del paso t. No del paso t+1 (que aún no existe)
    */

  exchange_halos(); //we update here prev_values, matrix(vector1D) that will be used in the stencil op

  int x, y;
  for(y = 1; y <= local_ny ; y++)     //until local_ny, not size_y global size of the matrix...
    {
      for(x = 1; x < size_x - 1 ; x++)
        {
          values[x + size_x * y] =
            alpha * ( prev_values[x - 1 + size_x * y] +
                      prev_values[x + 1 + size_x * y] +
                      prev_values[x + size_x * (y - 1)] +
                      prev_values[x + size_x * (y + 1)]) +
            (1.0 - 4.0 * alpha) * prev_values[x + size_x * y];
          if(convergence && (fabs(prev_values[x + size_x * y] - values[x + size_x * y]) > epsilon))
            {
              convergence = 0;
            }
        }
    }
  return convergence;
}

void display_matrix(int s, double t_usec)
{
   //in this function we put together the results from each process, so that rank0 can gather and print, and show metrics
    int base = size_y / nprocs;
    int rest = size_y % nprocs;

/* <-- uncomment to see result, begin 
    if (rank == 0) {
        global_values = malloc(size_y * size_x * sizeof(stencil_t));
        sendcounts = malloc(nprocs * sizeof(int));
        displs     = malloc(nprocs * sizeof(int));

        int offset = 0;
        for (int p = 0; p < nprocs; p++) {
            int rows_p = base + (p < rest ? 1 : 0);
            sendcounts[p] = rows_p * size_x;
            displs[p] = offset * size_x;
            offset += rows_p;
        }
    } else{
        sendcounts = NULL;
        displs = NULL;
    }

    MPI_Gatherv(
        &values[size_x],          // envío: filas reales
        local_ny * size_x,
        MPI_FLOAT,
        global_values,            // recepción (solo rank 0)
        sendcounts,
        displs,
        MPI_FLOAT,
        0,
        MPI_COMM_WORLD
    );
 <-- uncomment to not see result, end */

    if (rank == 0) {
        // imprimir tiempo, pasos, gflops
        printf("# steps = %d\n", s);
        printf("# time = %g usecs.\n", t_usec);
        printf("# gflops = %g\n", (6.0 * size_x * size_y * s) / (t_usec * 1000));
       
        //stencil_display_global(global_values);
        //free(global_values);
        //free(sendcounts);
        //free(displs);

    }

}


int main(int argc, char**argv)
{

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


  stencil_init_mpi();
  

  MPI_Barrier(MPI_COMM_WORLD);  // synchronize before measuring time

  if(rank == 0)
    printf("Number of MPI processes: %d\n", nprocs);

  struct timespec t1, t2;
  clock_gettime(CLOCK_MONOTONIC, &t1);

  int s, convergence_global = 0;
  for(s = 0; s < stencil_max_steps; s++)
    {
      int convergence_local = stencil_step_mpi();


      if(s % 100 == 0) {  //for each 100 iterations, we check convergence
        //All make a logic AND logic
        MPI_Allreduce(&convergence_local, &convergence_global, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
      }
      if (convergence_global) break;
    }

  MPI_Barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_MONOTONIC, &t2);
  const double t_usec = (t2.tv_sec - t1.tv_sec) * 1000000.0 + (t2.tv_nsec - t1.tv_nsec) / 1000.0;

 

  display_matrix(s, t_usec);
    


  stencil_free();
  MPI_Finalize();
  return 0;
}
