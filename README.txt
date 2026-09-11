STENCIL OPERATION PARALLELIZATION 



COMPILATION AND EXECUTION FOR THE MPI version -------------------------------------
OpenMPI paralellizes at the level of processes.

For the MPI pure version, write in the terminal:

sudo apt update
sudo apt install openmpi-bin libopenmpi-dev

make

mpirun ./stencil_mpi
(it will create n MPI process, being n the number of physical cores of your machine)




COMPILATION AND EXECUTION FOR THE OPENMP version -------------------------------------
OpenMP paralellizes at the level of threads.

For the OMP pure version, write in the terminal:


export OMP_NUM_THREADS=6 
(set the number of threads to six for example, or what is best for your machine)(put OMP_NUM_THREADS=n_threads, in this case 6, without spaces)

./stencil_mpi_omp



