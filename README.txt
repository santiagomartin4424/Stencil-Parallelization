# ⚡ Stencil Operation Parallelization

In this project presents the parallel performance optimization of a 2D Stencil computation using **MPI** (Process-level) and **OpenMP** (Thread-level), with its respective benchmarking study.


## Prerequisites & Setup

Before building the project, install the required dependencies:

```
sudo apt update
sudo apt install openmpi-bin libopenmpi-dev
```


## 🌐 How to execute the MPI Version (Process-Level Parallelism)

OpenMPI parallelizes the execution at the level of processes. 

### 1. Compilation

```bash
cd stencil_mpi
make
```

### 2. Execution. (n MPI processes will be created, being n the number of physical cores of your machine)

```
mpirun ./stencil_mpi
```



## 🧵 How to execute the OpenMP Version (Thread-Level Parallelism)

OpenMP parallelizes execution at the thread level.

### 1. Set Environment Variables. (set the target number of threads. Ensure there are no spaces around =. For example, you can use 6 threads)

```
cd stencil_omp
export OMP_NUM_THREADS=6 
```

### 2. Execution. Run the OpenMP executable.

```
./stencil_mpi_omp
```


