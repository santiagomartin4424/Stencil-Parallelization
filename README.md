# ⚡ Stencil Operation Parallelization

This project presents the parallel implementation of a 2D Stencil operation using **MPI** (Process-level) and **OpenMP** (Thread-level), with its respective benchmarking study.


## Setup

Once cloned the repository, install the required dependencies:

```
sudo apt update
sudo apt install openmpi-bin libopenmpi-dev
```


## How to execute the MPI Version (Process-Level Parallelism)

OpenMPI parallelizes the execution at the level of processes. 

### 1. Compilation

```
cd stencil_mpi
make
```

### 2. Execution.
Executing this, n MPI processes will be created, being n the number of physical cores of your machine.

```
mpirun ./stencil_mpi
```



## How to execute the OpenMP Version (Thread-Level Parallelism)

OpenMP parallelizes the execution at the thread level.

### 1. Set Environment Variables.
Set the target number of threads. Ensure there are no spaces around =. For example, you could use 6 threads:

```
cd stencil_omp
export OMP_NUM_THREADS=6 
```

### 2. Execution. Run the OpenMP executable.

```
./stencil_mpi_omp
```

