
---COMMENTS----------------------------------------------------

Comments from the first MPI version file are up to date and in english

However, in the other two files, comment are in spanish.



The comments that are in the first MPI file are the ones indicated to understand the code.








---COMPILATION AND EXECUTION-------------------------------------

MPI pure version:
make

mpirun -H miriel002:24,miriel003:24 -np 4 ./stencil
(two nodes, four processes per example)




OPENMP version:
export OMP_NUM_THREADS= 6   to set the number of threads, six for example.
./stencil




OPENMP + MPI version:
export OMP_NUM_THREADS= 6

and then run with mpirun

mpirun -H miriel002:24,miriel003:24 -np 4 ./stencil

