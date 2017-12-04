# MPI-Matrix-Vector-Multiplication-
MPI program for cross-multiplying a matrix by a vector in parallel.

In The folder :
MPIVMM.cpp                      #The program code using MPI collective communication functions.
MPIVMMP2P.cpp                   #The program code using MPI point-to-point communication functions. 
matrix.txt                      #4x4 input matrix.
matrix_16.txt                   #16x16 input matrix.
vector.txt                      #4 input vector.
vector_16.txt                   #16 input vector.
vector_matrix_result.txt        #4x4 output matrix if the default files were used.
vector_matrix_result_16.txt     #16x16 output matrix if the default files were used.

#How To Use:
1-Start by changing the N_DIM definition in the code to the wanted matrix dimensions.
2-replace the matrix file and change the name of the file in "the read_mat_from_file()" function.
3-compile the program using MPIC++ "C++ file".
4-use mpiexec -n [Multiple of matrix dimensions] executable_file_name.
