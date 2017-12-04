# MPI-Matrix-Vector-Multiplication-
MPI program for cross-multiplying a matrix by a vector in parallel.
<br>
<b>In The folder:</b><br>
MPIVMM.cpp                      #The program code using MPI collective communication functions.<br>
MPIVMMP2P.cpp                   #The program code using MPI point-to-point communication functions. <br>
matrix.txt                      #4x4 input matrix.<br>
matrix_16.txt                   #16x16 input matrix.<br>
vector.txt                      #4 input vector.<br>
vector_16.txt                   #16 input vector.<br>
vector_matrix_result.txt        #4x4 output matrix if the default files were used.<br>
vector_matrix_result_16.txt     #16x16 output matrix if the default files were used.<br>
<br>
<b>#How To Use:</b><br>
1-Start by changing the N_DIM definition in the code to the wanted matrix dimensions.<br>
2-replace the matrix file and change the name of the file in "the read_mat_from_file()" function.<br>
3-compile the program using MPIC++ "C++ file".<br>
4-use mpiexec -n [Multiple of matrix dimensions] executable_file_name.<br>
