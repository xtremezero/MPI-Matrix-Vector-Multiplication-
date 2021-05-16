# define N_DIM 4   //  confirmed
#include <assert.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <cstdlib> 
#include <stdio.h>
#include <mpi.h> 

using namespace std;

void read_mat_from_file(const char *s, int n_row, int n_col, double  *in_matrix);
void RowMatrixVectorMultiply(int dim, double *mat, double *vec, double *result);

int main( int argc, char *argv[]){
    int rank, size;                    //MPI RANK,SIZE
    
    MPI_Init (&argc, &argv);                //Initializations
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size); 

    if (N_DIM%size){  //Valid Communicator Size?
        MPI_Finalize();
        return(0);
    }
    
    double matrix_data[N_DIM][N_DIM];  //matrix
    double vector_data[N_DIM];         //vector
    double result[N_DIM] = {0.0};   //final result holder
 
    
    if (rank==0){
        read_mat_from_file("matrix.txt", N_DIM, N_DIM, (double *)matrix_data);   //Populating the Matrix
        read_mat_from_file("vector.txt", N_DIM, 1, vector_data);                 //Populating the Vector
    }
    RowMatrixVectorMultiply(N_DIM, (double *)matrix_data, vector_data,result);
    
    /* Printing the Matrix*/
    if (rank==0){
        printf("Matrix  :\n");
        for (int i=0;i<N_DIM;i++){
            for (int j=0;j<N_DIM;j++)
                printf("%.5f ", matrix_data[i][j]);
            printf("\n");
        }
        printf("Vector :\n");
        for (int i=0;i<N_DIM;i++)
            printf("%.5f ", vector_data[i]);
        printf("\n\n");
        
        printf("Vector :\n");
        for (int i=0;i<N_DIM;i++)
            printf("%.5f ", vector_data[i]);
        printf("\n\n");
        
        printf("Result :\n");
        for (int i=0;i<N_DIM;i++)
            printf("%.5f ", result[i]);
        printf("\n\n");
    }
    
    MPI_Finalize();
    return(0);
}

void RowMatrixVectorMultiply(int dim, double *matrix_data, double *vector_data,double *result){
    int rank,size;
    MPI_Status status;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size); 
    double* localresult = new double[dim / size]{};  //local result
    double matrix [dim][dim];   //local matrix
    double timer=MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    ///////////////////////////////////////////////////////////////Comm Start///////////////////////////////////////////////////////////////////
    //MPI_Scatter(matrix_data, (dim*dim)/size, MPI_DOUBLE, matrix, (dim*dim)/size, MPI_DOUBLE, 0, MPI_COMM_WORLD); //Scatter the Matrix
    if (rank==0){
        for (int i=0;i<dim/size;i++)
            for (int j=0;j<dim;j++)
                matrix[i][j] = matrix_data[j+(dim*i)];
        for (int i=1;i<size;i++){    
            MPI_Send(matrix_data+((dim*dim/size)*(i)),(dim*dim)/size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }
    else{
        MPI_Recv(matrix, dim*dim/size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    }
    /////////////////////////////////////////////////////////////////Comm End-Start/////////////////////////////////////////////////////////////
    //Broadcast
    if (rank == 0) {
        for (int i = 0; i < size; i++) {
            if (i != rank) {
                MPI_Send(vector_data, dim, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }
        }
    } else {
        MPI_Recv(vector_data, dim, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    
    //MY Broadcast - 2 processes only
    /*int limit = sqrt(size);
    for (int i=0;i<limit;i++){
        for (int j=0;j<=i;j++){
            if (rank == size-(size/pow(2,j))){
                MPI_Send(vector_data, dim, MPI_DOUBLE,rank+((size)/pow(2,i+1)), i, MPI_COMM_WORLD);
            }else{
                int receiver = (rank+1)-(size/pow(2,i));
                if (receiver>=0)
                MPI_Recv(vector_data, dim, MPI_DOUBLE, receiver, i, MPI_COMM_WORLD, &status);
            }
        }
    }*/
    ///////////////////////////////////////////////////////////////Comm End///////////////////////////////////////////////////////////////////
    //Calculate the results
    for (int i = 0;i<(dim/size);i++)
        for (int j = 0;j<dim;j++)
            localresult[i]+=vector_data[j]*matrix[i][j];
        
    /*cout <<"Debug INFO:\n";
    for (int i = 0;i<(N_DIM/size);i++)
        for (int j = 0;j<N_DIM;j++)
            cout << "am rank "<< rank<<" and my matrix"<<matrix[i][j]<<endl;
    
    
    for (int j = 0;j<N_DIM;j++)
        cout << "am rank "<< rank<<" and my vector is "<<vector_data[j]<<endl;
    
    
    
    for (int i = 0;i<(N_DIM/size);i++)
        cout << "am rank "<< rank<<" and my result is "<<localresult[i]<<endl;*/
    ///////////////////////////////////////////////////////////////Comm Start///////////////////////////////////////////////////////////////////
    //MPI_Gather(localresult, (dim)/size, MPI_DOUBLE, result, (dim)/size, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Gather the results
    if (rank==0){
        for (int i=0;i<dim/size;i++)
            result[i] = localresult[i];
        for (int i=1;i<size;i++){    
            MPI_Recv(result+((dim/size)*(i)),(dim)/size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    else{
        MPI_Send(localresult, dim/size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    
    }
    ///////////////////////////////////////////////////////////////Comm End///////////////////////////////////////////////////////////////////
    /*DEBUG INFO
    if (rank==0)
    for (int j = 0;j<N_DIM;j++)
        cout << "am g-Root and my result is "<<result[j]<<endl;*/
    timer = MPI_Wtime()-timer;
    if (rank==0)
        cout << "Time Needed for all ops = "<<timer<<endl;
}

void read_mat_from_file(const char *s, int n_row, int n_col, double  *in_matrix){
	std::ifstream fin(s);
	std::string line;
	double data_in;

//open the input stream
	if(!fin)
	{
		cout << "Unable to open " << s << " for reading.\n";
		exit (0);
	}

//	cout << " for file: " << s << "\n" << endl;
	for (int i = 0; i < n_row; i++)
	{
		std::getline(fin, line);
	    std::stringstream stream(line);
		for (int j = 0; j < n_col; j++)
		{
			stream >> data_in;       //now read the whitespace-separated floats
//			cout << " data_in = (" << data_in_real  << "," << data_in_imag << ") \n" << endl;
			*(in_matrix+(i*n_col)+j) = data_in;
		}  //  for (int j = 0; j < n_col; j++)
	}  //	for (int i = 0; i < n_row; i++)
}  // void read_n_complex_pair(string s, int n)

