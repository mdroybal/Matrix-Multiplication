/*
Monte Roybal
CS_577 Parallel and Distributed Programming
5-10-2018
Dr. Gil Gallegos
N x M Matrix Final
*/

/*Link Section*/
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"

/*Global Constant's Declarations and Definitions*/
#define N 5
#define M 4


/*Main Function*/
int main(void)
{
    /*Local Variable's Declarations and Definitions*/
    int i,j,k;
    double temp_vec;
    double *B_vec,*A_matrix,*A_trans,*A_star,*B_star;
    FILE *A_ptr,*B_ptr,*A_star_ptr,*B_star_ptr;
    
    /*Open Matrix and Vector Files for Reading*/
    open_file(&A_ptr,"A_matrix_final.dat","r");
    open_file(&B_ptr,"b_vector_final.dat","r");

    /*Open Matrix and Vector Result Files for Writing*/
    open_file(&A_star_ptr,"A_star_final.dat","w");
    open_file(&B_star_ptr,"b_star_final.dat","w");

    /*Two-Dimensional Flat Memory Allocation*/
    allocate_mem_2d(&A_matrix,N,M);
    allocate_mem_2d(&A_trans,M,N);
    allocate_mem_2d(&A_star,M,M);

    /*One-Dimensional Memory Allocation*/
    allocate_mem_1d(&B_vec,N);
    allocate_mem_1d(&B_star,M);

    /*Scan Data From File to B Vector*/
    for (i=0;i<N;i++)
    {
        fscanf(B_ptr,"%lf",&temp_vec);
        B_vec[i] = temp_vec;
    }
    print_vectors(B_vec,"B_vec",N);

    /*Scan Data From File to A Matrix*/
    for (i=0;i<N;i++)
    {
        for (j=0;j<M;j++)
        {
            fscanf(A_ptr,"%lf",&temp_vec);
            A_matrix[i* M +j] = temp_vec;
        }   
    }
    print_matrices(A_matrix,"A_matrix",N,M,M);

    /*Transpose A Matrix*/
    matrix_transpose(A_trans,A_matrix,N,M);
    print_matrices(A_trans,"A_trans",M,N,N);
    
    /*Calculate and Print A Star to File*/
    matrix_matrix_mult(A_star,A_trans,A_matrix,N,M);
    print_matrices(A_star,"A_star",M,M,M);  
    for (i=0;i<M;i++)
    {    
        for (j=0;j<M;j++)
        {
            fprintf(A_star_ptr,"%f ",A_star[i * M + j]);      
        }
        fprintf(A_star_ptr,"\n");
    }

    /*Calculate and Print B Star to File*/
    matrix_vector_mult(B_star,A_trans,B_vec,N,M);
    print_vectors(B_star,"B_star",M);
    for (i=0;i<M;i++)
    {
        fprintf(B_star_ptr, "%f ",B_star[i]);
    }

    /*Free Allocated Memory*/
    free(A_matrix);
    free(A_trans);
    free(A_star);
    free(B_vec);
    free(B_star);
    
    /*Close File Pointers*/
    fclose(A_ptr);
    fclose(B_ptr);
    fclose(A_star_ptr);
    fclose(B_star_ptr);

    /*Return 0 for Successful Execution*/
    return 0;
}