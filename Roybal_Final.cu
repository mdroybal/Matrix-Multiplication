/*
Monte Roybal
CS_577 Parallel and Distributed Programming
5-10-2018
Dr. Gil Gallegos
N x M Matrix Solution Kernels Final
*/

/*Link Section*/
#include <stdio.h>
#include <stdlib.h>

/*Global Constant's Declarations and Definitions*/
#define N 5
#define M 4


/*Matrix Vector Multiply Kernel*/
__global__ void DeviceMatrixVectorMult(double *c,double *a,double *b,int n)//Matrix x Vector Multiplication Kernel
{
        int j;
        double c_sum = 0.0;
        int idx = threadIdx.x;
        
        for (j=0;j<n;j++)
        {
            c_sum += a[idx * n + j] * b[j];
        }   
        c[idx] = c_sum;
}

/*Matrix Matrix Multiply Kernel*/
__global__ void DeviceMatrixMatrixMult(double *c,double *a,double *b,int n,int m)//Matrix x Vector Multiplication Kernel
{
        int j,k;
        double sum_m = 0.0;
        int idx = threadIdx.x;
        
        for (j=0;j<m;j++)
        {
            for (k=0;k<n;k++)
            {    
                sum_m += a[idx * n + k] * b[k * m + j];
            }
            c[idx * m + j] = sum_m;
            sum_m = 0.0;
        }   
}

/*Open File Function*/
void open_file(FILE **ptr, const char file[48],const char mode [8])
{
    *ptr = fopen(file, mode);
    if (file == NULL)
    {
    printf("No file found");
    }
}

/*Two-Dimensional Flat Memory Allocation Function*/
void allocate_mem_2d(double **arr, int n,int m)
{
    *arr = (double*)malloc(n * m * sizeof(double));
}

/*One-Dimensional Memory Allocation Function*/
void allocate_mem_1d(double **arr,int m)
{
    *arr = (double*)malloc(m * sizeof(double));
}

/*Print NxM Matrix Function*/
void print_matrices(double *matr,const char matr_name[50],int rows,int cols,int X)
{
    int i,j;
    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {    
            printf("%s[%d][%d] = %f ",matr_name,i,j,matr[i *X + j]);
        }
        printf("\n");
    }
    printf("\n");
}

/*Print Vector Function*/
void print_vectors(double *vec,const char vec_name[50],int n)
{
    int i;
    for (i=0;i<n;i++)
    {
        printf("%s[%d] = %f\n",vec_name,i,vec[i]);
    }
    printf("\n");
}

/*Matrix Transpose Function*/
void matrix_transpose(double *matr_trans,double *matr,int n,int m)
{
    int i,j;
    /*Transpose A Matrix*/
    for (i=0;i<m;i++)
    {
        for (j=0;j<n;j++)
        {
            matr_trans[i * n + j] = matr[j * m + i];
        }
    }
}


/*Main Function*/
int main(void)
{
    /*Local Variable's Declarations and Definitions*/
    int i,j;
    double temp_vec;
    double *B_vec,*A_matrix,*A_trans,*A_star,*B_star;
    double *dev_a_t,*dev_a_s,*dev_b_s,*dev_a_m,*dev_b_v;
    FILE *A_ptr,*B_ptr,*A_star_ptr,*B_star_ptr;

    /*Open Matrix and Vector Files for Reading*/
    open_file(&A_ptr,"A_matrix_final.dat","r");
    open_file(&B_ptr,"b_vector_final.dat","r");

    /*Open Matrix and Vector Result Files for Writing*/
    open_file(&A_star_ptr,"A_star_final.dat","w");
    open_file(&B_star_ptr,"b_star_final.dat","w");

    /*Host Memory Allocation*/
    allocate_mem_2d(&A_matrix,M,N);
    allocate_mem_2d(&A_trans,N,M);
    allocate_mem_2d(&A_star,M,M);
    allocate_mem_1d(&B_vec,N);
    allocate_mem_1d(&B_star,M);

    /*CUDA Memory Allocation*/
    cudaMalloc((void **) &dev_a_t, M*N*sizeof(double));
    cudaMalloc((void **) &dev_a_m, N*M*sizeof(double));
    cudaMalloc((void **) &dev_a_s, M*M*sizeof(double)); 
    cudaMalloc((void **) &dev_b_v, N*sizeof(double));
    cudaMalloc((void **) &dev_b_s, M*sizeof(double));

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

    /*CUDA Memory Copy from Host fo Device*/
    cudaMemcpy(dev_a_t, A_trans, M*N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_a_m, A_matrix, N*M*sizeof(double), cudaMemcpyHostToDevice);

    /*Matrix Matrix Multiply Kernel Call with 1 Block and N Threads*/
    DeviceMatrixMatrixMult<<<1,N>>>(dev_a_s,dev_a_t,dev_a_m,N,M);
    
    /*CUDA Memory Copy A Star from Device to Host*/
    cudaMemcpy(A_star,dev_a_s,M*M*sizeof(double), cudaMemcpyDeviceToHost);

    /*Print Kernel Calculated A Star to File*/    
    printf("\nFrom Device:\n");
    print_matrices(A_star,"A_star",M,M,M);
    for (i=0;i<M;i++)
    {    
        for (j=0;j<M;j++)
        {
            fprintf(A_star_ptr,"%f ",A_star[i * M + j]);      
        }
        fprintf(A_star_ptr,"\n");
    }

    /*CUDA Memory Copy from Host fo Device*/
    cudaMemcpy(dev_b_v, B_vec, N*sizeof(double), cudaMemcpyHostToDevice);

    /*Matrix Vector Multiply Kernel Call with 1 Block and N Threads*/    
    DeviceMatrixVectorMult<<<1,N>>>(dev_b_s,dev_a_t,dev_b_v,N);

    /*CUDA Memory Copy B Star from Device to Host*/
    cudaMemcpy(B_star,dev_b_s,M*sizeof(double), cudaMemcpyDeviceToHost);

    /*Print Kernel Calculated B Star to File*/ 
    printf("\nFrom Device:\n");
    print_vectors(B_star,"B_star",M);
    for (i=0;i<M;i++)
    {
        fprintf(B_star_ptr, "%f ",B_star[i]);
    }

    /*Free Host Allocated Memory*/
    free(A_matrix);
    free(A_trans);
    free(B_vec);
    free(A_star);
    free(B_star);

    /*Free CUDA Allocated Memory*/
    cudaFree(dev_a_t);
    cudaFree(dev_a_m);
    cudaFree(dev_a_s);
    cudaFree(dev_b_v);
    cudaFree(dev_b_s);
    
    /*Close File Pointers*/
    fclose(A_ptr);
    fclose(B_ptr);
    fclose(A_star_ptr);
    fclose(B_star_ptr);

    /*Return 0 for Successful Execution*/
    return 0;
}