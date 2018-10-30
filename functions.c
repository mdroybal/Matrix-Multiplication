/*
Monte Roybal
Function Definitions
*/

/*Link Section*/
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"


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


/*Matrix Matrix Multiply Function*/
void matrix_matrix_mult(double *h_c,double *h_a,double *h_b,int n,int m)
{
    int i,j,k;
    double sum_m = 0.0;
    for (i=0;i<m;i++)
    {
        for (j=0;j<m;j++)
        {
            for (k=0;k<n;k++)
            {
                sum_m += h_a[i * n + k] * h_b[k * m + j];
            }
            h_c[i * m + j] = sum_m;
            sum_m = 0;
        }
    }
}

/*Matrix Vector Multiply Function*/ 
void matrix_vector_mult(double *h_c,double *h_a,double *h_b,int n,int m)
{
    int i,j; 
    double sum_b = 0.0;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            sum_b = sum_b + h_a[i * n + j] * h_b[j];
        }
        h_c[i] = sum_b;
        sum_b = 0.0;
    }
}

/*Matrix Transpose Function*/
void matrix_transpose(double *matr_trans,double *matr,int n,int m)
{
    int i,j;
    for (i=0;i<m;i++)
    {
        for (j=0;j<n;j++)
        {
            matr_trans[i * n + j] = matr[j * m + i];
        }
    }
}
