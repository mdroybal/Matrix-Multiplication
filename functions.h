/*
Monte Roybal
Function's Declarations
*/

/*Link Section*/
#ifndef FUNCTIONS_H_ 
#define FUNCTIONS_H_

void open_file(FILE **ptr, const char file[48],const char mode [8]);//Open File Function Declaration
void allocate_mem_2d(double **arr,int n, int m);//2 Dimensional Memory Allocation Function Declaration
void allocate_mem_1d(double **arr,int m);//1 Dimensional Memory Allocation Function Declaration
void print_matrices(double *matr,const char matr_name[50],int rows,int cols,int X);//*Print NxM Matrix Function Declaration
void print_vectors(double *vec,const char vec_name[50],int n);//Print Vector Function Declaration
void matrix_matrix_mult(double *h_c,double *h_a,double *h_b,int n,int m);//Matrix Matrix Multiply Function Declaration
void matrix_vector_mult(double *h_c,double *h_a,double *h_b,int n,int m);//Matrix Vector Multiply Function Declaration
void matrix_transpose(double *matr_trans,double *matr,int n,int m);//Matrix Transpose Function Declaration

#endif