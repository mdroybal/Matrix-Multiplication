/*
Monte Roybal
CS_577 Parallel and Distributed Programming
5-2-2018
Dr. Gil Gallegos
Jacobi Kernel Solution with 4x4 and 9x9 A Matrices
*/

/*Link Section*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 4

/*Jacobi Method Algorithm Kernel*/
__global__ void jacobi_kernel(double *a,double *b,double *x_o,double *x_n,int n)
{       
	int j;
	double sum_x = 0.0;
	int idx = threadIdx.x;      
	for (j=0;j<n;j++)
	{
		if(idx!=j)
		{
	        sum_x += a[idx * n + j] * x_o[j];
		}
	}
	x_n[idx] = (b[idx]-sum_x)/a[idx * n + idx];
}

/*Assignment of X new Values to X old Kernel*/
__global__ void x_old_assignment(double *x_o,double *x_n)
{
	int idx = threadIdx.x;
	x_o[idx] = x_n[idx];
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
void allocate_mem_2d(double **arr, int n)
{
	*arr = (double*)malloc(n * n * sizeof(double));
}

/*One-Dimensional Memory Allocation Function*/
void allocate_mem_1d(double **arr,int m)
{
	*arr = (double*)malloc(m * sizeof(double));
}

/*Main Function*/
int main(void)
{
	/*Local Variable's Declaration and Definition*/
	int i,j,k,epsilon_size;
	int count = 0;
	double epsilon[5] = {0.1,0.01,0.001,0.0001,0.00001};
	double sum_error,temp_vec;
	double error_max = 999;
	double *A_matrix;
	double *B_vec,*x_new,*x_old,*x_initial;
	double *dev_a, *dev_b, *dev_xo, *dev_xn;//, *dev_xoa, *dev_xna;
	FILE *A_ptr,*B_ptr;
    
    /*Calculate Epsilon Size*/
	epsilon_size = sizeof(epsilon)/sizeof(double);

    /*Open 4x4 A and B Matrix/Vector Data*/
	open_file(&A_ptr,"A_star_final.dat","r");
    open_file(&B_ptr,"b_star_final.dat","r");		
	
    /*Host Memory Allocation*/
    allocate_mem_2d(&A_matrix,N);    
    allocate_mem_1d(&B_vec,N);
    allocate_mem_1d(&x_old,N);
    allocate_mem_1d(&x_new,N);
    allocate_mem_1d(&x_initial,N);

    /*CUDA Memory Allocation*/
	cudaMalloc((void **) &dev_a, N*N*sizeof(double)); 
	cudaMalloc((void **) &dev_b, N*sizeof(double));
	cudaMalloc((void **) &dev_xo, N*sizeof(double));
	cudaMalloc((void **) &dev_xn, N*sizeof(double));

	/*Scan Data from File to B Vector and Initialize X old*/
	for (i=0;i<N;i++)
	{
		fscanf(B_ptr,"%lf",&temp_vec);
        B_vec[i] = temp_vec;
        x_old[i] = x_initial[i];
    }

	/*Scan Data from File Pointer to A Matrix and Initialize X old*/
	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			fscanf(A_ptr,"%lf",&temp_vec);
            A_matrix[i*N+j] = temp_vec;
		}	
	}

    /*CUDA Memory Copy from Host fo Device*/
    cudaMemcpy(dev_a, A_matrix, N*N*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_b, B_vec, N*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_xo, x_old, N*sizeof(double), cudaMemcpyHostToDevice);

    /*Loop Over Each Epsilon Value for Algorithm Execution*/
    for(k=0;k<epsilon_size;k++)
    {
	    /*Loop Until Error Max is Less Than or Equal to Epsilon*/
	    while (error_max > epsilon[k])
	    {
	    	/*Jacobi Method Algorithm Kernel Call with 1 Block and N Threads */
	        jacobi_kernel<<<1,N>>>(dev_a,dev_b,dev_xo,dev_xn,N);        
	        /*CUDA Memory Copy X new from Device to Host*/
            cudaMemcpy(x_new,dev_xn,N*sizeof(double), cudaMemcpyDeviceToHost);
	        
	        sum_error = 0.0; 
	        /*Calculate a Summation of the Average Errors*/
	        for (i=0;i<N;i++)
	        {
	        	sum_error = sum_error + ((x_new[i]-x_old[i])*(x_new[i]-x_old[i]));       	
	        }
	        
	        /*CUDA Memory Copy X new from Host to Device*/
	        cudaMemcpy(dev_xn, x_new, N*sizeof(double), cudaMemcpyHostToDevice);
	    	/*X new values to X old Assignment Kernel Call with 1 Block and N Threads */
	        x_old_assignment<<<1,N>>>(dev_xo,dev_xn);
			/*CUDA Memory Copy X old from Device to Host*/
			cudaMemcpy(x_old, dev_xo, N*sizeof(double), cudaMemcpyDeviceToHost);     
	        /*CUDA Memory Copy of X old for Next Iteration from Host to Device*/
         	cudaMemcpy(dev_xo, x_old, N*sizeof(double), cudaMemcpyHostToDevice);
	        /*Calculate the Square Root of Error Summations for Stopping Criterion*/
	        error_max = sqrt(sum_error);
	        /*Increment Count to Keep Track of Algorithm Iterations*/
	        count += 1;
	    }  

        /*Print Out X new Values from Jacobi Algorithm Calculations*/ 
		for (i=0;i<N;i++)
		{
	        printf("x_new[%d] = %f \n",i,x_new[i]);
	    }
	    printf("Converged in %d Iterations with %lf Epsilon\n\n",count,epsilon[k]);	    
    }  

    /*Free One-Dimensional CUDA Allocated Memory*/
	cudaFree(dev_a);
	cudaFree(dev_b);
	cudaFree(dev_xo); 
	cudaFree(dev_xn);  
    
    /*Free One-Dimensional Host Allocated Memory*/
    free(A_matrix);
    free(B_vec);
    free(x_initial);
    free(x_old);
    free(x_new);
    
    /*Close File Pointers*/
    fclose(A_ptr);
    fclose(B_ptr);

    /*Return 0 for Successful Execution*/
    return 0;
}
