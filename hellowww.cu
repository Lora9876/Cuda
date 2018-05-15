
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include <time.h>
#include "device_functions.h"
#include "cuda.h"
#include<cuda_runtime.h>
#define SUBMATRIX_SIZE 16384



            

__global__ void VecAdd(float* A, float* B, float* C, int N)
{		float m;
 		int n; 
 		float *addr; 
		int idx = blockDim.x * blockIdx.x + threadIdx.x;
			
	
    __syncthreads();
	
			if (idx<10000)
				for(int i=0; i<10000; i++)
				{
					C[idx]= A[idx]*B[i];
					 
				}
 							
	
 

}
// CPU Host code
int main(int argc, char *argv[])
{
	
 

int N =10000;
size_t arraybytes = N * sizeof(float);
	size_t arraybytes1 = 720*16384 *sizeof(float);
	size_t l=720*sizeof(float);
// Allocate input vectors h_A and h_B in host memory
float* h_A = (float*)malloc(arraybytes);
float* h_B = (float*)malloc(arraybytes);
float* h_C = (float*)malloc(arraybytes1); 
	float* result=(float*)malloc(l); 
	
	for(int i=0; i<10000; i++)
	{ h_A[i]=1; h_B[i]=1;  }
	h_A[0]=5; h_B[1] =3; 
float* d_A; cudaMalloc(&d_A, arraybytes);
float* d_B; cudaMalloc(&d_B, arraybytes);
float* d_C; cudaMalloc(&d_C, arraybytes1);
// Copy arrays from host memory to device memory
cudaMemcpy(d_A, h_A, arraybytes, cudaMemcpyHostToDevice);
cudaMemcpy(d_B, h_B, arraybytes, cudaMemcpyHostToDevice);
// Invoke kernel

	int thr=512;
	int blocksInGrid=32; 
	
VecAdd<<<blocksInGrid, thr>>>(d_A, d_B, d_C, N);
// Copy result from device memory to host memory
// h_C contains the result in host memory
cudaMemcpy(h_C, d_C, arraybytes, cudaMemcpyDeviceToHost);
	
	for(int i=0; i<720*16384; i++)
	{	result[i%720]+= d_C[i]; } 
		/*
		for(int i=0; i<720*16384; i++)
		printf("%f ", result[i]);   */
// Free device memory
cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
	cudaFree(h_A); cudaFree(h_B); cudaFree(h_C);
// Free host memory ...
	
}










/*
__global__ void angles(volatile float *a0, volatile float *b1, volatile float *histi)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  //int idy  = threadIdx.y + blockDim.y * blockIdx.y;
  float m;
  __shared__ unsigned int shared[20];
    		// za prvu petlju ocistis uvek
    			
	
   	

  if(idx<20)
	 
	 histi[idx]=  5;
             

	
}
             
             
  int main()
             {
               
               float *a, *b,*a0, *b0,*tmp, *tmp1;
               a= (float*) malloc(20* sizeof(float));
               b= (float*) malloc(20*sizeof(float));
                 tmp= (float*) malloc(20*sizeof(float));
	  a0= (float*) malloc(20* sizeof(float));
               b0= (float*) malloc(20*sizeof(float));
                 tmp1= (float*) malloc(20*sizeof(float));
                
               for(int i=0; i<20;i++)
               { a[i]= i+1; b[i]=i+2; tmp[i]=0;}
            
               cudaMemcpy(a0, a, 20* sizeof(float), cudaMemcpyHostToDevice );
               cudaMemcpy(b0, b,20* sizeof(float), cudaMemcpyHostToDevice );
                dim3 grid, block;
    
               grid.x = 1024; 
                	

                          block.x = 1; 
                angles<<<block, grid>>>(a0, b0, tmp1);
           	   cudaMemcpy(tmp, tmp1, 20*sizeof(float), cudaMemcpyDeviceToHost);
               
               for(int i=0; i<20;i++)
                 printf("%d ", tmp[i]); 
           //    free(a0);
	 // free(a); free(b); free(b0); free(tmp); free(tmp1); 
             return EXIT_SUCCESS;
             
    }
  */
