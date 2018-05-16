
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
#include <time.h>


            

__global__ void VecAdd(float* A, float* B, int* C, int* D,int N,int sum)
{		float m;
 		int n; 
 		
		int idx = blockDim.x * blockIdx.x + threadIdx.x;
		int idy = blockIdx.y*blockDim.y+threadIdx.y;	
		__shared__ int mn[720];
 			if(threadIdx.x==0)
				for(int i=0; i<720; i++)
					mn[i]=0; 
    __syncthreads();
	
 
 			if(idx<sum && idy<sum)
			
			{	for(int i=0;i<sum; i++)
			{		m=A[sum*idx+i]*B[i*sum+idy];
					n=int(m); 
			 	D[sum*idx+i]=n;
					//mn[n]++;}}
			}}
 			
 							
	__syncthreads();

    if(threadIdx.x==0 && threadIdx.y=0)
    {
        for(int i=0;i<10000;i++)
	{
          n=D[i];
		C[n]++; 
	}
    }
 

}
// CPU Host code
int main(int argc, char *argv[])
{
	
 

int N =10000;
size_t arraybytes = N * sizeof(float);
	size_t arraybytes1 = 720*16384 *sizeof(int);
	size_t l=100*sizeof(int);
// Allocate input vectors h_A and h_B in host memory
float* h_A = (float*)malloc(arraybytes);
float* h_B = (float*)malloc(arraybytes);
int* h_C = (int*)malloc(arraybytes); 
	int* h_D = (int*)malloc(arraybytes); 
	int* result=(int*)malloc(l); 
	
	for(int i=0; i<10000; i++)
	{ h_A[i]=1; h_B[i]=1;  }
	h_A[0]=5; h_B[1] =3; 
float* d_A; cudaMalloc(&d_A, arraybytes);
float* d_B; cudaMalloc(&d_B, arraybytes);
int* d_C; cudaMalloc(&d_C, arraybytes);
	int* d_D; cudaMalloc(&d_D, arraybytes);
// Copy arrays from host memory to device memory
cudaMemcpy(d_A, h_A, arraybytes, cudaMemcpyHostToDevice);
cudaMemcpy(d_B, h_B, arraybytes, cudaMemcpyHostToDevice);
// Invoke kernel
int NN=50; 
	
	clock_t start, end;
	 dim3 threadsPerBlock(NN, NN);
    dim3 blocksPerGrid(1, 1);
        if (NN*NN > 512){
            threadsPerBlock.x = 512;
            threadsPerBlock.y = 512;
            blocksPerGrid.x = ceil(double(NN)/double(threadsPerBlock.x));
            blocksPerGrid.y = ceil(double(NN)/double(threadsPerBlock.y));
        }
     double cpu_time_used;
     
     start = clock();
    cudaMemset(d_C,0,arraybytes);
	 cudaMemset(d_D,0,arraybytes);
 	VecAdd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C,d_D, N,NN);

cudaMemcpy(h_C, d_C, arraybytes, cudaMemcpyDeviceToHost);
	
	for(int i=0; i<10000; i++)
	{	result[i%100]+= h_C[i]; } 

		
	
	
	end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("%f\n", cpu_time_used); 
		for(int i=0; i<720; i++)
			//if(result[i]>0)
		printf("%d ", result[i]);   
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
