#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include <time.h>
#include "device_functions.h"
#include "cuda.h"
#include<cuda_runtime.h>
#include <time.h>

using namespace std;

__global__ void angles(volatile float *a0, volatile float *b0,   volatile int *hist, volatile int* hist_r, volatile int* hist_s)

{
	int angle;
    int idx = blockIdx.x * blockDim.x + threadIdx.x; 

   
    __shared__ int mn[720], r[720], s[720];
    if(threadIdx.x==0)
    {
        for (int i=0;i<720;i++)
	{ mn[i] = 0; r[i]=0;s[i]=0;} 
    }
    __syncthreads();


    if (idx<10000)
    {
      
        for(int i=0; i<10000; i++)
        	{
            angle= (int)(a0[idx]*b0[i]); 
             atomicAdd(&mn[angle],1);
                }
	for(int i=idx+1; i<10000; i++)
	{
	     angle= (int)(a0[idx]*a0[i]); 
             atomicAdd(&r[angle],1);
	     angle= (int)(b0[idx]*b0[i]); 
             atomicAdd(&s[angle],1);
	
	}
               

            }
        
    

    __syncthreads();

    if(threadIdx.x==0)
    {
        for(int i=0;i<720;i++)
	{ hist[i+(blockIdx.x*720)]=mn[i]; hist_r[i+(blockIdx.x*720)]=r[i]; hist_s[i+(blockIdx.x*720)]=s[i];}
    }

}

int main(int argc, char *argv[])
{
	
 

int N =10000;
size_t arraybytes = N * sizeof(float);
	size_t arraybytes1 = 20*720 *sizeof(int);
	size_t l=720*sizeof(int);
// Allocate input vectors h_A and h_B in host memory
float* h_A = (float*)malloc(arraybytes);
float* h_B = (float*)malloc(arraybytes);
int* h_C = (int*)malloc(arraybytes1); 
int* h_D = (int*)malloc(arraybytes1); 	
int* h_E = (int*)malloc(arraybytes1); 	
	int* result=(int*)malloc(l); 
	int* result_r=(int*)malloc(l); 
	int* result_s=(int*)malloc(l); 
	for(int i=0; i<10000; i++)
	{ h_A[i]=1.0; h_B[i]=1.0;  }
	h_A[0]=5.0; h_B[1] =3.0; 
float* d_A; cudaMalloc(&d_A, arraybytes);
float* d_B; cudaMalloc(&d_B, arraybytes);
int* d_C; cudaMalloc(&d_C, arraybytes1);
	int* d_D; cudaMalloc(&d_D, arraybytes1);
	int* d_E; cudaMalloc(&d_E, arraybytes1);
// Copy arrays from host memory to device memory
cudaMemcpy(d_A, h_A, arraybytes, cudaMemcpyHostToDevice);
cudaMemcpy(d_B, h_B, arraybytes, cudaMemcpyHostToDevice);
// Invoke kernel
	
	clock_t start, end;
    int threadsPerBlock=736;
    int blocksPerGrid=15; 
     double cpu_time_used;
     
     start = clock();
    cudaMemset(d_C,0,arraybytes1);
	cudaMemset(d_D,0,arraybytes1);
	cudaMemset(d_E,0,arraybytes1);
	
 	angles<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C,d_D,d_E);

      cudaMemcpy(h_C, d_C, arraybytes1, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_D, d_D, arraybytes1, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_E, d_E, arraybytes1, cudaMemcpyDeviceToHost);
	
	for(int i=0; i<720*20; i++)
	{	result[i%720]+= h_C[i];result_r[i%720]+=h_D[i];result_s[i%720]+=h_E[i];} 

		
	
	
	end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("%f\n", cpu_time_used); 
		for(int i=0; i<720; i++)
		{printf("%d ", result[i]);   }
	printf("\n druga\n " ) ; 
	for(int i=0; i<720; i++)
		{printf("%d ", result_r[i]);   }
	printf("\n treca\n " ) ; 
	for(int i=0; i<720; i++)
		{printf("%d ", result_s[i]);   }

cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
	cudaFree(h_A); cudaFree(h_B); cudaFree(h_C);cudaFree(d_D);cudaFree(h_D);cudaFree(d_E);cudaFree(h_E);

	
}








    
  
