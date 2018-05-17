
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


            

__global__ void VecAdd(float* A, float* B, int* C,int N,int sum)
{		float m;
 		int n; 
 		
		int idx = blockDim.x * blockIdx.x + threadIdx.x;
		int idy = blockIdx.y*blockDim.y+threadIdx.y;	
 		
 
 			if(idx<sum && idy<sum)
			
			{	for(int i=0;i<sum; i++)
			{		m=A[sum*idx+i]*B[i*sum+idy];
					n=int(m); 
			 		C[idx*10000+ idy*sum+i]=n;
					//mn[n]++;}}
			}}
 			
 						
 //__syncthreads();
 
 

}
// CPU Host code
int main(int argc, char *argv[])
{
	
 

int N =10000;
	int angle; 
size_t arraybytes = N * sizeof(float);
	size_t arraybytes1 = 720*16384 *sizeof(int);
	size_t l=720*sizeof(int);
// Allocate input vectors h_A and h_B in host memory
float* h_A = (float*)malloc(arraybytes);
float* h_B = (float*)malloc(arraybytes);
int* h_C = (int*)malloc(N*arraybytes); 
	
	int* result=(int*)malloc(l); 
	
	for(int i=0; i<10000; i++)
	{ h_A[i]=1; h_B[i]=1;  }
	h_A[0]=5; h_B[1] =3; 
float* d_A; cudaMalloc(&d_A, arraybytes);
float* d_B; cudaMalloc(&d_B, arraybytes);
int* d_C; cudaMalloc(&d_C, N*arraybytes);
	
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
    cudaMemset(d_C,0,100*arraybytes);
	
 	VecAdd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, N,NN);

cudaMemcpy(h_C, d_C, arraybytes, cudaMemcpyDeviceToHost);
	
	result[0] = h_C[0] ; 
	result[1]= h_C[1]; 
	result[2]=h_C[3]; 
/*	for(int i=0; i<N*N; i++)
	{	result[0]= h_C[i]; //angle= h_C[i]; result[angle]++; } */

		
	
	
	end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("%f\n", cpu_time_used); 
		for(int i=0; i<3; i++)
			//if(result[i]>0)
		printf("%d ", result[i]);   
// Free device memory
cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
	cudaFree(h_A); cudaFree(h_B); cudaFree(h_C);
// Free host memory ...
	
}








