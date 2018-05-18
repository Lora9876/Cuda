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

__global__ void angles(volatile float *a0, volatile float *b0,   volatile int *hist)

{
    int idx = blockIdx.x * blockDim.x + threadIdx.x; // This should range to SUBMATRIX_SIZE

   
    __shared__ int mn[720];
    if(threadIdx.x==0)
    {
        for (int i=0;i<720;i++)
            mn[i] = 0;
    }
    __syncthreads();


    if (idx<10000)
    {
      
        for(int i=0; i<10000; i++)
        {
            angle= (int)(a0[idx]*b0[i]); 
            
                }

                atomicAdd(&mn[angle],1);

            }
        }
    

    __syncthreads();

    if(threadIdx.x==0)
    {
        for(int i=0;i<720;i++)
            hist[i+(blockIdx.x*720)]=mn[i];
    }

}

int main(int argc, char *argv[])
{
	
 

int N =10000;
	int angle; 
size_t arraybytes = N * sizeof(float);
	size_t arraybytes1 = 20*720 *sizeof(int);
	size_t l=720*sizeof(int);
// Allocate input vectors h_A and h_B in host memory
float* h_A = (float*)malloc(arraybytes);
float* h_B = (float*)malloc(arraybytes);
int* h_C = (int*)malloc(arraybytes1); 
	
	int* result=(int*)malloc(l); 
	
	for(int i=0; i<10000; i++)
	{ h_A[i]=1.0; h_B[i]=1.0;  }
	h_A[0]=5.0; h_B[1] =3.0; 
float* d_A; cudaMalloc(&d_A, arraybytes);
float* d_B; cudaMalloc(&d_B, arraybytes);
int* d_C; cudaMalloc(&d_C, arraybytes1);
	
// Copy arrays from host memory to device memory
cudaMemcpy(d_A, h_A, arraybytes, cudaMemcpyHostToDevice);
cudaMemcpy(d_B, h_B, arraybytes, cudaMemcpyHostToDevice);
// Invoke kernel
int NN=100; 
	
	clock_t start, end;
    int threadsPerBlock=512;
    int blocksPerGrid=20; 
	/* dim3 threadsPerBlock(128, 128);
    dim3 blocksPerGrid(1, 1);
        /*if (NN*NN > 512){
            threadsPerBlock.x = 512;
            threadsPerBlock.y = 512;
            blocksPerGrid.x = ceil(double(NN)/double(threadsPerBlock.x));
            blocksPerGrid.y = ceil(double(NN)/double(threadsPerBlock.y));
        }*/
     double cpu_time_used;
     
     start = clock();
    cudaMemset(d_C,0,arraybytes1);
	
 	angles<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C);

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








    
  
