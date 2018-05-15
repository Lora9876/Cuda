
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include <time.h>
#include<cuda_runtime.h>

__global__ void VecAdd(float* A, float* B, float* C, int N)
{
int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j= threadIdx.y; 
	
		C[i]=C[i]; 

}
// CPU Host code
int main(int argc, char *argv[])
{
	
	 int numBlocks;        // Occupancy in terms of active blocks
    int blockSize = 32;

    // These variables are used to convert occupancy to warps
    int device;
    cudaDeviceProp prop;
    int activeWarps;
    int maxWarps;

    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);
    
    

   // activeWarps = numBlocks * blockSize / prop.warpSize;
   // maxWarps = prop.maxThreadsPerMultiProcessor / prop.warpSize;

   printf("%d\n", prop.warpSize); 
	printf("%d\n", prop.maxThreadsPerMultiProcessor); 
    

int N =1024;
size_t arraybytes = N * sizeof(float);
// Allocate input vectors h_A and h_B in host memory
float* h_A = (float*)malloc(arraybytes);
float* h_B = (float*)malloc(arraybytes);
float* h_C = (float*)malloc(arraybytes); 
	for(int i=0; i<1024; i++)
	{ h_A[i]=i; h_B[i]=i+1;  }
float* d_A; cudaMalloc(&d_A, arraybytes);
float* d_B; cudaMalloc(&d_B, arraybytes);
float* d_C; cudaMalloc(&d_C, arraybytes);
// Copy arrays from host memory to device memory
cudaMemcpy(d_A, h_A, arraybytes, cudaMemcpyHostToDevice);
cudaMemcpy(d_B, h_B, arraybytes, cudaMemcpyHostToDevice);
// Invoke kernel
/*dim3 thr,blocksInGrid;	
// thr.x = 256;
	thr.y=256; 
 blocksInGrid.x = 1;*/
	dim3 thr(32,32), blocksInGrid(1);
	
VecAdd<<<blocksInGrid, thr>>>(d_A, d_B, d_C, N);
// Copy result from device memory to host memory
// h_C contains the result in host memory
cudaMemcpy(h_C, d_C, arraybytes, cudaMemcpyDeviceToHost);
	
	for(int i=0; i<1024; i++)
	{printf("%f  ", h_A[i]); 
		printf("%f\n", h_C[i]); h_C[i]=0;  }
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
