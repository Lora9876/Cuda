
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

void getDeviceDiagnostics(int tot_gals, int n_coords){

 ////////////////////////////////////////////////////////////////////////////
    // Now get the info from the device.
    ////////////////////////////////////////////////////////////////////////////
   
        printf("\n------ CUDA device diagnostics ------\n\n");

        
        int nx = SUBMATRIX_SIZE;
        int ncalc = nx * nx;
        int gpu_mem_needed = int(tot_gals * sizeof(float)) * n_coords; // need to allocate ra, dec.
        printf("Requirements: %d calculations and %d bytes memory on the GPU \n\n", ncalc, gpu_mem_needed);

        int deviceCount = 0;
        cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
        if (error_id != cudaSuccess) {
            printf( "cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id) );
        }
        // This function call returns 0 if there are no CUDA capable devices.
        if (deviceCount == 0)
            printf("There is no device supporting CUDA\n");
        else
            printf("Found %d CUDA Capable device(s)\n", deviceCount);


        int dev=0;
        for (dev = 0; dev < deviceCount; ++dev) {
            cudaDeviceProp deviceProp;
            cudaGetDeviceProperties(&deviceProp, dev);
            printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

            printf("  Total amount of global memory:                 %.0f MBytes (%llu bytes)\n",
                    (float)deviceProp.totalGlobalMem/1048576.0f, (unsigned long long) deviceProp.totalGlobalMem);


            printf("  Warp size:                                     %d\n", deviceProp.warpSize);
            printf("  Maximum number of threads per block:           %d\n", deviceProp.maxThreadsPerBlock);
            printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
                    deviceProp.maxThreadsDim[0],
                    deviceProp.maxThreadsDim[1],
                    deviceProp.maxThreadsDim[2]);
            printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
                    deviceProp.maxGridSize[0],
                    deviceProp.maxGridSize[1],
                    deviceProp.maxGridSize[2]);

            // does this device have enough capcacity for the calculation?
            printf("\n*************\n");

            // check memory
            if((unsigned long long) deviceProp.totalGlobalMem < gpu_mem_needed) printf(" FAILURE: Not eneough memeory on device for this calculation! \n");
            else
            {
                printf("Hurrah! This device has enough memory to perform this calculation\n");

                // check # threads

                int threadsPerBlock = deviceProp.maxThreadsPerBlock; // maximal efficiency exists if we use max # threads per block.
                int blocksPerGrid = int(ceil(ncalc / threadsPerBlock)); // need nx*nx threads total
                if(deviceProp.maxThreadsDim[0] >blocksPerGrid) printf("FAILURE: Not enough threads on the device to do this calculation!\n");
                else
                {
                    printf("Hurrah! This device supports enough threads to do this calculation\n");
                    // how many kernels can we run at once on this machine?
                    int n_mem = floor(deviceProp.totalGlobalMem / float(gpu_mem_needed));
                    int n_threads = floor(threadsPerBlock * deviceProp.maxThreadsDim[0]*deviceProp.maxThreadsDim[1] / float(ncalc) ); // max # threads possible?

                    printf("%d %d  \n",  n_threads, deviceProp.maxThreadsDim[0]);

                    int max_kernels = 0;
                    n_mem<n_threads ? max_kernels = n_mem : max_kernels = n_threads;

                    printf(" you can run %d kernels at a time on this device without overloading the resources \n", max_kernels);
                }
            }

        }

        printf("\n------ End CUDA device diagnostics ------\n\n");
    }

__global__ void VecAdd(float* A, float* B, float* C, int N)
{		float m;
 		int n; 
 		float *addr; 
		int idx = blockDim.x * blockIdx.x + threadIdx.x;
			__shared__ float sab[720]; 
			
 			if(threadIdx.x==0)
			{	
					for(int i=0; i<720; i++) sab[i]=0; 
			}	
	
    __syncthreads();
	
			if (idx<10000)
				for(int i=0; i<10000; i++)
				{
					m= A[idx]*B[i];
					n= int(m);
					sab[n]=sab[n]+1; 
				}
 							
	 __syncthreads();
 		if(threadIdx.x==0)
   	 {
        for(int i=0;i<720;i++)
            C[i+(blockIdx.x*720)]=sab[i];
    }
	
	
	

}
// CPU Host code
int main(int argc, char *argv[])
{
	
 getDeviceDiagnostics(20000,2); 
    

/*int N =10000;
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
/*dim3 thr,blocksInGrid;	
// thr.x = 256;
	thr.y=256; 
 blocksInGrid.x = 1;
	//dim3 thr(1024), blocksInGrid(100);
	int thr=512;
	int blocksInGrid=32; 
	
VecAdd<<<blocksInGrid, thr>>>(d_A, d_B, d_C, N);
// Copy result from device memory to host memory
// h_C contains the result in host memory
cudaMemcpy(h_C, d_C, arraybytes, cudaMemcpyDeviceToHost);
	
	for(int i=0; i<720*8192; i++)
	{	result[i%720]+= d_C[i]; } 
		
		for(int i=0; i<720*8192; i++)
		printf("%f ", result[i]);   
// Free device memory
cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
	cudaFree(h_A); cudaFree(h_B); cudaFree(h_C);*/
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
