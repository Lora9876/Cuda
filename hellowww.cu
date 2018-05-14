
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include <time.h>
#include<cuda_runtime.h>

__global__ void angles(volatile float *a0, volatile float *b1, volatile float *histi)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy  = threadIdx.y + blockDim.y * blockIdx.y;
  float m;
  __shared__ unsigned int shared[20];
    		// za prvu petlju ocistis uvek
    			if((threadIdx.x==0) && (threadIdx.y==0))
   			 {
       			 for (int i=0;i<20;i++)
         		   shared[i] = 0;
   			 }
	
   	 __syncthreads();

  if(idx<20 && idy<20)
	 m=  a0[idx] + b1[idy];
             
    __syncthreads();
	shared[0]+=idx; 
}
             
             
  int main()
             {
               
               float *a, *b,*a0, *b0,*tmp, *tmp1;
               a= (float*) malloc(20* sizeof(float));
               b= (float*) malloc(20*sizeof(float));
                 tmp= (float*) malloc(20*sizeof(float));
                
               for(int i=0; i<20;i++)
               { a[i]= i+1; b[i]=i+2; tmp[i]=0;}
               cudaMalloc((void **) &a0, 20* sizeof(float));
               cudaMalloc((void **) &b0, 20* sizeof(float));
                cudaMalloc((void **) &tmp1, 20* sizeof(float));
               cudaMemset(a0,0,20* sizeof(float));
                cudaMemset(b0,0,20* sizeof(float));
               cudaMemset(tmp1,0,20* sizeof(float));
               cudaMemcpy(a0, a, 20* sizeof(float), cudaMemcpyHostToDevice );
               cudaMemcpy(b0, b,20* sizeof(float), cudaMemcpyHostToDevice );
                dim3 grid, block;
    
               grid.x = 1024; 
                	 grid.y=1024;

                          block.x = 1; 
                angles<<<grid,block>>>(a0, b0, tmp1);
           	   cudaMemcpy(tmp, tmp1, 20*sizeof(float), cudaMemcpyDeviceToHost);
               
               for(int i=0; i<20;i++)
                 printf("%d ", tmp[i]); 
               
             return EXIT_SUCCESS;
             
    }
  
