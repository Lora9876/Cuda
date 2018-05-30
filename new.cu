
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include <vector>
#include "vector_types.h"
#include<unistd.h>
#include <time.h>
#include "device_functions.h"
#include "cuda.h"
#include<cuda_runtime.h>
#include <time.h>
#define fix1 3.14/(60*180)


using namespace std;

__global__ void angles(volatile float *a0, volatile float *b0, volatile float *a1, volatile float*b1, volatile int *hist, volatile int* hist_r, volatile int* hist_s)

{
	int idx= blockIdx.x * blockDim.x + threadIdx.x; 
	
	float ac, bb0,aa0,sb0,cb0; 
    int angle;  float fix2=57;
    	bb0=b0[idx]; 
     sb0=sin(bb0); cb0=cos(bb0);
    __shared__ int mn[720], r[720], s[720];
   if(threadIdx.x==0 )
    {
        for (int i=0;i<720;i++)
	{ mn[i] = 0; r[i]=0;s[i]=0;} 
    }
    __syncthreads();


   if ( idx<100000)
    {
      
        for(int i=0; i<100000; i++)
        	{
		   
           ac= acosf(sb0*sin(b1[i])+ cb0*cos(b1[i])*cos((a1[i]-a0[idx])));
		ac= (ac*fix2/0.25); 
	
		angle=(int) ac; 
		  atomicAdd(&mn[angle],1);
		
		}
		
	   for(int i=idx+1; i<100000;i++)
	    { ac= acosf(sb0*sin(b0[i])+ cb0*cos(b0[i])*cos((a0[i]-a0[idx])));
	    ac= (ac*fix2/0.25); 
            angle=(int) ac; 
            atomicAdd(&r[angle],1);
	     
          ac= acosf((sin(b1[idx])*sin(b1[i]))+ cos(b1[idx])*cos(b1[i])*cos((a1[idx]-a1[i])));
            ac= (ac*fix2/0.25); 
	    angle=(int) ac; 
            atomicAdd(&s[angle],1);

                }
	
           }
	

    __syncthreads();

      if(threadIdx.x==0)
    {
        for(int i=0;i<720;i++)
	{  
		
		hist[i+(blockIdx.x*720)]=mn[i]; hist_r[i+(blockIdx.x*720)]=r[i]; hist_s[i+(blockIdx.x*720)]=s[i];}
    }

}

int main(int argc, char *argv[])
{
FILE *real_g; FILE *synthetic_g;
int galaxies_r, galaxies_s; 
clock_t start, end;
start = clock();
	real_g = fopen("data_100k_arcmin.txt","r");
    	synthetic_g = fopen("flat_100k_arcmin.txt","r");	
	fscanf(real_g, "%d", &galaxies_r);
	fscanf(synthetic_g,  "%d", &galaxies_s);
	

int N =100000;
	int xx=137; 
size_t arraybytes = N * sizeof(float);
	size_t arraybytes1 =xx *720 *sizeof(int);
	size_t l=720*sizeof(int);
	size_t l1=720*sizeof(float);
// Allocate input vectors h_A and h_B in host memory
float* h_A = (float*)malloc(arraybytes);
float* h_B = (float*)malloc(arraybytes);
float* h_A1 = (float*)malloc(arraybytes);
float* h_B1 = (float*)malloc(arraybytes);
int* h_C = (int*)malloc(arraybytes1); 
int* h_D = (int*)malloc(arraybytes1); 	
int* h_E = (int*)malloc(arraybytes1); 	
	int* result=(int*)malloc(l); 
	int* result_r=(int*)malloc(l); 
	int* result_s=(int*)malloc(l); 
	float* final=(float*)malloc(l1); 
	for(int i=0; i<galaxies_r; i++)
    {
       
        fscanf(real_g, "%e %e", &h_A[i], &h_B[i]);
       fscanf(synthetic_g, "%e %e", &h_A1[i], &h_B1[i]);
	h_A[i]=h_A[i]*fix1; h_A1[i]=h_A1[i]*fix1; h_B[i]=h_B[i]*fix1; h_B1[i]=h_B1[i]*fix1; }
    fclose(real_g);
	 fclose(synthetic_g);	
	
float* d_A; cudaMalloc(&d_A, arraybytes);
float* d_B; cudaMalloc(&d_B, arraybytes);
float* d_A1; cudaMalloc(&d_A1, arraybytes);
float* d_B1; cudaMalloc(&d_B1, arraybytes);
int* d_C; cudaMalloc(&d_C, arraybytes1);
	int* d_D; cudaMalloc(&d_D, arraybytes1);
	int* d_E; cudaMalloc(&d_E, arraybytes1);
// Copy arrays from host memory to device memory
cudaMemcpy(d_A, h_A, arraybytes, cudaMemcpyHostToDevice);
cudaMemcpy(d_B, h_B, arraybytes, cudaMemcpyHostToDevice);
cudaMemcpy(d_A1, h_A1, arraybytes, cudaMemcpyHostToDevice);
cudaMemcpy(d_B1, h_B1, arraybytes, cudaMemcpyHostToDevice);
// Invoke kernel
	dim3 threadsPerBlock(736);
	
 
    dim3 blocksPerGrid(xx); 
     double cpu_time_used;
     
    
    cudaMemset(d_C,0,arraybytes1);
	cudaMemset(d_D,0,arraybytes1);
	cudaMemset(d_E,0,arraybytes1);
	
		angles<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B,d_A1, d_B1, d_C,d_D,d_E);

      cudaMemcpy(h_C, d_C, arraybytes1, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_D, d_D, arraybytes1, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_E, d_E, arraybytes1, cudaMemcpyDeviceToHost);
	
	for(int i=0; i<720*xx; i++)
	{	result[i%720]+= h_C[i];result_r[i%720]+=h_D[i];result_s[i%720]+=h_E[i];} 

		
	for(int i=0;i<720;i++)
		final[i]=(float) ((float)(result_r[i]-2*result[i]+result_s[i])/(float) result_s[i]); 
	
	end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	
	double brk=0; 
	printf("%f\n", cpu_time_used); 
	
	for(int i=0; i<720; i++)
		{
		 
		printf("%f ", final[i]);
	}
	
		
	

cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
	cudaFree(h_A); cudaFree(h_B); cudaFree(h_C);cudaFree(d_D);cudaFree(h_D);cudaFree(d_E);cudaFree(h_E);
	cudaFree(d_A1);cudaFree(h_A1);cudaFree(d_B1);cudaFree(h_B1);

	
}








    
  
