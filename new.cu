
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
	
	float ac, bb0,sb1,sb0,cb0,k,bb1,ssb1,cb1,ccb1,pom,minus,factorial, factorials,fb,fsb ,k1, pom1,ccd; 
    int angle;  float fix2=57;
    	 bb0=b0[idx];  bb1=b1[idx]; ssb1=sin(bb1); 
     		sb0=sin(bb0); cb0=cos(bb0); cb1=cos(bb1);
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
		  k=b1[i]; k1=a0[i]-a0[idx]
		  pom=k*k; pom1=k1*k1;
		  minus=-1;  
		  sb1=k; ccb1=1;ccd1=1; 
		  fb=3;
		  fsb=4;
		  factorial=2; factorials=6;
		for(int i=0; i<5; i++)
		{ 
			ccb1=ccb1+minus*pom/factorial; 
			sb1=sb1+minus*k*pom/factorials; 
			ccd=ccd+minus*pom1/factorial; 
			factorial=factorial*fb*(fb+1); fb+=2;
			factorials=factorials*fsb*(fsb+1);  fsb+=2; minus=minus*(-1); pom=pom*k*k;  
		}
		 /* sb1=k-k*k*k/6 + k*k*k*k*k/120- k*k*k*k*k*k*k/5040+k*k*k*k*k*k*k*k*k/362880-k*k*k*k*k*k*k*k*k*k*k/39916800;
		  ccb1=1-k*k/2+k*k*k*k/24-k*k*k*k*k*k/720+k*k*k*k*k*k*k*k/40320-k*k*k*k*k*k*k*k*k*k/3628800;*/
          	  ac= acosf(sb0*sb1+ cb0*ccb1*ccd);
		  ac= (ac*fix2/0.25); 
	
		  angle=(int) ac; 
		  atomicAdd(&mn[angle],1);
		
		}
		
	   for(int i=idx+1; i<100000;i++)
	    { 		
		   k=b0[i]; 
		   pom=k*k; 
		   minus=-1;  
		   sb1=k; ccb1=1; 
		   fb=3;
		   fsb=4;
		   factorial=2; factorials=6;
		for(int i=0; i<5; i++)
		{ 
			ccb1=ccb1+minus*pom/factorial; 
			sb1=sb1+minus*k*pom/factorials; 
			factorial=factorial*fb*(fb+1); fb+=2;
			factorials=factorials*fsb*(fsb+1);  fsb+=2; minus=minus*(-1); pom=pom*k*k;  
		}
		 /*  sb1=k-k*k*k/6 + k*k*k*k*k/120- k*k*k*k*k*k*k/5040+k*k*k*k*k*k*k*k*k/362880-k*k*k*k*k*k*k*k*k*k*k/39916800;
		   ccb1=1-k*k/2+k*k*k*k/24-k*k*k*k*k*k/720+k*k*k*k*k*k*k*k/40320-k*k*k*k*k*k*k*k*k*k/3628800;*/
		   ac= acosf(sb0*sb1+ cb0*ccb1*cos((a0[i]-a0[idx])));
	           ac= (ac*fix2/0.25); 
                   angle=(int) ac; 
                   atomicAdd(&r[angle],2);
		   
	     	   k=b1[i];
		   pom=k*k; 
		  minus=-1;  
		  sb1=k; ccb1=1;
		  fb=3;
		  fsb=4;
		  factorial=2; factorials=6;
		for(int i=0; i<5; i++)
		{ 
			ccb1=ccb1+minus*pom/factorial; 
			sb1=sb1+minus*k*pom/factorials; 
			factorial=factorial*fb*(fb+1); fb+=2;
			factorials=factorials*fsb*(fsb+1);  fsb+=2; minus=minus*(-1); pom=pom*k*k;  
		}
		 /*  sb1=k-k*k*k/6 + k*k*k*k*k/120- k*k*k*k*k*k*k/5040+k*k*k*k*k*k*k*k*k/362880-k*k*k*k*k*k*k*k*k*k*k/39916800;
		   ccb1=1-k*k/2+k*k*k*k/24-k*k*k*k*k*k/720+k*k*k*k*k*k*k*k/40320-k*k*k*k*k*k*k*k*k*k/3628800;*/
                   ac= acosf((ssb1*sb1)+ cb1*ccb1*cos((a1[idx]-a1[i])));
                   ac= (ac*fix2/0.25); 
	           angle=(int) ac; 
                   atomicAdd(&s[angle],2);

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
	int xx=782; 
size_t arraybytes = N * sizeof(float);
	size_t arraybytes1 =xx *720 *sizeof(int);
	size_t arraybytes11 =xx *720 *sizeof(float);
	size_t l=720*sizeof(int);
	size_t l1=720*sizeof(float); // ovo proveri
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
	float* final = (float*)malloc(l1); 
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
	int* d_result; cudaMalloc(&d_result, l);
	int* d_result_r; cudaMalloc(&d_result_r, l);
	int* d_result_s; cudaMalloc(&d_result_s, l);
// Copy arrays from host memory to device memory
cudaMemcpy(d_A, h_A, arraybytes, cudaMemcpyHostToDevice);
cudaMemcpy(d_B, h_B, arraybytes, cudaMemcpyHostToDevice);
cudaMemcpy(d_A1, h_A1, arraybytes, cudaMemcpyHostToDevice);
cudaMemcpy(d_B1, h_B1, arraybytes, cudaMemcpyHostToDevice);
// Invoke kernel
	dim3 threadsPerBlock(128);
	dim3 threadsPerBlock1(736);
	dim3 blocksize2(1);
 
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

	result_r[0]=result_r[0]+100000; result_s[0]=result_s[0]+100000; 

	final[0]=(float) ((float)(result_r[0]-2*result[0]+result_s[0]+200000)/(float)(100000+ result_s[0]));	
	for(int i=1;i<720;i++)
		final[i]=(float) ((float)(result_r[i]-2*result[i]+result_s[i])/(float) result_s[i]);
	
	
	end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	
	double brk=0; 
	printf("%f\n", cpu_time_used); 
	
	for(int i=0; i<720; i++)
		{
		 
		printf( "%f ", final[i]);
	}
	
		
	

cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
	cudaFree(h_A); cudaFree(h_B); cudaFree(h_C);cudaFree(d_D);cudaFree(h_D);cudaFree(d_E);cudaFree(h_E);
	cudaFree(d_A1);cudaFree(h_A1);cudaFree(d_B1);cudaFree(h_B1); 
	
}








    
  
