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



using namespace std;

__global__ void angles(volatile float *a0, volatile float *b0, volatile float *a1, volatile float*b1, volatile int *hist, volatile int* hist_r, volatile int* hist_s)

{
	int idx= blockIdx.x * blockDim.x + threadIdx.x; 
	
	float ac;//721? koliko puta ucitavas i gde  da mnozis...zasto float proveri koliko imas preracunavanja
    int angle; float fix1=3.14/(60*180); float fix2=57;
    
   
    __shared__ int mn[180], mn1[180], mn2[180], mn3[180], r[180], r1[180], r2[180], r3[180], s[180], s1[180], s2[180], s3[180];
   if(threadIdx.x==0 )
    {
        for (int i=0;i<720;i++)
	{ mn[i] = 0; r[i]=0;s[i]=0;  mn1[i] = 0; r1[i]=0;s1[i]=0; mn2[i] = 0; r2[i]=0;s2[i]=0; mn3[i] = 0; r3[i]=0;s3[i]=0;
	 } 
    }
    __syncthreads(); 
 

   if ( idx<100000)
    {
      
        for(int i=0; i<100000; i++)
        	{
		   
            ac= acosf((sin(b0[idx]*fix1)*sin(b1[i]*fix1))+ cos(b0[idx]*fix1)*cos(b1[i]*fix1)*cos((a1[i]-a0[idx])*fix1));
		ac= (ac*fix2/0.25); 
	
		angle=(int) ac; 
		if(ac>539) 
		  atomicAdd(&mn3[angle-540],1);
		else if(ac>359)
			atomicAdd(&mn2[angle-360],1);
		else if(ac>179)
			atomicAdd(&mn1[angle-179],1);
		else 
			atomicAdd(&mn[angle],1);
		}
		
	   for(int i=idx+1; i<100000;i++)
	    {  ac= acosf((sin(b0[idx]*fix1)*sin(b0[i]*fix1))+ cos(b0[idx]*fix1)*cos(b0[i]*fix1)*cos((a0[i]-a0[idx])*fix1));
	    ac= (ac*fix2/0.25); 
            angle=(int) ac; 
            if(ac>539) 
		  atomicAdd(&r3[angle-540],1);
		else if(ac>359)
			atomicAdd(&r2[angle-360],1);
		else if(ac>179)
			atomicAdd(&r1[angle-179],1);
		else 
			atomicAdd(&r[angle],1);
		
	     
            ac= acosf((sin(b1[idx]*fix1)*sin(b1[i]*fix1))+ cos(b1[idx]*fix1)*cos(b1[i]*fix1)*cos((a1[idx]-a1[i])*fix1));
            ac= (ac*fix2/0.25); 
	    angle=(int) ac; 
            if(ac>539) 
		  atomicAdd(&s3[angle-540],1);
		else if(ac>359)
			atomicAdd(&s2[angle-360],1);
		else if(ac>179)
			atomicAdd(&s1[angle-179],1);
		else 
			atomicAdd(&s[angle],1);
		

                }
	
           }
	

    __syncthreads();

      if(threadIdx.x==0)
    {
        for(int i=0;i<180;i++)
	{  
		hist[i+(blockIdx.x*720)]=mn[i]; hist_r[i+(blockIdx.x*720)]=r[i]; hist_s[i+(blockIdx.x*720)]=s[i];
		hist[i+180+(blockIdx.x*720)]=mn1[i]; hist_r[i+180+(blockIdx.x*720)]=r1[i]; hist_s[i+180+(blockIdx.x*720)]=s1[i];
		hist[i+360+(blockIdx.x*720)]=mn2[i]; hist_r[i+360+(blockIdx.x*720)]=r2[i]; hist_s[i+360+(blockIdx.x*720)]=s2[i];
		hist[i+540+(blockIdx.x*720)]=mn3[i]; hist_r[i+540+(blockIdx.x*720)]=r3[i]; hist_s[i+540+(blockIdx.x*720)]=s3[i];}
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
	int xx=530; 
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
       fscanf(synthetic_g, "%e %e", &h_A1[i], &h_B1[i]);}
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
	dim3 threadsPerBlock(192);
	
 
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
	
		/*for(int i=0; i<720; i++)
			brk+=(double) result[i]; 
	printf("%f\n ", brk);
	brk=0;
	for(int i=0; i<720; i++)
			brk+=result_s[i]; 
	printf("%f\n ", brk);
	brk=0;
	for(int i=0; i<720; i++)
			brk+=result_r[i]; 
	//printf("%f\n ", brk);
	brk=0;
		//{printf("%f ", final[i]);   }
	/*printf("\n druga\n " ) ; 
	for(int i=0; i<720; i++)
		{printf("%d ", result_r[i]);   }
	printf("\n treca\n " ) ; 
	for(int i=0; i<720; i++)
		{printf("%d ", result_s[i]);   }*/
	

cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
	cudaFree(h_A); cudaFree(h_B); cudaFree(h_C);cudaFree(d_D);cudaFree(h_D);cudaFree(d_E);cudaFree(h_E);
	cudaFree(d_A1);cudaFree(h_A1);cudaFree(d_B1);cudaFree(h_B1);

	
}








    
  
