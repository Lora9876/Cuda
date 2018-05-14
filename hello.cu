#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include <time.h>

#include<cuda_runtime.h>
#define PI 3.14159265
const int val= 180.0 / PI;
using namespace std; 
//const int N = 16; 
//const int blocksize = 16; 
//const int  SUBMATRIX_SIZE=16384 ;
const int thread= 256; 
const float fix =1/60 ; 
const int bins=720; 
__global__ void angles(volatile float *a0, volatile float *b0, volatile float *a1, volatile float *b1, int xind, int yind, int max_x, int max_y, volatile int *histi)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;// ovo proveri
	float angle; 
   	idx+=xind; 
	float arccos[180]={170.565836,167.031761,164.312253,162.013086,159.981700,158.139924,156.441357,154.855789,153.362253,151.945493,150.594000,149.298836,148.052901,146.850438,145.686705,144.557743,143.460201,142.391218,141.348325,140.329374,139.332483,138.355992,137.398427,136.458473,135.534951,134.626797,133.733048,132.852831,131.985347,131.129867,130.285721,129.452292,128.629010,127.815347,127.010813,126.214954,125.427343,124.647584,123.875305,123.110156,122.351810,121.599957,120.854305,120.114580,119.380521,118.651879,117.928421,117.209924,116.496173,115.786968,115.082114,114.381427,113.684729,112.991851,112.302629,111.616908,110.934536,110.255369,109.579267,108.906094,108.235721,107.568021,106.902872,106.240155,105.579755,104.921560,104.265462,103.611353,102.959132,102.308697,101.659949,101.012792,100.367132,99.722875,99.079931,98.438211,97.797626,97.158091,96.519520,95.881828,95.244934,94.608753,93.973206,93.338211,92.703689,92.069559,91.435743,90.802162,90.168737,89.535391,88.902044,88.268620,87.635039,87.001223,86.367093,85.732570,85.097575,84.462028,83.825848,83.188953,82.551261,81.912690,81.273155,80.632570,79.990850,79.347906,78.703650,78.057989,77.410832,76.762084,76.111649,75.459428,74.805320,74.149221,73.491026,72.830626,72.167909,71.502760,70.835060,70.164687,69.491514,68.815412,68.136245,67.453873,66.768152,66.078930,65.386052,64.689354,63.988667,63.283813,62.574608,61.860858,61.142360,60.418902,59.690260,58.956201,58.216476,57.470825,56.718972,55.960625,55.195477,54.423197,53.643438,52.855827,52.059968,51.255435,50.441772,49.618489,48.785060,47.940914,47.085434,46.217951,45.337733,44.443985,43.535830,42.612308,41.672354,40.714790,39.738299,38.741408,37.722457,36.679563,35.610580,34.513038,33.384076,32.220343,31.017880,29.771945,28.476782,27.125288,25.708528,24.214992,22.629425,20.930857,19.089081,17.057696,14.758528,12.039020,8.504946,0.000005};

	
 	__shared__ unsigned int shared[bins];
    		// za prvu petlju ocistis uvek
    			if(threadIdx.x==0)
   			 {
       			 for (int i=0;i<bins;i++)
         		   shared[i] = 0;
   			 }
	
   	 __syncthreads();

	
      //provera
			for(int i=yind; i<max_y; i++)
       			
			{ angle = sin(b0[idx]*fix) *sin(b1[i]*fix) + cos(b0[idx]*fix) * cos(b1[i]*fix) * cos(fix*a0[idx]*-fix*a1[i]);
				shared[int(arccos[int(angle)])]++; 	
				
			//nadji nacin da atomic add proradi :D
			//atomicAdd(&shared[int(angle)],1); 
				}
   	 __syncthreads();

    if(threadIdx.x==0)
    {
        for(int i=0;i<bins;i++)
            histi[i]=shared[i];
    }



}


 void read_the_files()
{
	//reading files 1 and 2 
	FILE *real_g; FILE *synthetic_g;
	int galaxies_r, galaxies_s; 
	float *a0, *a1, *b0, *b1;
	real_g = fopen("data_100k_arcmin.txt","r");
    	synthetic_g = fopen("flat_100k_arcmin.txt","r");	
	 fscanf(real_g, "%d", &galaxies_r);
	 fscanf(synthetic_g,  "%d", &galaxies_s);
	
	
	a0= (float*) malloc(galaxies_r* sizeof(float));
	b0= (float*) malloc(galaxies_r* sizeof(float)); 
	a1= (float*) malloc(galaxies_s* sizeof(float)); 
	b1= (float*) malloc(galaxies_s* sizeof(float)); 
	for(int i=0; i<galaxies_r; i++)
    {
       
        fscanf(real_g, "%e %e", &a0[i], &b0[i]);
       fscanf(synthetic_g, "%e %e", &a1[i], &b1[i]);
    }	
	/*for(int i=0; i<galaxies_r; i++)
	{
		a0[i]=a1[i]=b0[i]=b1[i]=0.5; 
	
	}*/
	 fclose(real_g);
	 fclose(synthetic_g);
//for(int i=0; i<galaxies_r; i++) printf("%f", a0[i]);
	 
    dim3 grid, block;
    
    grid.x = 1024; 
	 grid.y=1024;
	// grid.y=1024; 
    block.x = 1; 
	 float *aa1, *bb1, *aa0, *bb0; 
	 
    cudaMalloc((void **) &aa0, galaxies_r* sizeof(float));
    cudaMalloc((void **) &bb0, galaxies_r* sizeof(float));

    cudaMalloc((void **) &aa1, galaxies_s* sizeof(float));
    cudaMalloc((void **) &bb1, galaxies_s* sizeof(float) );

    // dovoljno memorije?
    

    // Initialize array to all 0's
    cudaMemset(aa0,0,galaxies_r* sizeof(float));
    cudaMemset(bb0,0,galaxies_r* sizeof(float));
    cudaMemset(aa1,0,galaxies_s* sizeof(float));
    cudaMemset(bb1,0,galaxies_s* sizeof(float));

    cudaMemcpy(aa0, a0, galaxies_r* sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy(bb0, b0,galaxies_r* sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy(aa1, a1, galaxies_s* sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy(bb1, b1,galaxies_s* sizeof(float), cudaMemcpyHostToDevice );

    int x, y;

	 /*  int num_x = galaxies_r/ SUBMATRIX_SIZE;
    int num_y = galaxies_s / SUBMATRIX_SIZE;
    
    // Take care of edges of matrix.
    if (galaxies_r%SUBMATRIX_SIZE != 0)
    {
        num_x ++;
    }
    if (galaxies_s%SUBMATRIX_SIZE != 0)
    {
        num_y ++;
     */
	 
	 //preparing the histogram array 
	 int *hist, *histi , *tmp; 
	 
   
    int size_h_bytes = 720*sizeof(int);

    hist = (int*)malloc(size_h_bytes);
    memset(hist, 0, size_h_bytes);

   
    cudaMalloc((void **) &tmp, (size_h_bytes));
    cudaMemset(tmp, 0, size_h_bytes);

    unsigned long  *hist_array;

    int hist_array_size = 720 * sizeof(unsigned long);
    hist_array =  (unsigned long*)malloc(hist_array_size);
  
    memset(hist_array,0,hist_array_size); 
	 cudaMemset(tmp, 0,size_h_bytes);
	 
	   angles<<<grid,block>>>(aa0, bb0, aa1, bb1, 0, 0, 512, 512, tmp);
            cudaMemcpy(hist, tmp, size_h_bytes, cudaMemcpyDeviceToHost);
	 
	 for(int i=0; i<720; i++)
		printf("%d ", hist[i]);
	 
	 
    free(a1);
    free(b1);
    free(a0);
    free(b0);

    cudaFree(aa1);
    cudaFree(aa0);  
    cudaFree(bb0);
    cudaFree(bb1);  
    cudaFree(tmp);

 }
 //prepration for the kernel
	 


int main()
{
	float alpha1= 4646.98;
	float b1= 3749.51;
	float a2=4644.35; 
	float b2=3749.52;
	
	float theta1= acos(sin(b1)*sin(b2) + cos(b1)*cos(b2) *cos(alpha1-a2));
	
	printf("%f\n", b1);
	printf("%f\n", theta1);
	
	
	 clock_t start, end;
     double cpu_time_used;
     
     start = clock();
    
 
		
 	//read_the_files(); 
	// do some calculations
 	end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("%f", cpu_time_used); 
	return EXIT_SUCCESS;
}


