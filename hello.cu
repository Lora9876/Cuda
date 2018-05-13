#include <stdio.h>
#include <math.h>  
//const int N = 16; 
//const int blocksize = 16; 
//const int  SUBMATRIX_SIZE=16384 ;
const int thread= 256; 
const int bins=720; 
__global__ void angles(volatile float *a0, volatile float *b0, volatile float *a1, volatile float *b1, int xind, int yind, int max_x, int max_y, volatile int *histi)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x; // ovo proveri
	float angle; 
    idx += xind;
 __shared__ int shared[bins];
    // za prvu petlju ocistis uvek
    if(threadIdx.x==0)
    {
        for (int i=0;i<bins;i++)
            shared[i] = 0;
    }
    __syncthreads();

	 if (idx<max_x)
    {
       
      //provera

        for(int i=yind; i<max_y; i++)
        {
           
                
		
                angle = sin(b0[idx]) *sin(b1[i]) + cos(b0[idx]) * cos(b1[i]) * cos(a0[idx]-a1[0]);
		//how to put angle
               
		atomicAdd(&shared[angle],1); 
	
	}
	
	
	

          

    __syncthreads();

    if(threadIdx.x==0)
    {
        for(int i=0;i<bins;i++)
            histi[i]=shared[i];
    }



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
       
        fscanf(real_g, "%f %f", &a0[i], &b0[i]);
       fscanf(synthetic_g, "%f %f", &a1[i], &b1[i]);
    }		    
//for(int i=0; i<galaxies_r; i++) printf("%f", a0[i]);
	 
    dim3 grid, block;
    
    grid.x = thread; 
    block.x = 1024; 
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

	 {  int num_x = galaxies_r/ SUBMATRIX_SIZE;
    int num_y = galaxies_s / SUBMATRIX_SIZE;
    
    // Take care of edges of matrix.
    if (galaxies_r%SUBMATRIX_SIZE != 0)
    {
        num_x ++;
    }
    if (galaxies_s%SUBMATRIX_SIZE != 0)
    {
        num_y ++;
    }
	 }
	 //preparing the histogram array 
	 int *hist, *histi , *tmp;

   
    int size_h_bytes = 720*sizeof(int);

    hist = (int*)malloc(size_h_bytes);
    memset(hist, 0, size_h_bytes);

   
    cudaMalloc((void **) &histi, (size_h_bytes));
    cudaMemset(histi, 0, size_h_bytes);

    unsigned long  *hist_array;

    int hist_array_size = 720 * sizeof(unsigned long);
    hist_array =  (unsigned long*)malloc(hist_array_size);
  
    memset(hist_array,0,hist_array_size); 
	 cudaMemset(tmp,0,size_h_bytes);
	   angles<<<grid,block>>>(aa0, bb0,aa1, bb1, 0, 0, 512, 512, tmp);
            cudaMemcpy(hist, tmp, size_h_bytes, cudaMemcpyDeviceToHost);
	 
	 for(int i=0; i<720; i++)
		printf("%d", hist[i]);
 }
 //prepration for the kernel
	 
}
//__global__ 
/*void hello(char *a, char *b) 
{
	a[threadIdx.x] += b[threadIdx.x];
}
 */
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
int main()
{
	float alpha1= 4646.98;
	float b1= 3749.51;
	float a2=4644.35; 
	float b2=3749.52;
	
	float theta1= acos(sin(b1)*sin(b2) + cos(b1)*cos(b2) *cos(alpha1-a2));
	
	printf("%f\n", b1);
	printf("%f\n", theta1);
	getDeviceDiagnostics(20000, 2);
 	read_the_files(); 
	/*char *ad;
	int *bd;
	const char csize = N*sizeof(char);
	const char isize = N*sizeof(char);
 
	printf("%s", a);
 
	cudaMalloc( (void**)&ad, csize ); 
	cudaMalloc( (void**)&bd, isize ); 
	cudaMemcpy( ad, a, csize, cudaMemcpyHostToDevice ); 
	cudaMemcpy( bd, b, isize, cudaMemcpyHostToDevice ); 
	
	dim3 dimBlock( blocksize, 1 );
	dim3 dimGrid( 1, 1 );
	hello<<<dimGrid, dimBlock>>>(ad, bd);
	cudaMemcpy( a, ad, csize, cudaMemcpyDeviceToHost ); 
	cudaFree( ad );
	cudaFree( bd );
	
	printf("%s\n", a);*/
	return EXIT_SUCCESS;
}


