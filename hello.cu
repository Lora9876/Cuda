#include <stdio.h>
#include <math.h>  
const int N = 16; 
const int blocksize = 16; 
const int  SUBMATRIX_SIZE=16384 ;
const int thread= 256; 
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
    
    grid.x = 8192/thread; 
    block.x = SUBMATRIX_SIZE/grid.x; 
	 
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

    int num_x = galaxies_r/ SUBMATRIX_SIZE;
    int num_y = galaxies_s / SUBMATRIX_SIZE;
    
    // Take care of edges of matrix.
    if (galaxies_r%SUBMATRIX_SIZE != 0)
    {
        num_submatrices_x ++;
    }
    if (galaxies_s%SUBMATRIX_SIZE != 0)
    {
        num_submatrices_y ++;
    }

	 //preparing the histogram array 
	 int *hist, *dev_h;

    int size_h = SUBMATRIX_SIZE * thread;
    int size_h_bytes = size_h*sizeof(int);

    hist = (int*)malloc(size_h_bytes);
    memset(hist, 0, size_h_bytes);

   
    cudaMalloc((void **) &dev_h, (size_h_bytes));
    cudaMemset(dev_h, 0, size_h_bytes);

    unsigned long  *hist_array;

    int hist_array_size = threads * sizeof(unsigned long);
    hist_array =  (unsigned long*)malloc(hist_array_size);
    printf("Size of histogram array: %d bytes\n",hist_array_size);
    memset(hist_array,0,hist_array_size); 
	 
 //prepration for the kernel
	 
}
//__global__ 
/*void hello(char *a, char *b) 
{
	a[threadIdx.x] += b[threadIdx.x];
}
 */
int main()
{
	float alpha1= 4646.98;
	float b1= 3749.51;
	float a2=4644.35; 
	float b2=3749.52;
	
	float theta1= acos(sin(b1)*sin(b2) + cos(b1)*cos(b2) *cos(alpha1-a2));
	
	printf("%f\n", b1);
	printf("%f\n", theta1);
	
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


