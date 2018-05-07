#include <stdio.h>
#include <math.h>  
const int N = 16; 
const int blocksize = 16; 
 
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

void read_the_files()
{
	FILE *real_g, FILE *synthetic_g;
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
for(int i=0; i<galaxies_r; i++) printf("%f", a0[i]); 
}
