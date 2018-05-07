#include <stdio.h>
#include <math.h>  
const int N = 16; 
const int blocksize = 16; 
 
__global__ 
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
	
	float theta1= arccos(sin(b1)*sin(b2) + cos(b1)*cos(b2) *cos(alpha1-a2));
	printf("%f\n", theta1);
 
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
