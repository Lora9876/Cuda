#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include <cuda_runtime.h>

#define PI 3.14159265

const int val = 180.0 / PI;
const int thread = 256;
const float fix = 1 / 60 * PI / 180;
const int bins = 720;

using namespace std;

__global__ void angles(volatile float *a0, volatile float *b0, volatile float *a1, volatile float *b1, int xind, int yind, int max_x, int max_y, volatile int *histi)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	float angle;
	idx += xind;
	
	__shared__ unsigned int shared[bins];

	// Always clear the first iteration
	if (threadIdx.x == 0) {
		for (int i = 0; i < bins; i++) {
			shared[i] = 0;
		}
	}
	
	__syncthreads();

	// Check
	for (int i = yind; i < yind + max_y; i++) {
		angle = acosf((sin(b0[idx]) * sin(b1[i]) + cos(b0[idx]) * cos(b1[i]) * cos(a0[idx] * -a1[i])) * 1 / val);
		shared[int(angle)]++;

		// TODO: Make atomicAddwork :D
		// atomicAdd(&shared[int(angle)], 1);
	}

	__syncthreads();

	if (threadIdx.x == 0) {
		for (int i = 0; i < bins; i++) {
			histi[i] = shared[i];
		}
	}
}

void read_the_files()
{
	// Reading files 1 and 2
	FILE *real_g;
	FILE *synthetic_g;
	int galaxies_r, galaxies_s;
	float *a0, *a1, *b0, *b1;

	real_g = fopen("data_100k_arcmin.txt", "r");
	synthetic_g = fopen("flat_100k_arcmin.txt", "r");
	fscanf(real_g, "%d", &galaxies_r);
	fscanf(synthetic_g, "%d", &galaxies_s);
	
	a0 = (float*) malloc(galaxies_r * sizeof(float));
	b0 = (float*) malloc(galaxies_r * sizeof(float));
	a1 = (float*) malloc(galaxies_s * sizeof(float));
	b1 = (float*) malloc(galaxies_s * sizeof(float));

	for (int i = 0; i < galaxies_r; i++) {
		fscanf(real_g, "%e %e", &a0[i], &b0[i]);
		fscanf(synthetic_g, "%e %e", &a1[i], &b1[i]);

		a0[i] *= fix;
		b0[i] *= fix;
		a1[i] *= fix;
		b1[i] *= fix;
	}

	fclose(real_g);
	fclose(synthetic_g);

	dim3 grid, block;

	grid.x = 1024;
	grid.y = 1024;
	block.x = 1;
	float *aa1, *bb1, *aa0, *bb0;

	cudaMalloc((void **) &aa0, galaxies_r * sizeof(float));
	cudaMalloc((void **) &bb0, galaxies_r * sizeof(float));

	cudaMalloc((void **) &aa1, galaxies_s * sizeof(float));
	cudaMalloc((void **) &bb1, galaxies_s * sizeof(float));

	// Do we have enough memory?

	// Initialize array to all 0's
	cudaMemset(aa0, 0, galaxies_r * sizeof(float));
	cudaMemset(bb0, 0, galaxies_r * sizeof(float));
	cudaMemset(aa1, 0, galaxies_s * sizeof(float));
	cudaMemset(bb1, 0, galaxies_s * sizeof(float));

	cudaMemcpy(aa0, a0, galaxies_r * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(bb0, b0, galaxies_r * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(aa1, a1, galaxies_s * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(bb1, b1, galaxies_s * sizeof(float), cudaMemcpyHostToDevice);

	int x, y;

	// Preparing the histogram array
	int *hist, *histi, *tmp;

	int size_h_bytes = 720 * sizeof(int);

	hist = (int*)malloc(size_h_bytes);
	memset(hist, 0, size_h_bytes);

	cudaMalloc((void **) &tmp, (size_h_bytes));
	cudaMemset(tmp, 0, size_h_bytes);

	unsigned long *hist_array;

	int hist_array_size = 720 * sizeof(unsigned long);
	hist_array = (unsigned long*)malloc(hist_array_size);

	memset(hist_array, 0, hist_array_size);
	cudaMemset(tmp, 0, size_h_bytes);

	for (int i = 0; i < 9; i++) {
		y = i * 512 * 512;

		for (int j = 0; j < 9; j++) {
			x = 512 * 512 * i;

			angles<<<grid,block>>>(aa0, bb0, aa1, bb1, x, y, 512, 512, tmp);
			cudaMemcpy(hist, tmp, size_h_bytes, cudaMemcpyDeviceToHost);

			for(int i = 0; i < 720; i++) {
				hist_array[i] += hist[i];
			}
		}
	}

	for(int i = 0; i < 720; i++) {
		printf("%d ", hist_array[i]);
	}

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

// Preparation for the kernel

int main()
{
	float alpha1 = 4646.98;
	float b1 = 3749.51;
	float a2 = 4644.35;
	float b2 = 3749.52;

	float theta1 = acos(sin(b1) * sin(b2) + cos(b1) * cos(b2) * cos(alpha1 - a2));

	printf("%f\n", b1);
	printf("%f\n", theta1);

	clock_t start, end;
	double cpu_time_used;

	start = clock();
	read_the_files();
	// Do some calculations
	end = clock();

	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("%f", cpu_time_used);

	return EXIT_SUCCESS;
}
