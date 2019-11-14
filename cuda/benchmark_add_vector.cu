#include <stdio.h>
#include <unistd.h>
#include <time.h>

const int BLOCKSIZE = 128; 

__global__
void vecAddKernel(float* A, float* B, float* C, int n) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j;
    float k;

    for (j = 0; j < 10000000000; j++) {
        k = 3.14 * 2.2223;
    }
}


void vecAdd(float* A, float* B, float* C, int n) {
    int size = n * sizeof(float);
    float *d_A, *d_B, *d_C;

    cudaMalloc((void **) &d_A, size);
    cudaMemcpy(d_A, A, size, cudaMemcpyHostToDevice);
    cudaMalloc((void **) &d_B, size);
    cudaMemcpy(d_B, B, size, cudaMemcpyHostToDevice);
    cudaMalloc((void **) &d_C, size);

    //Call Kernel function to add vector
    time_t a = time(NULL);
    vecAddKernel<<<ceil(1.0*n/BLOCKSIZE), BLOCKSIZE>>>(d_A, d_B, d_C, n);
    time_t b = time(NULL);
    printf("Time diff is %ld\n", b - a);


    cudaMemcpy(C, d_C, size, cudaMemcpyDeviceToHost);

    cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
}


int main () {
    float A[5] = {1, 2, 3, 4, 5};
    float B[5] = {5, 4, 3, 2, 1};
    float C[5];

    vecAdd(A, B, C, 5);

}
