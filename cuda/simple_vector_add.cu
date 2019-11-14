#include <stdio.h>

const int BLOCKSIZE = 128; 

__global__
void vecAddKernel(float* A, float* B, float* C, int n) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i < n) C[i] = A[i] + B[i];
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
    vecAddKernel<<<ceil(1.0*n/BLOCKSIZE), BLOCKSIZE>>>(d_A, d_B, d_C, n);

    cudaMemcpy(C, d_C, size, cudaMemcpyDeviceToHost);

    cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
}


int main () {
    float A[5] = {1, 2, 3, 4, 5};
    float B[5] = {5, 4, 3, 2, 1};
    float C[5];

    vecAdd(A, B, C, 5);

    printf("Elements: %f %f %f %f %f\n", C[0], C[1], C[2], C[3], C[4]);
}
