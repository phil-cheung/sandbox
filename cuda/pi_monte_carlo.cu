#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <math.h>

const int GRIDSIZE  = 4;
const int BLOCKSIZE = 128; 

const int N_PER_THREAD = 1000000;

__global__
void PiKernel(int* N) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    curandState_t state;
    int my_n = 0;

    curand_init(1234, i, 0, &state);
    float r1, r2;
    for (int j = 0; j < N_PER_THREAD; j++) {
        //float r1 = rand() / RAND_MAX;
        //float r2 = rand() / RAND_MAX;
        //r1 = (curand(&state) % 10000) / 10000.0;
        //r2 = (curand(&state) % 10000) / 10000.0;

        r1 = curand_uniform(&state);
        r2 = curand_uniform(&state);
 
        if (sqrt(r1*r1 + r2*r2) <= 1.0) {
            my_n++;
        }
    }

    //Watch out for overflow
    atomicAdd(N, my_n);    
}


int main () {
    int N = 0;
    int *d_N;

    cudaMalloc((void**) &d_N, sizeof(int));
    cudaMemcpy(d_N, &N, sizeof(int), cudaMemcpyHostToDevice);
    PiKernel<<<GRIDSIZE, BLOCKSIZE>>>(d_N);
    cudaMemcpy(&N, d_N, sizeof(int), cudaMemcpyDeviceToHost);
    float pi = 4.0 * N / GRIDSIZE / BLOCKSIZE / N_PER_THREAD;

    printf("N is %i\n", N);
    printf("Pi is %f\n", pi);
}
