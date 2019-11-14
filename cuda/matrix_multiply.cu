/* This code parallelize matrix multiplication using the tiled method */
#include <stdio.h>

const int WIDTH = 16;           /* Width of square matrix */
const int TILE_WIDTH = 8;      /* width of 2D Blocks */

__global__ void MatrixMulKernel(float* d_M, float* d_N, float* d_P, int Width) {
    __shared__ float Mds[TILE_WIDTH][TILE_WIDTH];
    __shared__ float Nds[TILE_WIDTH][TILE_WIDTH];

    int bx = blockIdx.x;   int by = blockIdx.y;
    int tx = threadIdx.x;  int ty = threadIdx.y;

    /* Identify the rorw and column of the d_P element to work on */
    int Row = by * TILE_WIDTH + ty;
    int Col = bx * TILE_WIDTH + tx;

    float Pvalue = 0;

    /* Loop over the d_M and d_N tiles required to compute the d_P elements */
    for (int m = 0; m < Width/TILE_WIDTH; ++m) {
        /* Collaborative loading of d_M and d_N tiles into shared memory */
        /* After each iteration, the tile for M shifts right along the row */
        /* direction, and the tile for N shift down along the column direction */
        Mds[tx][ty] = d_M[Row*Width + m*TILE_WIDTH + tx];
        Nds[tx][ty] = d_N[(m*TILE_WIDTH + ty)*Width + Col];
        __syncthreads();

        for (int k = 0; k < TILE_WIDTH; ++k)
            Pvalue += Mds[tx][k] * Nds[k][ty];
        __syncthreads();
     }
     d_P[Row*Width + Col] = Pvalue;
}

void matrixMultiply(float* M, float* N, float* P, int width) {
    int size = width * width * sizeof(float);
    float *d_M, *d_N, *d_P;

    cudaMalloc((void **) &d_M, size);
    cudaMalloc((void **) &d_N, size);
    cudaMalloc((void **) &d_P, size);

    cudaMemcpy(d_M, M, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_N, N, size, cudaMemcpyHostToDevice);

    /* Calcuate block sizes */
    int numBlocks = ceil(1.0*width/TILE_WIDTH);
    dim3 dimGrid(numBlocks, numBlocks, 1);
    dim3 dimBlock(TILE_WIDTH, TILE_WIDTH, 1);

    MatrixMulKernel<<<dimGrid, dimBlock>>>(d_M, d_N, d_P, width);

    cudaMemcpy(P, d_P, size, cudaMemcpyDeviceToHost);

    cudaFree(d_M); cudaFree(d_N); cudaFree(d_P);
}

void printMatrix(float* M) {
    for (int i = 0; i < WIDTH; i++) {
        for (int j = 0; j < WIDTH; j++) {
            printf("%3.1f ", M[i*WIDTH + j]);
        }
        printf("\n");
    }
}

int main () {
    //float A[5] = {1, 2, 3, 4, 5};
    float M[WIDTH*WIDTH], N[WIDTH*WIDTH], P[WIDTH*WIDTH];

    /* Initialize values */
    int index;
    for (int i = 0; i < WIDTH; i++)
       for (int j = 0; j < WIDTH; j++) {
           index = i*WIDTH + j;
           if (i == j) {
               M[index] = 1.0;
               N[index] = 1.0;
           } else {
               M[index] = 0.0;
               N[index] = 0.0;
           }
       }

    matrixMultiply(M, N, P, WIDTH);

    printf("Result of M x N:\n");
    printMatrix(M);
    printMatrix(N);
    printMatrix(P);
    //printf("Elements: %f %f %f %f %f\n", C[0], C[1], C[2], C[3], C[4]);
}
