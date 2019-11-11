// This is a OpenMP parallel code that calculates solution to Laplace's 
// equation, using Jacobi's iterative relaxation method
// 
#include <stdio.h>
#include <omp.h>
#include <math.h>

int main (int argc, char *argv[]) {

double m_pi=3.14159265358979323846;
int n=1500;	                    //size of matrix
double dx;                          //grid physical dimension
double m[n+2][n+2];                 //Main matrix
double new[n][n];                   //Matrix holding updated values
int ithread,iter, row,col,nproc,i,j;
double norm, mynorm;
int conv;
nproc=8;
conv=0;
iter=0;
norm=0;

// Physical dimension of each grid
dx=1./(double) (n+1) ;

// Initialize matrix M
for (i=1;i<=n;++i) {
	for (j=1;j<=n;++j)
		m[i][j]=0.;
	}

// Apply boundary condition to edges of M
for (i=0;i<n+2;++i) {
	m[0][i]=sin(M_PI*(double)i*dx);
	m[n+1][i]=m[0][i]*exp(-M_PI);
	}
for (i=0;i<n+2;++i) {
	m[i][0]=0.;
	m[i][n+1]=0.;
	}
for (i=1;i<=n;++i) {
	for (j=1;j<=n;++j)
	  m[i][j]=0.25*(m[i-1][j]+m[i+1][j]+m[i][j+1]+m[i][j-1]);
	}


// Main parallel loop
// Repeat Jacobi iterations on all inner matrix elements until the solution
// converges. The solution's convergence is determined by checking if the sum
// over all the deltas of every matrix element's value between current and 
// previous iterations.
#pragma omp parallel shared(n,m,new,iter,nproc,conv,norm) private(ithread,row,col,mynorm)
{
    ithread=omp_get_thread_num();
    while (!conv) {
        mynorm=0.;
        // Each OMP thread is responsible for reading a few assigned rows of matrix m, and 
        // then writes updated values to shared matrix "new"
        for (row=ithread*(n/nproc); row<(ithread+1)*(n/nproc); ++row) {
            for (col=0; col<n; ++col) {
                new[row][col] = .25*(m[row][col+1] + m[row+1][col] + m[row+1][col+2] + m[row+2][col+1]);
                mynorm += pow(new[row][col] - m[row+1][col+1],2);
            }
        }

        // Require atomic update
        #pragma omp atomic
        norm+=mynorm;

        // Wait for all threads to finish update
        #pragma omp barrier

        if (ithread==0) {
            // Only main thread does final convergence check
            ++iter;
            if (norm < 1e-10)   conv=1;
            printf("norm %d %f \n",ithread,norm);
            norm=0;
        }

        // All threads help copy data from matrix "new" back to matrix m
        for (row=ithread*(n/nproc); row<(ithread+1)*(n/nproc); ++row){
            for (col=0;col<n;++col)
                m[row+1][col+1]=new[row][col];
        }

        // Wait for all threads to finish update
        #pragma omp barrier
	}
}

// Back to Main thread
// Print out entire matrix
printf("Results with %d iters \n", iter);
for (row=0;row<=n+1;++row) {
	for (col=0;col<=n+1;++col)
		printf("%f  ",m[row][col]);
	printf("\n");
	}

}
