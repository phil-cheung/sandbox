// This is a serial code (non-parallelized) calculates solution to Laplace's 
// equation, using Jacobi's iterative relaxation method
// 
#include <stdio.h>
#include <omp.h>
#include <math.h>

int main (int argc, char *argv[]) {

double m_pi=3.14159265358979323846;
int n=160;	             //size of matrix
double dx;                   //grid physical dimension
double m[n+2][n+2];          //Main matrix
double new[n][n];            //Matrix holding updated values
int ithread,iter, row,col,nproc,i,j;
double norm, mynorm;
int conv;
conv=0;
iter=0;
norm=0;

// Physical dimension of each grid
dx=1./(double) (n+1) ;

// Initialize matrix M
for (i=1; i<=n; ++i) {
    for (j=1; j<=n; ++j)
        m[i][j]=0.;
}

// Apply boundary condition to edges of M
for (i=0; i<n+2; ++i) {
    m[0][i]=sin(M_PI*(double)i*dx);
    m[n+1][i]=m[0][i]*exp(-M_PI);
}
for (i=0;i<n+2;++i) {
    m[i][0]=0.;
    m[i][n+1]=0.;
}
for (i=1; i<=n; ++i) {
    for (j=1; j<=n; ++j)
        m[i][j] = 0.25*( m[i-1][j] + m[i+1][j] + m[i][j+1] + m[i][j-1]);
}

// Main loop
// Repeat Jacobi iterations on all inner matrix elements until the solution
// converges. The solution's convergence is determined by checking if the sum
// over all the deltas of every matrix element's value between current and 
// previous iterations.
{
    while (!conv) {
        norm=0.;
        for (row=0; row<n;++row) {
            for (col=0;col<n;++col) {
                // m is the main matrix, "new" is a temp matrix used to store results of next
                // iteration
                new[row][col]=.25*(m[row][col+1]+m[row+1][col]+m[row+1][col+2]+m[row+2][col+1]);
                norm+=pow(new[row][col]-m[row+1][col+1],2);
            }
        }
        ++iter;
        if (norm < 1e-10)   conv=1;
        printf("norm %d %f \n",ithread,norm);
        norm=0;

        // Copy elements back from "new" matrix to "m" matrix
        for (row=0; row<n;++row){
            for (col=0;col<n;++col)
                m[row+1][col+1]=new[row][col];
            }
        }
}

// Print out entire matrix
printf("Results with %d iters \n", iter);
for (row=0;row<=n+1;++row) {
	for (col=0;col<=n+1;++col)
		printf("%f  ",m[row][col]);
	printf("\n");
	}

}
