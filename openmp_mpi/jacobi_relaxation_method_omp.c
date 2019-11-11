#include <stdio.h>
#include <omp.h>
#include <math.h>

int main (int argc, char *argv[]) {

double m_pi=3.14159265358979323846;
int n=1500;	//size of matrix
double m[n+2][n+2],new[n][n],dx;
int ithread,iter, row,col,nproc,i,j;
double norm, mynorm;
int conv;
nproc=8;
conv=0;
iter=0;
norm=0;

dx=1./(double) (n+1) ;
for (i=1;i<=n;++i) {
	for (j=1;j<=n;++j)
		m[i][j]=0.;
	}

for (i=0;i<n+2;++i) {
	m[0][i]=sin(M_PI*(double)i*dx);
	m[n+1][i]=m[0][i]*exp(-M_PI);
	printf("aaaa %f \n",m[n+1][i]);
	}
for (i=0;i<n+2;++i) {
	m[i][0]=0.;
	m[i][n+1]=0.;
	}

for (i=1;i<=n;++i) {
	for (j=1;j<=n;++j)
	  m[i][j]=0.25*(m[i-1][j]+m[i+1][j]+m[i][j+1]+m[i][j-1]);
//	printf("fff %f",m[2][0]);
	}

//printf("55 %f \n",m[5][5]);

#pragma omp parallel shared(n,m,new,iter,nproc,conv,norm) private(ithread,row,col,mynorm)
{
	ithread=omp_get_thread_num();
	while (!conv) {
		mynorm=0.;
		for (row=ithread*(n/nproc); row<(ithread+1)*(n/nproc);++row) {
			for (col=0;col<n;++col) {
				new[row][col]=.25*(m[row][col+1]+m[row+1][col]+m[row+1][col+2]+m[row+2][col+1]);
//				printf("test %d %f    ",ithread, new[row][col]);
				mynorm+=pow(new[row][col]-m[row+1][col+1],2);
				}
			}
		#pragma omp atomic
		norm+=mynorm;
		#pragma omp barrier
		if (ithread==0) {
			++iter;
			if (norm < 1e-10)   conv=1;
			printf("norm %d %f \n",ithread,norm);
			norm=0;

		}
//		printf ("Thread %d done \n", ithread);
		for (row=ithread*(n/nproc); row<(ithread+1)*(n/nproc);++row){
			for (col=0;col<n;++col)
				m[row+1][col+1]=new[row][col];
		}
		#pragma omp barrier
	}
}

printf("Results with %d iters \n", iter);
for (row=0;row<=n+1;++row) {
	for (col=0;col<=n+1;++col)
		printf("%f  ",m[row][col]);
	printf("\n");
	}

}
