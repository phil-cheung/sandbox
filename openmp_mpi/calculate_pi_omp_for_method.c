#include <stdio.h>
#include <omp.h>
#include <math.h>

int main (int argc, char *argv[]) {

int n, nthread, ithread, iter, chunk,i,j,npair,nterms,ncore;
double result, timei, timef;

n=10;
chunk=300;
for (npair=2; npair<=20100; npair*=2)
{
result=0.0;
timei=omp_get_wtime();
n=10000000/npair;
for (i=1;i<=n;i+=1)
{
//{
#pragma omp parallel private(ithread,ncore)
 {
   ithread=omp_get_thread_num();
   nthread=omp_get_num_threads();
//   nthread=-1;
// printf ("thread %d num threads: %d \n ",ithread, nthread);
// ncore=14; nthread=1222;
#pragma omp parallel for schedule(dynamic,chunk)shared(npair) reduction(+:result)
	for (j=0;j<npair; j=j+1)
	  result=result+ 1./(4.*j+1)-1./(4.*j+3);
 }
}
timef=omp_get_wtime();
	printf("Pi = %f \n", result*4./n);
	printf("Npair= %d, Time= %f s \n", npair,(timef-timei)/n);
}
}
