#include <stdio.h>
#include <omp.h>
#include <math.h>

void calculate(int ithread, int npair, double* temp)
{
  int j;
  //  *temp=*temp+ithread;
	for (j=ithread*npair;j<(ithread+1)*npair; j=j+1)
	{
	  *temp=*temp+1.0/(2.*2.*j+1.)-1./(2.*2.*j+3.);
	printf("check%d %f \n",ithread,*temp);
	}
}

int main (int argc, char *argv[]) {

int n, nthread, ithread, iter, chunk,i,j,npair,nterms;
double result, timei, timef,temp;

npair=100000*2;
result=0.0;
//for (i=1;i<=10;i+=1)
//{
#pragma omp parallel private(temp)  shared(npair) reduction(+:result)
	{
	temp=0.;
//	nthread=omp_get_num_threads();
	nthread=8;
	ithread=omp_get_thread_num();
	//	printf("thread %d \n",ithread);
	//	printf("npair %d \n",npair);
	nterms=npair/nthread;
	//	printf ("nterms %d \n",nterms);
	calculate(ithread,nterms,&temp);
	result=result+temp;
	}
	printf("Pi = %f \n", result*4.);

// }
}

