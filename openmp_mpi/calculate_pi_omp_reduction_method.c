// Calculate Pi using Taylor series expansion
// This code uses OpenMP's reduction method to sum over the partial sums
// done by each OMP thread
#include <stdio.h>
#include <omp.h>
#include <math.h>

void calculate(int ithread, int npair, double* temp)
{
    int j;
    for (j=ithread*npair; j < (ithread+1)*npair; j=j+1) {
        *temp = *temp + 1.0/(2.*2.*j+1.) - 1./(2.*2.*j+3.);
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
    // nthread=omp_get_num_threads();
    nthread=8;
    ithread=omp_get_thread_num();
    // printf("thread %d \n",ithread);
    // printf("npair %d \n",npair);
    nterms=npair/nthread;
    // printf ("nterms %d \n",nterms);

    // Each thread calls function calculate() which does partial summing and
    // stores into each thread's own variable "temp"
    calculate(ithread, nterms, &temp);

    // Gathering contributions from each thread into variable "result"
    result=result+temp;
}

// Back to Main thread
printf("Pi = %f \n", result*4.);
}

