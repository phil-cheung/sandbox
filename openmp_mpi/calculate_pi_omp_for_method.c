// Calculate Pi using Taylor series expansion
// This code uses OpenMP's parallel-for loop to automatically distribute partial sum
// caclulations to the threads
#include <stdio.h>
#include <omp.h>
#include <math.h>

int main (int argc, char *argv[]) {

int n, nthread, ithread, iter, chunk,i,j,npair,nterms,ncore;
double result, timei, timef;

n=10;
chunk=300;

// We add two partial sum terms at a time, so NPAIR is the total num of pairs of partial sum
// to do in the outer-for-loop. For each value of NPAIR, we have an inner-for-loop to repeat
// doing the summing by int(10,000,000 / n) times so we can measure the avg run time
for (npair=2; npair<=20100; npair*=2) {
    result=0.0;
    timei=omp_get_wtime();
    n=10000000/npair;

    // Here we use standard C for-loop to repeat the OMP code block, in order to re-establish
    // the OMP parallel environment so we can repeat the calculation multiple times for run
    // time measurement purpose
    for (i=1; i<=n; i+=1) {
        #pragma omp parallel private(ithread,ncore)
        {
            ithread=omp_get_thread_num();
            nthread=omp_get_num_threads();

            // Distribute work using OMP parallel for-loop
            #pragma omp parallel for schedule(dynamic,chunk) shared(npair) reduction(+:result)
            for (j=0;j<npair; j=j+1)
	         result=result+ 1./(4.*j+1)-1./(4.*j+3);
        }
    }

    // Back to Main thread work
    timef=omp_get_wtime();
    printf("Pi = %f \n", result*4./n);
    printf("Npair= %d, Time= %f s \n", npair,(timef-timei)/n);
}

}
