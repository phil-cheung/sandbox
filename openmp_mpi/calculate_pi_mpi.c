// Calculate Pi using Taylor series expansion
// This code uses MPI to distribute partial sum caclulations to the nodes, and 
// then do MPI_Reduce to get the total
#include "mpi.h"
#include <stdio.h>
#include <math.h>

int main(int argc,char *argv[])
{
    int size,rank;
    int npair,j,i,iter,n,nprocs;
    double partialsum, sum,timei,timef;

    // Standard MPI initialization
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Main loop
    // In the outer for-loop, N is the total num of series terms we want to distribute to the 
    // nodes. For each value of N, we have the inner for-loop to have each node repeat the same
    // summing of partial sums by int(10,000,000 / n) times, this is to get a more accurate 
    // measurement of run time
    for (n=size; n<1000000; n*=3) {
        npair=n/size;
        sum=0.;
        partialsum=0.;
        timei=MPI_Wtime();
        iter=10000000/n;

        // Each node does a chunk of partial sums below, and is repeated by "iter" num of
        // times to get a time avg
        for (i=0; i<iter; i+=1) {
            for (j=rank*npair; j<(rank+1)*npair; j=j+1)
                partialsum=partialsum+ 1./(2.*2.*(double)j+1.)-1./(4.*(double)j+3.);

            // Intentionally do MPI_Reduce muliple times, so we can include its overhead in
            // run time measurement
            MPI_Reduce(&partialsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        timef=MPI_Wtime();

        // Have only main node report statistics
        if (rank==0) {
            printf("nprocs= %d Npair= %d Time = %f , Pi = %f \n",size, npair, (timef-timei)/iter,4.*sum/(double)iter);
        }
    }

    MPI_Finalize();

    return 0;
}
