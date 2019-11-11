#include "mpi.h"
#include <stdio.h>
#include <math.h>

int main(int argc,char *argv[])
{
    int size,rank;
    int npair,j,i,iter,n,nprocs;
    double partialsum, sum,timei,timef;

//    for (nprocs=2;nprocs<65;nprocs*=2)
//    {
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (n=size;n<1000000;n*=3)
    {
    npair=n/size;
    sum=0.;
    partialsum=0.;
    timei=MPI_Wtime();
    iter=10000000/n;
    for (i=0;i<iter;i+=1)
    {
    for (j=rank*npair;j<(rank+1)*npair; j=j+1)
	partialsum=partialsum+ 1./(2.*2.*(double)j+1.)-1./(4.*(double)j+3.);
    MPI_Reduce(&partialsum,&sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    }
    timef=MPI_Wtime();
    if (rank==0) {
	printf("nprocs= %d Npair= %d Time = %f , Pi = %f \n",size, npair, (timef-timei)/iter,4.*sum/(double)iter);
    }
    }
    MPI_Finalize();

    return 0;
}
