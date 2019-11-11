// This is a MPI + OpenMP hybrid parallel code that calculates solution to Laplace's 
// equation, using Jacobi's iterative relaxation method
#include "mpi.h"
#include <omp.h>
#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[])
{
    // Variables
    double M_PI=3.14159265358979323846;
    int size,rank,ithread,nthread;
    int i,j,conv,a,b,iter,row,col;
    double timei,timef,norm,mynorm,tempnorm,dx;
    MPI_Request req[4],req2[4];
    MPI_Status status[4],status2[4];

    // Standard MPI initializations
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // Other variables
    int nrow=8;                    //Each node holds 8 rows
    int n=nrow*size;               //Dimension of full matrix is n x n;

    // The full matrix is broken down into chunks of rows, and each node has their own
    // matrices below to hold the chunk of rows they are responsible for
    double m[nrow+2][n];           //Main matrix,   row 0 & nrow+1 are buffers for adjacent nodes data
    double new[nrow+2][n];         //Update matrix, row 0 & nrow+1 are buffers for adjacent nodes data

    printf("size is %d \n",size);

    // Initialize some variables
    conv=0;
    iter=0;
    dx=1./(double) (n-1) ;     //Physical grid dimension
    for (row=1; row<nrow+1; ++row){
        m[row][0]=0.;   new[row][0]=0;
        m[row][n-1]=0.; new[row][n-1]=0;
        for (col=1; col<n-1; ++col) {
            m[row][col]=0.;
        }
    }

    // Have only Node 0 set boundary conditions on top row of its own matrix
    if (rank==0) {
        for (col=0; col<n; ++col){
            m[0][col]=sin(M_PI*(double)(col)*dx);
            new[0][col]=m[0][col];
        }
    }

    // Have only last node set boundary conditions on bottom row of its own matrix
    if (rank==size-1){
        for (col=0; col<n; ++col){
            m[nrow+1][col]=sin(M_PI*(double)col*dx)*exp(-M_PI);
            new[nrow+1][col]=m[nrow+1][col];
        }

    }
    if (rank==0) a=1;
        else a=0;
    if (rank==(size-1)) b=-1;
        else b=0;

    // Main loop
    // Each node has its own matrix m holding "nrows" of rows (index ranges from 1 to nrow), plus 
    // one top row above and one bottom row below, totally "nrows + 2" rows. The top and bottom rows
    // (with row index 0 & nrow + 1) are used as buffer to store the rows exchanged with adjacent Nodes
    while (!conv) {
        norm=0.;

        // For each node, get the first non-buffer row of data (i.e. row 1) from next Node, and the 
        // bottom non-buffer row of data (i.e. row "nrow") from previous Node, and store them in "my" 
        // own buffer rows (i.e. rows 0 & nrow + 1)
        //
        // Make sure post receive before send
        if (rank<(size-1)) MPI_Irecv(&m[nrow+1][0],n,MPI_DOUBLE,rank+1,rank+1,MPI_COMM_WORLD,&req[0]);
        if (rank>0)        MPI_Irecv(&m[0][0],n,MPI_DOUBLE,rank-1,rank-1,MPI_COMM_WORLD,&req[1]);
        if (rank>0)        MPI_Isend(&m[1][0],n,MPI_DOUBLE,rank-1,rank,MPI_COMM_WORLD,&req[2]);
        if (rank<(size-1)) MPI_Isend(&m[nrow][0],n,MPI_DOUBLE,rank+1,rank,MPI_COMM_WORLD,&req[3]);

        // We only need to wait for the non-blocking MPI_Irecv calls above. Because the inner nodes
        // exchanged data with 2 neighbors, there are two wait handles. Similarly, the edge nodes
        // only has 1 wait handle.
        if ((rank<size-1)&&(rank>0)) MPI_Waitall(2,req,status);
        else if (rank==0)      MPI_Waitall(1,req,status);
        else if (rank==size-1) MPI_Waitall(1,&req[1],status);

        mynorm=0.;
        // Here we invoke OpenMP to utilize multiple threads in each MPI node
        #pragma omp parallel shared(m,new,nrow,n) private(ithread,nthread,row,col,tempnorm) reduction(+:mynorm)
        {
            tempnorm=0.;
            nthread=omp_get_num_threads();
            ithread=omp_get_thread_num();
            // printf("i %d n-thread %d \n",ithread,nthread);
            // Read from matrix m and store updated values in matrix "new"
            for (row=1+(nrow/nthread)*ithread; row<(nrow/nthread)*(ithread+1)+1;++row) {
                for (col=1; col<n-1; ++col) {
                    new[row][col] = .25*(m[row][col+1] + m[row+1][col] + m[row][col-1] + m[row-1][col]);
                    tempnorm+=pow(new[row][col]-m[row][col],2);
                    if (new[row][col]>100000.) printf("alarm!! %d %d %f %f\n", row,col,new[row][col],m[row][col]);
                }
            }
            mynorm+=tempnorm;
        }

        // Now we repeat the same steps but this time reading from matrix "new" and store updated 
        // values in matrix m. This way saves CPU cycles from having to copy data back to matrix m
        //
        // Again, start with exchanging rows with adjacent nodes, but this time using matrix "new"
        if (rank<(size-1)) MPI_Irecv(&new[nrow+1][0],n,MPI_DOUBLE,rank+1,rank+1,MPI_COMM_WORLD,&req2[0]);
        if (rank>0)        MPI_Irecv(&new[0][0],n,MPI_DOUBLE,rank-1,rank-1,MPI_COMM_WORLD,&req2[1]);
        if (rank>0)        MPI_Isend(&new[1][0],n,MPI_DOUBLE,rank-1,rank,MPI_COMM_WORLD,&req2[2]);
        if (rank<(size-1)) MPI_Isend(&new[nrow][0],n,MPI_DOUBLE,rank+1,rank,MPI_COMM_WORLD,&req2[3]);

        if ((rank<size-1)&&(rank>0)) MPI_Waitall(2,req2,status2);
        else if (rank==0)      MPI_Waitall(1,req2,status2);
        else if (rank==size-1) MPI_Waitall(1,&req2[1],status2);


        norm=0.;
        mynorm=0.;
        #pragma omp parallel shared(m,new,nrow,n) private(ithread,nthread,row,col,tempnorm) reduction(+:mynorm)
        {
            tempnorm=0;
            nthread=omp_get_num_threads();
            ithread=omp_get_thread_num();
            // Just as before, do Jacobi update here on matrix m instead
            for (row=1+(nrow/nthread)*ithread; row<(nrow/nthread)*(ithread+1)+1;++row) {
                for (col=1; col<n-1; ++col) {
                    m[row][col] = .25*(new[row][col+1] + new[row+1][col] + new[row][col-1] + new[row-1][col]);
                    if (m[row][col]>100000.) printf("alarm!! %d %d \n", row,col);   
                    tempnorm+=pow(new[row][col]-m[row][col],2);
                }
            }
            // Check convergence only on 2nd Jacobi update
            mynorm+=tempnorm;
        }
 
        // Increment iter by 2, because we did two Jacobi updates, one on matrix m and one on
        // matrix "new"

        iter=iter+2;
        MPI_Allreduce(&mynorm,&norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if (norm < 1e-10)   conv=1;
			
        if (rank==0) printf("rank %d norm %f iter %d\n",rank,norm,iter);
    } // End of while-loop


    // Print out results
    for (i=0;i<size;++i){
        if (rank==i){
            for (row=1; row<nrow+1; ++row) {
                printf("rank %d ",rank);
                for (col=0; col<n; ++col) {
                    printf("%f ",m[row][col]);
                }
		printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
		
    MPI_Finalize();
    return 0;
}

