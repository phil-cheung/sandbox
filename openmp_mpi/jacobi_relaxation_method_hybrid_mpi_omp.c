// Make sure post receive before send to avoid delay
#include "mpi.h"
#include <omp.h>
#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[])
{

//	double M_PI=3.14159265358979323846;
	int size,rank,nrow=8,ithread,nthread;
	int i,j,conv,a,b,iter,row,col;
	double timei,timef,norm,mynorm,tempnorm,dx;
	MPI_Request req[4],req2[4];
	MPI_Status status[4],status2[4];

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	int n=nrow*size;	//matrix[n][n];
	double m[nrow+2][n],new[nrow+2][n];	//2 more rows for storage of data above and below

	printf("size is %d \n",size);
// initialize here
	conv=0;
	iter=0;
	dx=1./(double) (n-1) ;
	for (row=1;row<nrow+1;++row){
		m[row][0]=0.;new[row][0]=0;
		m[row][n-1]=0.;new[row][n-1]=0;
		for (col=1;col<n-1;++col) {
			m[row][col]=0.;
//			printf("%d %f ", rank,m[row][col]);
		}
	}
	if (rank==0) {
		for (col=0;col<n;++col){
			m[0][col]=sin(M_PI*(double)(col)*dx);
			new[0][col]=m[0][col];
//                        printf("\n top %d %f \n", rank,m[1][col]);
		}

	}
	if (rank==size-1){
		for (col=0;col<n;++col){
			m[nrow+1][col]=sin(M_PI*(double)col*dx)*exp(-M_PI);
			new[nrow+1][col]=m[nrow+1][col];
 //                       printf("\n bottom  %d %f \n ", rank,m[nrow][col]);
		}

	}
	if (rank==0) a=1;
	else a=0;
	if (rank==(size-1)) b=-1;
	else b=0;

//	printf("checkm  \n",m);
// main loop
	while (!conv) {
		norm=0.;

		if (rank<(size-1)) MPI_Irecv(&m[nrow+1][0],n,MPI_DOUBLE,rank+1,rank+1,MPI_COMM_WORLD,&req[0]);
		if (rank>0) MPI_Irecv(&m[0][0],n,MPI_DOUBLE,rank-1,rank-1,MPI_COMM_WORLD,&req[1]);
		if (rank>0) MPI_Isend(&m[1][0],n,MPI_DOUBLE,rank-1,rank,MPI_COMM_WORLD,&req[2]);
                if (rank<(size-1)) MPI_Isend(&m[nrow][0],n,MPI_DOUBLE,rank+1,rank,MPI_COMM_WORLD,&req[3]);

		if ((rank<size-1)&&(rank>0)) MPI_Waitall(2,req,status);
		else if (rank==0) MPI_Waitall(1,req,status);
		else if (rank==size-1) MPI_Waitall(1,&req[1],status);

		mynorm=0.;
		#pragma omp parallel shared(m,new,nrow,n) private(ithread,nthread,row,col,tempnorm) reduction(+:mynorm)
		{
		tempnorm=0.;
		nthread=omp_get_num_threads();
		ithread=omp_get_thread_num();
//		printf("i %d n-thread %d \n",ithread,nthread);
		for (row=1+(nrow/nthread)*ithread; row<(nrow/nthread)*(ithread+1)+1;++row) {
			for (col=1;col<n-1;++col) {
				new[row][col]=.25*(m[row][col+1]+m[row+1][col]+m[row][col-1]+m[row-1][col]);
//				printf("checkupdate %d %d %d %f  %f  ",rank,row,col, new[row][col],m[row][col]);
//				if (new[row][col]>100000.) printf("alarm!! %d %d %d\n", row,col,iter);
				tempnorm+=pow(new[row][col]-m[row][col],2);
                              if (new[row][col]>100000.) printf("alarm!! %d %d %f %f\n", row,col,new[row][col],m[row][col]);
				}
			}
			mynorm+=tempnorm;
		}
//                printf ("rank %d mynorm %f \n",rank,mynorm);

//		MPI_Barrier(MPI_COMM_WORLD);		//implicit barrier from parallel region
                if (rank<(size-1)) MPI_Irecv(&new[nrow+1][0],n,MPI_DOUBLE,rank+1,rank+1,MPI_COMM_WORLD,&req2[0]);
                if (rank>0) MPI_Irecv(&new[0][0],n,MPI_DOUBLE,rank-1,rank-1,MPI_COMM_WORLD,&req2[1]);
                if (rank>0) MPI_Isend(&new[1][0],n,MPI_DOUBLE,rank-1,rank,MPI_COMM_WORLD,&req2[2]);
                if (rank<(size-1)) MPI_Isend(&new[nrow][0],n,MPI_DOUBLE,rank+1,rank,MPI_COMM_WORLD,&req2[3]);

                if ((rank<size-1)&&(rank>0)) MPI_Waitall(2,req2,status2);
                else if (rank==0) MPI_Waitall(1,req2,status2);
                else if (rank==size-1) MPI_Waitall(1,&req2[1],status2);


		norm=0.;
                mynorm=0.;
                #pragma omp parallel shared(m,new,nrow,n) private(ithread,nthread,row,col,tempnorm) reduction(+:mynorm)
		{
		tempnorm=0;
                nthread=omp_get_num_threads();
                ithread=omp_get_thread_num();
                for (row=1+(nrow/nthread)*ithread; row<(nrow/nthread)*(ithread+1)+1;++row) {
//                        m[row][0]=new[row][0]; m[row][n-1]=new[row][n-1];
			for (col=1;col<n-1;++col) {
				m[row][col]=.25*(new[row][col+1]+new[row+1][col]+new[row][col-1]+new[row-1][col]);
                                if (m[row][col]>100000.) printf("alarm!! %d %d \n", row,col);   
  //                              printf("checkupdate2 %d %d %d %f  %f  ",rank,row,col, new[row][col],m[row][col]);
				tempnorm+=pow(new[row][col]-m[row][col],2);
				}
			}
			mynorm+=tempnorm;
		}
//		MPI_Barrier(MPI_COMM_WORLD);
		iter=iter+2;
//		printf ("rank %d mynorm %f \n",rank,mynorm);
		MPI_Allreduce(&mynorm,&norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		if (norm < 1e-10)   conv=1;
			
		if (rank==0) printf("rank %d norm %f iter %d\n",rank,norm,iter);

	}

	for (i=0;i<size;++i){
		if (rank==i){
		for (row=1; row<nrow+1;++row) {
			printf("rank %d ",rank);
			for (col=0;col<n;++col) {
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


//		for (row=rank*nrow+a+1; row<(rank+1)*nrow+b-1;++row) {
//			for (col=1;col<n-1;++col) {
//				new[row][col]=.25*(m[row][col+1]+m[row+1][col]+m[row+1][col+2]+m[row+2][col+1]);
//				printf("test %d %f    ",ithread, new[row][col]);
//				mynorm+=pow(new[row][col]-m[row+1][col+1],2);
//				}
//			}
