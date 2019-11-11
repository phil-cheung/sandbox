// Simple code to compare performance of CBLAS or MLK library

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cblas.h"
//#include "mkl.h"

// time_t time(time_t *time);

int main(int argc,char *argv[])
{
    time_t timef, timei;
    int n=600000;
    int m=1000000;
    double a[n],b[n];         //consists of vector of 1,2,..,40 components
    double avgtime,dot;
    int done = 0, myid, numprocs, i,j,k,index;
    
    for (i=0;i<=n-1;i+=1)
      {
	a[i]=20.2;
	b[i]=10.23;
      }

    for (i=1;i<=n;i*=2)
      {
	         timei=time(NULL);
		 //		 scanf("%d",&dot);
	    k=10000000000/i;
	    for (j=0;j<k;j+=1)
	    {
	      	 	dot=cblas_ddot(i,&a[1],1,&b[1],1);
	     }
	timef=time(NULL);
	printf("It takes %e seconds to do a %d vector dotproduct\n",difftime(timef,timei), i+1);
      }
 
    return 0;
}
