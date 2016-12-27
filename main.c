#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <lapacke.h>
//#include "../headers/functions.h"
#include "functions.h"

int main(int argc, char ** argv)
{
	int m, n;
	int i,j,k=0;

	if(argc != 3)
	{
		perror("USAGE: ./prr <m> <Matrix size>");
	}else{
		m = atoi(argv[1]);
		n = atoi(argv[2]);
	}

	double * v = (double *)malloc(sizeof(double)*n);
	double * y = (double *)malloc(sizeof(double)*n);
	double * b = (double *)malloc(sizeof(double)*m);
	double ** A = (double **)malloc(sizeof(double*)*n);
	double ** V = (double **)malloc(sizeof(double*)*n);
	double ** C = (double **)malloc(sizeof(double*)*m);
	double * c = (double *)malloc(sizeof(double)*m*m);
	double * lap = (double *)malloc(sizeof(double)*m*m);
	double lambda1,lambda2,delta;

	for(i=0;i<n;i++)
	{
		A[i] = (double *)malloc(sizeof(double)*n);		
		V[i] = (double *)malloc(sizeof(double)*n);		
	}
	for(i=0;i<m;i++)
		C[i] = (double *)malloc(sizeof(double)*m);		

	//symMatInit(A,n);
	//matPrint("A",A,n);
	A[0][0]=9.;A[0][1]=1.;A[0][2]=-2.;A[0][3]=1.;A[1][1]=8.;A[1][2]=-3.;A[1][3]=-2.;A[2][2]=7.;A[2][3]=-1.;A[3][3]=6.;
	for(i=0;i<n;i++)
		for(j=i;j<n;j++)
			A[j][i]=A[i][j];
	matPrint("A",A,n);
	v[0]=1;
	for(i=1;i<n;i++)
		v[i]=0.0;
	c[k] = 1;
	y=matVectProd(A,n,n,v,n);		

	for(i=0;i<n;i++)	
	{
		V[i][0]=v[i];
		V[i][1]=y[i];	
	}

	for(j=0;j<m-1;j++)
	{
		for(i=0;i<n;i++)	
		{	
			v[i]=V[i][j];
			y[i]=V[i][j+1];	
		}

		k++;	
		c[k]=dotProd(v,y,n);
		k++;
		c[k]=dotProd(y,y,n);
		y=matVectProd(A,n,n,y,n);

		for(i=0;i<n;i++)
			V[i][j+2]=y[i];
	}
		
	for(i=0;i<n;i++)
		v[i]=V[i][k];
	y=matVectProd(A,n,n,v,n);
	for(i=0;i<n;i++)
		V[i][n-1]=y[i];
	//vectPrint("v",v,n);
	//vectPrint("y",y,n);
	//matPrint("V",V, n);
	
	for(i=0;i<n;i++)	
	{	
		v[i]=V[i][n-2];
		y[i]=V[i][n-3];	
	}
	k++;
	c[k]=dotProd(v,y,n);
	//printf("%lf\n",c[k]);	
	b[m-1]=c[k];			
	for(i=0;i<m;i++)
		for(j=0;j<m;j++)
		{
			C[i][j]=c[j+i];
			C[j][i]=C[i][j];
		}
	/*printf("\n");
	matPrint("C",C,m);*/

	j=m;
	for(i=0;i<m;i++)
	{
		b[i]=-c[j];
		j++;
	}
	//matPrint("C",C,m);
	lap=matToVect(C,m);
	vectPrint("b",b,m);
	lapack_int ipiv[m];
	lapack_int info = LAPACKE_dgesv( LAPACK_COL_MAJOR, m, m, lap, m, ipiv, b, m );	
	vectPrint("b",b,m);
	return 0;
}
