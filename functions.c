#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "../headers/functions.h"
#include "functions.h"

double vectNorm(double * v, int size)
{
	double norm= 0.;
	double sum = 0.;
	int i;

	for (i = 0; i < size; i++)
		sum += (v[i] * v[i]);
	
	norm = sqrt(sum);

	return norm;
}

void symMatInit(double ** A, int n)
{
	int i,j;
	
	for(i=0;i<n;i++)
		for(j=i;j<n;j++)
		{
			A[i][j]=i*j+2*n;
			A[j][i]=A[i][j];
		}
}

void matPrint(const char * string, double ** A, int n)
{
	int i,j;
	printf("\nMatrix name %s\n",string);
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			printf("%lf ",A[i][j]);
	
		printf("\n");
	}
}

void vectPrint(const char * string, double * v, int n)
{
	int i;
	printf("\nVector name %s\n",string);
	for( i = 0; i < n; i++)
		printf("%lf ",v[i]);
	printf("\n");
}

double * matToVect(double ** A, int n)
{
	double * res = (double *)malloc(sizeof(double)*n*n);
	int i,j;

	for( i = 0; i < n; i++)
		for( j = 0; j < n; j++)
			res[(i*n)+j] = A[i][j];
	
	return res;
}

double dotProd(double *v, double *w, int t)
{
	double res = 0.;
	int i;

	for( i = 0; i < t; i++)
		res += v[i] * w[i];
	
	return res;
}

double * vectProd(double *v, double *w, int t)
{
	double * res = (double *)malloc(sizeof(double)*t);
	int i;

	for( i = 0; i < t; i++)
		res[i] = v[i] * w[i];
	
	return res;
}

double * matVectProd( double ** A, int nbl, int nbc, double * v, int t)
{
	double * res = calloc ( nbl, sizeof( double));
	double somme = 0.;
	int i,j;

	for( i = 0; i < nbl; i++)
	{
		somme = 0.;
		
		for( j = 0; j < nbc; j++)
			somme += v[j] * A[i][j];
			
		res[i] = somme;
	}
	
	return res;
}
