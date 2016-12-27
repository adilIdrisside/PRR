#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double vectNorm(double * v, int size); /* Norme du vecteur */

void symMatInit(double ** A, int n); /* Matrix initialisation */

void matPrint(const char * string, double ** A, int n); /* Matrix printing */

void vectPrint(const char * string, double * v, int n); /* Vector printing */

double * matToVect(double ** A, int n); /* Matrix linearisation */

double dotProd(double *v, double *w, int t); /* Produit scalaire */

double * vectProd(double *v, double *w, int t); /* Produit vectoriel */

double * matVectProd( double ** A, int nbl, int nbc, double * v, int t); /* Produit matrice-vecteur */

#endif //FUNCTIONS_H_INCLUDED
