#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include "disk.h"
#include "global.h"





double *invert_tridiagonal_problem(double *x, double * dd, double * du, double * dl, double * b)
{
  /*
    Simple routine to invert a general tridiagonal matrix A and solve for Ax = b
    Inputs: 
           -double *dd: (pointer to) vector of length M containing
	   diagonal elements of A
           -double *du: (pointer to) vector of length M containing 
	   upper diagonal elements of A 
           -double *dl: (pointer to) vector of length M containing 
	   lower diagonal elements of A
	   -double *b: (pointer to) vector on the RHS of equation
	   -double *x: (pointer to) vector to be solved for
   */

  x = (double *)malloc((M) *sizeof(double));
  /* First, determine the L, U matrices */
  double u[M], l[M];
  l[0] = 0;
  u[M-1] = 0;

  u[0] = dd[0];
  for (int k=1; k < M-1; k++)
    {
      l[k] = dl[k] / u[k-1];
      u[k] = dd[k] - l[k] / du[k-1];
    }

  /* Second, solve for y in Ly=b */
  double y[M];
  y[0] = b[0];
  for (int k=1; k < M; k++)
    y[k] = b[k] - l[k] * y[k-1];
  
  /* Finally, solve for x in Ux=y */
  x[M-1] = y[M-1] / u[M-1];
  for (int k = M-1; k>-1 ; k--)
    {
      x[k] = (y[k] - du[k] * x[k+1]) / u[k];
      //printf("x=%g\n",x[k]);
    }
  return x;
}


void invert_tridiagonal_CN(double alpha, double * b, double * x, int M)
{
  /*
    Simple routine to invert CRANK-NICOLSON type tridiagonal matrix of
    fixed computational diffusion coefficient alpha and solve for Ax = b
    Inputs: 
           -double alpha: scalar in the Crank-Nicolson scheme
	   -double *b: (pointer to) vector on the RHS of equation
	   -double *x: (pointer to) vector to be solved for
	   -int M: length of arrays
   */

  
  /* First, determine the L, U matrices */
  double u[M], l[M];
  l[0] = 0;
  u[M-1] = 0;

  u[0] = 2 * (1 + alpha);
  for (int k=1; k < M; k++)
    {
      l[k] = -alpha / u[k-1];
      u[k] = 2*(1 + alpha) + l[k] / alpha;
    }

  /* Second, solve for y in Ly=b */
  double y[M];
  y[0] = b[0];
  for (int k=1; k < M; k++)
    y[k] = b[k] - l[k] * y[k-1];
  
  /* Finally, solve for x in Ux=y */
  x[M-1] = y[M-1] / u[M-1];
  for (int k = M-2; k>-1 ; k--)
    x[k] = (y[k] + alpha * x[k+1]) / u[k];
}
