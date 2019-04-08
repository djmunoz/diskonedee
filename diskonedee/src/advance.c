#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>


#include "disk.h"
#include "global.h"




void advance(double * Q,  struct grid *G, double dt)
{

  //advance_euler(Q, G, dt);

  advance_cranknicolson(Q, G, dt);
  
  
  return;

}


void advance_euler(double * Q, struct grid *G, double dt)
{
  Q[0] = eval_inner_boundary(Q,G);
  Q[M+1] = eval_outer_boundary(Q,G);
  double Qjminus = Q[0];
  double Qj;
  for (int j = 1; j < M+1; j++)
    {
      Qj= Q[j] + dt * (G->vals[j].A_coeff_val * Q[j+1] + G->vals[j].B_coeff_val * Q[j] + G->vals[j].C_coeff_val * Qjminus);
      Qjminus = Q[j];
      Q[j] = Qj;
      //Q[j] += dt * (G->vals[j].A_coeff_val * Q[j+1] + G->vals[j].B_coeff_val * Q[j] + G->vals[j].C_coeff_val * Q[j-1]);
            
    }
  
  
  return;

}

void advance_cranknicolson(double * Q, struct grid *G, double dt)
{
  double Qin = 0; //eval_inner_boundary(Q,G);
  double Qout = 0; // eval_outer_boundary(Q,G);
  /* Fill the RHS vector and diagonals */
  double b[M];

  b[0] = -G->vals[0].cn_upper_diag * Q[1] + (2.- G->vals[0].cn_middle_diag)* Q[0];
  b[M-1] = (2.- G->vals[M-1].cn_middle_diag)* Q[M-1] -G->vals[M-1].cn_lower_diag * Q[M-2];

  printf("b[0]=%g zeta0=%g, eta0=%g\n",b[0],(2.- G->vals[0].cn_middle_diag),G->vals[M-1].cn_middle_diag);
  
  gsl_matrix *m = gsl_matrix_alloc(M, M);
  gsl_vector *x = gsl_vector_alloc(M);
  for (int j = 1; j < M-1; j++)
    {
      gsl_matrix_set(m,j,j,  G->vals[j].cn_middle_diag);
      gsl_matrix_set(m,j,j+1, G->vals[j].cn_upper_diag);
      gsl_matrix_set(m,j,j-1, G->vals[j].cn_lower_diag);
      b[j] = -G->vals[j].cn_upper_diag * Q[j+1] + (2.- G->vals[j].cn_middle_diag)* Q[j] -G->vals[j].cn_lower_diag * Q[j-1];
      gsl_vector_set(x, j, b[j]);
    }
  gsl_matrix_set(m,0,0,G->vals[0].cn_middle_diag);
  gsl_matrix_set(m,0,1,G->vals[0].cn_upper_diag);
  gsl_vector_set(x, 0, b[0]);
  gsl_matrix_set(m,M-1,M-1,G->vals[M-1].cn_middle_diag);
  gsl_matrix_set(m,M-1,M-2,G->vals[M-1].cn_lower_diag);
  gsl_vector_set(x, M-1, b[M-1]);


  /*
  for (int j = 0; j < M; j++)
    for (int i = 0; i < M; i++) 
	printf ("m(%d,%d) = %g\n", i, j, 
		gsl_matrix_get (m, i, j));
  */
  gsl_vector *Qnew = gsl_vector_alloc(M);


  gsl_matrix *inverse = invert_a_matrix(m);

  
  gsl_blas_dgemv(CblasNoTrans, 1.0, inverse, x, 0.0, Qnew);
  
  //double *Qnew = invert_tridiagonal_problem(Qnew, dd, du, dl, b);
  //Q = Qnew;
  //invert_tridiagonal_problem(Q, dd, du, dl, b);
  for (int j = 0; j < M; j++)
      Q[j] = gsl_vector_get(Qnew,j);

  //invert_tridiagonal_problem(Q, dd, du, dl, b);  
  return;

}
