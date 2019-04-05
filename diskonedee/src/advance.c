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
  double b[M], dd[M], du[M], dl[M];

  for (int j = 1; j < M-1; j++)
    {
      b[j] = dt* 0.5 * G->vals[j].A_coeff_val* Q[j+1] + (1. + dt * 0.5 * G->vals[j].B_coeff_val)* Q[j] + dt * 0.5  * G->vals[j].C_coeff_val* Q[j-1];
      dd[j] = 1. - dt * 0.5 * G->vals[j].B_coeff_val;
      du[j] = -dt*0.5 * G->vals[j].A_coeff_val;
      dl[j] = -dt*0.5 * G->vals[j].C_coeff_val;
    }
  b[0] =  dt* 0.5 * G->vals[0].A_coeff_val* Q[1] + (1. + dt * 0.5 * G->vals[0].B_coeff_val)* Q[0]+
    dt * 0.5 * G->vals[1].C_coeff_val * Qin +  Qin * dt * 0.5 * G->vals[1].C_coeff_val;

  b[M-1] = (1. + dt * 0.5 * G->vals[M-1].B_coeff_val)* Q[M] + dt * 0.5  * G->vals[M-1].C_coeff_val* Q[M-1]+
    Qout * dt * 0.5 * G->vals[M-1].A_coeff_val + Qout * dt * 0.5 * G->vals[M-1].A_coeff_val;
  dl[0] = 0;
  dl[M-1]= -dt*0.5 * G->vals[M-1].C_coeff_val;
  du[0] = -dt*0.5 * G->vals[0].A_coeff_val;
  du[M-1] = 0;
  dd[0]= 1. - dt * 0.5 * G->vals[0].B_coeff_val;
  dd[M-1]= 1. - dt * 0.5 * G->vals[M-1].B_coeff_val;


  
  gsl_matrix *m = gsl_matrix_alloc(M, M);
  gsl_vector *x = gsl_vector_alloc(M);
  for (int j = 1; j < M-1; j++)
    {
      gsl_matrix_set(m,j,j, dd[j]);
      gsl_matrix_set(m,j,j+1, du[j]);
      gsl_matrix_set(m,j,j-1, dl[j]);
      gsl_vector_set(x, j, b[j]);
    }

  
  gsl_matrix_set(m,0,0,dd[0]);
  gsl_matrix_set(m,0,1,du[0]);
  gsl_vector_set(x, 0, b[0]);
  
  gsl_matrix_set(m,M-1,M-1,dd[M-1]);
  gsl_matrix_set(m,M-1,M-2,dl[M-1]);
  gsl_vector_set(x, M-1, b[M-1]);

  /*
  for (int j = 1; j < M-1; j++)
    for (int i = 0; i < M-1; i++) 
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
