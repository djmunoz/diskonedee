#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

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
  Q[0] = eval_inner_boundary(Q,G);
  Q[M+1] = eval_outer_boundary(Q,G);
  /* Fill the RHS vector and diagonals */
  double b[M], dd[M], du[M], dl[M];

  for (int j = 1; j < M+1; j++)
    {
      b[j] = dt*(0.5 * G->vals[j].A_coeff_val* Q[j+1] + (-1. + 0.5 * G->vals[j].B_coeff_val)* Q[j] + 0.5  * G->vals[j].C_coeff_val* Q[j-1]);
      dd[j] = dt*(1. - 0.5 * G->vals[j].B_coeff_val);
      du[j] = -dt*0.5 * G->vals[j].A_coeff_val;
      dl[j] = -dt*0.5 * G->vals[j].C_coeff_val;
      printf("b=%g\n",b[j]);
    }
  
  double *Qnew = invert_tridiagonal_problem(Qnew, dd, du, dl, b);
  for (int j = 1; j < M+1; j++)
    Q[j] = Qnew[j-1];
  
  //invert_tridiagonal_problem(Q, dd, du, dl, b);  
  return;

}
