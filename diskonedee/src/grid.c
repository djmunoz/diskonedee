#include <stdio.h>
#include <math.h>


#include "disk.h"
#include "global.h"

/*
struct grid
{
  int Mpoints;
  double vals[];
  //double Alpha[];
  //double Beta[];
};
*/
struct grid *initGrid(struct grid *G, double R1, double R2, int Npoints)
{
  G = malloc( sizeof(*G) +  (M ) * sizeof(struct gridvals));
  
  G->Npoints = Npoints;
  double deltalog = (log10(Rmax) - log10(Rmin))/Npoints;
  double delta_inv, func_jplus, func_jminus, func_j; 
  for (int j=0; j < M; j++)
    {
      G->vals[j].center_val = Rmin * pow(10, (j*1. + 0.5) * deltalog);
      G->vals[j].side_val = Rmin * pow(10, j * deltalog);
    }
  
  for (int j=1; j < M-1; j++)
    {
      delta_inv = 1.0 / (G->vals[j+1].side_val - G->vals[j].side_val);
      func_jplus = eval_nu_func(G->vals[j+1].center_val)  * eval_omegaprime_func(G->vals[j+1].center_val) * G->vals[j+1].center_val * G->vals[j+1].center_val;
      func_j = eval_nu_func(G->vals[j].center_val)  * eval_omegaprime_func(G->vals[j].center_val) * G->vals[j].center_val * G->vals[j].center_val;
      func_jminus = eval_nu_func(G->vals[j-1].center_val)  * eval_omegaprime_func(G->vals[j-1].center_val) * G->vals[j-1].center_val * G->vals[j-1].center_val;

      

      G->vals[j].A_coeff_val = -1.0/eval_lprime_func(G->vals[j+1].side_val) * delta_inv * func_jplus /(G->vals[j+1].center_val - G->vals[j].center_val);
      G->vals[j].B_coeff_val = (1.0/eval_lprime_func(G->vals[j+1].side_val)/(G->vals[j+1].center_val - G->vals[j].center_val) +  1.0/eval_lprime_func(G->vals[j].side_val)/(G->vals[j].center_val - G->vals[j-1].center_val)) * func_j * delta_inv;
      G->vals[j].C_coeff_val =  -1.0/eval_lprime_func(G->vals[j].side_val) * delta_inv * func_jminus /(G->vals[j].center_val - G->vals[j-1].center_val);

      G->vals[j].cn_upper_diag = -dt*0.5 * G->vals[j].A_coeff_val;
      G->vals[j].cn_middle_diag = 1. - dt * 0.5 * G->vals[j].B_coeff_val;
      G->vals[j].cn_lower_diag = -dt*0.5 * G->vals[j].C_coeff_val;
    }

  if (params.BoundaryConditionType == 0)
    {
      G->vals[0].cn_upper_diag =  -dt*0.5 * (-1.0/eval_lprime_func(G->vals[1].side_val) * delta_inv * func_jplus /(G->vals[1].center_val - G->vals[0].center_val););
      G->vals[0].cn_middle_diag = 1. - dt * 0.5 * ((1.0/eval_lprime_func(G->vals[1].side_val)/(G->vals[j+1].center_val - G->vals[j].center_val) +  1.0/eval_lprime_func(G->vals[0].side_val)/(G->vals[j].center_val - G->vals[j-1].center_val)) * func_j * delta_inv);
	G->vals[0].cn_lower_diag = 0;

    }
 
   return G;
}

double *init_quant(double *Q, struct grid *G)
{
  //read_from_file()
  Q = (double *)malloc((M) *sizeof(double));
  double radius;
  for (int j=0; j < M; j++) 
    {
      radius = G->vals[j].center_val;
      Q[j] = exp(-(radius-30)*(radius-30)/20.0);
    }

  return Q;
}



