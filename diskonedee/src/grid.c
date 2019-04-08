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

  double Rc, Rm, Rp, lprimem, lprimep;
  for (int j=0; j < M; j++)
    {
      if (j == M-1)
	delta_inv = 1.0 / (G->vals[j].center_val - G->vals[j-1].center_val);
      else
	delta_inv = 1.0 / (G->vals[j+1].side_val - G->vals[j].side_val);
      if ((j > 0) && (j < M-1))
	{
	  Rc = G->vals[j].center_val;
	  Rp = G->vals[j+1].center_val;
	  Rm = G->vals[j-1].center_val;
	  delta_inv = 1.0 / (G->vals[j+1].side_val - G->vals[j].side_val);
	  lprimep = eval_lprime_func(G->vals[j+1].side_val);
	  lprimem = eval_lprime_func(G->vals[j].side_val);
	}
      else if (j == 0)
	{
	  Rc = G->vals[0].center_val;
	  Rp = G->vals[1].side_val;
	  Rm = G->vals[0].side_val;
	  delta_inv = 1.0 / (G->vals[1].side_val - Rmin);
	  lprimep = eval_lprime_func(G->vals[1].side_val);
	  lprimem = eval_lprime_func(Rmin);
	}
      else
	{
	  Rc = G->vals[M-1].center_val;
	  Rp = G->vals[M-1].side_val;
	  Rm = G->vals[M-2].side_val;
	  delta_inv = 1.0 / (Rmax - G->vals[M-1].side_val);
	  lprimep = eval_lprime_func(Rmax);
	  lprimem = eval_lprime_func(G->vals[M-1].side_val);
	}
      
      func_jplus = eval_nu_func(Rp)  * eval_omegaprime_func(Rp) * Rp * Rp;
      func_j = eval_nu_func(Rc)  * eval_omegaprime_func(Rc) * Rc * Rc;
      func_jminus = eval_nu_func(Rm)  * eval_omegaprime_func(Rm) * Rm * Rm;
      

      G->vals[j].A_coeff_val = -1.0/lprimep * delta_inv * func_jplus /(Rp - Rc);
      G->vals[j].B_coeff_val = (1.0/lprimep/(Rp - Rc) +  1.0/lprimem/(Rc - Rm)) * func_j * delta_inv;
      G->vals[j].C_coeff_val =  -1.0/lprimem * delta_inv * func_jminus /(Rc - Rm);

      G->vals[j].cn_upper_diag = -dt*0.5 * G->vals[j].A_coeff_val;
      G->vals[j].cn_middle_diag = 1. - dt * 0.5 * G->vals[j].B_coeff_val;
      G->vals[j].cn_lower_diag = -dt*0.5 * G->vals[j].C_coeff_val;
    }
  
  if (params.BoundaryConditionType == 1)
    {
      Rc = G->vals[0].center_val;
      Rp = G->vals[1].center_val;
      G->vals[0].cn_middle_diag = 1 - 0.5 / (G->vals[1].side_val - Rmin) * eval_nu_func(Rc) / Rc *
	(Rc * Rc * Rc * eval_omegaprime_func(Rc) / eval_lprime_func(G->vals[0].side_val) / (Rp - Rc) +
	 Rmin * Rmin * Rmin * eval_omegaprimeprime_func(Rmin) / eval_lprime_func(Rmin) +
	 3 * Rmin * Rmin * eval_omegaprime_func(Rmin) / eval_lprime_func(Rmin));
      G->vals[0].cn_upper_diag = 0.5 / (G->vals[1].side_val - Rmin) * eval_nu_func(Rp) / Rp *
	(Rp * Rp * Rp * eval_omegaprime_func(Rp) / eval_lprime_func(G->vals[0].side_val) / (Rp - Rc));

      G->vals[0].cn_lower_diag = 0;

      G->vals[M-1].cn_upper_diag = 0;
      G->vals[M-1].cn_middle_diag = 0.25 - 0.5 * eval_g_func(G->vals[M-1].side_val) / (G->vals[M-1].center_val - G->vals[M-2].center_val);
      G->vals[M-1].cn_lower_diag = 0.25 + 0.5 * eval_g_func(G->vals[M-1].side_val) / (G->vals[M-1].center_val - G->vals[M-2].center_val);
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



