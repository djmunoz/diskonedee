#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>


#include "disk.h"
#include "global.h"


double eval_omegaprime_func(double R)
{
  double GM = 1.0;
  double omegaprime = -1.5 * sqrt(GM) / R / R / sqrt(R);
  return omegaprime;

}

double eval_omegaprimeprime_func(double R)
{
  double GM = 1.0;
  double omegaprimeprime = 3.75 * sqrt(GM) / R / R / R / sqrt(R);
  return omegaprimeprime;

}

double eval_nu_func(double R)
{
  double GM = 1.0;
  double nu = params.AlphaCoefficient * params.VerticalAspectRatio * params.VerticalAspectRatio * sqrt(GM) * pow(R, 1.5 - params.TempProfileIndex);
  return nu;

}

double eval_nuprime_func(double R)
{
  double GM = 1.0;
  double nuprime = (1.5 - params.TempProfileIndex) * params.AlphaCoefficient * params.VerticalAspectRatio * params.VerticalAspectRatio * sqrt(GM) * pow(R, 0.5 - params.TempProfileIndex);
  return nuprime;

}

double eval_l_func(double R)
{
  double GM = 1.0;
  double l = sqrt(GM) * sqrt(R);
  return l;

}

double eval_lprime_func(double R)
{
  double GM = 1.0;
  double lprime = 0.5 * sqrt(GM) / sqrt(R);
  return lprime;

}

double eval_g_func(double R)
{
  double g = -1.0/(eval_lprime_func(R) / eval_l_func(R) -
		   eval_nuprime_func(R)/ eval_nu_func(R) -
		   2./ R -
		   eval_omegaprimeprime_func(R) / eval_omegaprime_func(R));
  return g;

}

double blip(double R)
{
  double R0,w, Omega;
  R0 = 5.0;
  w = 0.5;
  Omega= 1./sqrt(R0)/R0;
  return Omega*exp(-(R-R0)*(R-R0)/2/w/w);
}

double eval_beta_func(double R)
{
  return 0; //blip(R);
}


double eval_gamma_func(double R)
{
  return 0.01 * blip(R);
}



double eval_inner_boundary(double * Q, struct grid *G)
{
  double radius = Rmin; //G->vals[0].center_val;
  double coeff = (eval_lprime_func(radius) / eval_l_func(radius) -
		  (eval_nuprime_func(radius) * radius * radius * eval_omegaprime_func(radius) +
		   eval_nu_func(radius) * 2 * radius * eval_omegaprime_func(radius) +
		   eval_nu_func(radius) * radius * radius * eval_omegaprimeprime_func(radius))/
		  (eval_nu_func(radius) * radius * radius * eval_omegaprime_func(radius)));
  double Q0 = (Q[2] - Q[1]) / (G->vals[2].center_val - G->vals[1].center_val) / coeff;
  
  Q0=0;
  return Q0;
}


double eval_outer_boundary(double * Q, struct grid *G)
{
  
  double Q1 = Q[M-2] + (Q[M-2] - Q[M-3]) / (G->vals[M-2].center_val - G->vals[M-3].center_val);
  //Q1=0;
  return Q1;
}
