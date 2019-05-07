#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>


#include "disk.h"
#include "global.h"


double eval_omega_func(double R)
{
  double GM = 1.0;
  double omega;
  if (params.Softening == 0)
      omega = sqrt(GM) / R / sqrt(R);
  else
    {
      double h = params.Softening * 2.8;
      double u, fac, h_inv, h3_inv;
      /*Use spline softening*/
      if(R >= h)
	fac = 1 / (R * R * R);
      else
	{
	  h_inv = 1.0 / h;
	  h3_inv = h_inv * h_inv * h_inv;
	  u = R * h_inv;
	  if(u < 0.5)
	    fac = h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac = h3_inv * (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
	}
      
      omega = sqrt(GM * fac);
    }
  
  return omega;
}

double eval_omegaprime_func(double R)
{
  double GM = 1.0;
  double omegaprime;
  if (params.Softening == 0)
      omegaprime  = -1.5 * sqrt(GM) / R / R / sqrt(R); 
  else
    {
      double h = params.Softening * 2.8;
      double u, fac, facprime, h_inv, h3_inv, h4_inv;
      /*Use spline softening*/
      if(R >= h)
	{
	  fac = 1 / (R * R * R);
	  facprime = -3. / (R * R * R * R);
	}
      else
	{
	  h_inv = 1.0 / h;
	  h3_inv = h_inv * h_inv * h_inv;
	  h4_inv = h3_inv * h_inv;
	  u = R * h_inv;
	  if(u < 0.5)
	    {
	      fac = h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	      facprime = -h4_inv * (76.8 * u - 96. * u * u);
	    }
	  else
	    {
	      fac = h3_inv * (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
	      facprime = -h4_inv * (-76.8 * u - 0.2 / u / u / u / u + 48.  + 32. * u * u);
	    }
	}

      omegaprime = 0.5 * sqrt(GM / fac) * facprime;
    }

  
  return omegaprime;
}

double eval_omegaprimeprime_func(double R)
{
  double GM = 1.0;
  double omegaprimeprime = 3.75 * sqrt(GM) / R / R / R / sqrt(R);
  return omegaprimeprime;

}

double eval_h_func(double R)
{
  double h = params.VerticalAspectRatio * pow(R, 0.5 - 0.5 * params.TempProfileIndex);
  return h;

}
double eval_nu_func(double R)
{
  double GM = 1.0;
  double nu;
  if (params.Softening == 0)
      nu = params.AlphaCoefficient * params.VerticalAspectRatio * params.VerticalAspectRatio * sqrt(GM) * pow(R, 1.5 - params.TempProfileIndex);
  else
    nu = params.AlphaCoefficient * params.VerticalAspectRatio * params.VerticalAspectRatio * pow(R, 3 - params.TempProfileIndex) * eval_omega_func(R);
  
  return nu;

}

double eval_nuprime_func(double R)
{
  double GM = 1.0;
  double nuprime;
  if (params.Softening == 0)
    nuprime = (1.5 - params.TempProfileIndex) * params.AlphaCoefficient * params.VerticalAspectRatio * params.VerticalAspectRatio * sqrt(GM) * pow(R, 0.5 - params.TempProfileIndex);
  else
    nuprime =  (3.0 - params.TempProfileIndex) * params.AlphaCoefficient * params.VerticalAspectRatio * params.VerticalAspectRatio * \
      pow(R, 1 - params.TempProfileIndex) * eval_omega_func(R) +\
      params.AlphaCoefficient * params.VerticalAspectRatio * params.VerticalAspectRatio * pow(R, 3 - params.TempProfileIndex) * eval_omegaprime_func(R);

  return nuprime;

}

double eval_l_func(double R)
{
  double GM = 1.0;
  double l;
  if (params.Softening == 0)
        l = sqrt(GM) * sqrt(R);
  else
    l = R * R * eval_omega_func(R);
  
  return l;

}

double eval_lprime_func(double R)
{
  double GM = 1.0;
  double lprime;
  if (params.Softening == 0)
    lprime = 0.5 * sqrt(GM) / sqrt(R);
  else
    lprime =  R * (2 * eval_omega_func(R) + R * eval_omegaprime_func(R));
    
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
  R0 = 2.8;
  w = 0.4;
  Omega= 1./sqrt(R0)/R0;
  return Omega*exp(-(R-R0)*(R-R0)/2/w/w);
}

double eval_beta_func(double R)
{
  return 0; //blip(R);
}


double eval_gamma_func(double R)
{
  return 0.005 * blip(R);
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
