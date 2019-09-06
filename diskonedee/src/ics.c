#include <stdio.h>
#include <math.h>


#include "disk.h"
#include "global.h"
#include "gsl_sf_erf.h"

/*
  Some simple initial condition profiles
*/

double get_profile_powerlaw_truncated(double R)
{

  return 0.01*pow(R,-0.5) * (gsl_sf_erf (R - 6.) + 1)* R;

}

double get_profile_powerlaw(double R)
{

  return 0.1*pow(R,-0.5) * R;

}


double get_profile_deltafunc(double R)
{

  return R * exp(-(R-20)*(R-20)/10.0);

}

double get_profile_similarity_cavity(double R)
{
  double sigma0 = 9.115351581034054e-05;
  double Rcav = 2.0;
  double xi = 4.0;
  double Rc = 12.0;
  double gamma = 0.5;
  double sigma = sigma0 * pow(R/Rc,-gamma)* exp(-pow(R/Rc,2-gamma)) * exp(-pow(Rcav/R,xi));
  if (sigma < 1.e-18)
      sigma = 1.e-18;	
  return R * sigma;

}
