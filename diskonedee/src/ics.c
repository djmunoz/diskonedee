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

  return 0.01*pow(R,-0.5) * (gsl_sf_erf (R - 3.) + 1)* R;

}

double get_profile_powerlaw(double R)
{

  return 0.1*pow(R,-0.5) * R;

}


double get_profile_deltafunc(double R)
{

  return exp(-(R-30)*(R-30)/20.0);

}
