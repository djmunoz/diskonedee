#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "disk.h"
#include "global.h"



int main(int argc, char **argv)
{
  int Noutput=30;  
  double t= 0;
  N = timemax/dt;
  
  params.AlphaCoefficient = 0.06;
  params.VerticalAspectRatio = 0.1;
  params.TempProfileIndex = 1.0;
  params.BoundaryConditionType = 1;
    
  struct grid *Grid = initGrid(Grid, Rmin, Rmax, M);
  double *lambda=init_quant(lambda, Grid);

  char *name;
  if (params.BoundaryConditionType == 0)
    {
      name = "output_zerotorque.txt";
    }
  else if (params.BoundaryConditionType == 1)
    {
      name = "output_zamfb.txt";
    }
  else
    {
      name = "output.txt";
    }
  
  FILE *output = fopen(name, "w");
  write_header(Grid,output);
  
  /* Main loop */
  for (int j=0; j < N; j++) 
    {
      if (j % (N/Noutput) == 0) 
	write_to_file(lambda,Grid,output);
      printf("j=%d, t=%g\n",j,t);
      advance(lambda, Grid, dt);
      t += dt;

    }
  free(lambda);
  free(Grid);
  fclose(output);
}



