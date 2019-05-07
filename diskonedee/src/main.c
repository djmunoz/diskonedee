#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "disk.h"
#include "global.h"



int main(int argc, char **argv)
{
  
  read_params();
  
  int Noutput=40;  
  double t= 0;

  printf("TimeMax=%g\n",params.TimeMax);
  printf("AlphaCoefficient=%g\n",params.AlphaCoefficient);
  printf("VerticalAspectRatio=%g\n", params.VerticalAspectRatio);
  printf("TempProfileIndex=%g\n",params.TempProfileIndex);
  
  
  N = params.TimeMax/dt;

  
  struct grid *Grid = init_grid(Grid);
  double *lambda=init_quant(lambda, Grid);
  init_grid_terms(lambda,Grid);
  init_boundaries(lambda,Grid);
  
  char *name;
  if (params.InnerBoundaryConditionType == 0)
    {
      name = "output_zerotorque.txt";
    }
  else if (params.InnerBoundaryConditionType == 1)
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
	write_to_file(lambda,Grid,output,t);
      printf("j=%d, t=%g\n",j,t);
      advance(lambda, Grid, dt);
      t += dt;

    }
  free(lambda);
  free(Grid);
  fclose(output);
}



