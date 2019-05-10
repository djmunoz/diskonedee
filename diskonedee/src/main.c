#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "disk.h"
#include "global.h"



int main(int argc, char **argv)
{
  
  read_params();
  
  int Noutput=50;  
  double t= 0;

  N = params.TimeMax/params.dt;
  dt = params.dt;
  M = params.Ngrid;
  
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



