#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "disk.h"
#include "global.h"



int main(int argc, char **argv)
{
    
  double t= 0;
  double dt;
  dt = timemax/N;

  struct grid *Grid = initGrid(Grid, Rmin, Rmax, M);

  double * lambda;
  init_quant(lambda);//, Grid);
  
  
  /* Main loop */
  for (int j=0; j < N; j++) 
    {
      t += dt;
      printf("j=%d, t=%g\n",j,t);
      //advance(lambda, Grid, dt);
      //write_to_file(W,t);
    }
  free(lambda);
  free(Grid);
}



