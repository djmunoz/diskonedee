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
  G = malloc( sizeof(*G) + 3 * M *sizeof(double));

  G->Npoints = Npoints;
  double deltalog = log(Rmax) - log(Rmin);
  double centervals, halfvals;
  for (int k=0; k < M; k++)
    {
      centervals = Rmin * pow(10, k * deltalog);
      halfvals = Rmin * pow(10, 0.5 * k * deltalog);
      G->vals[k].radius_val = centervals;
      G->vals[k].alpha_val = 1.0;
      G->vals[k].beta_val = 1.0;
    }

  return G;
}

void init_quant(double *Q)//, struct grid *G)
{
  //read_from_file()
  Q = (double *)malloc(M *sizeof(double));
  for (int k=0; k < M; k++) 
    {
      Q[k] = exp(-1);
    }
}



