#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include "disk.h"


static int M = 50;
static int N = 100;
static double timemax = 5.0;
  
void main()
{

  double t, dt;
  dt = timemax/N;
  
  double * W;
  init(W);
  
  /* Main loop */
  for (int j=0; j < N; j++) 
    {
      t += dt * j;
      //advance(W, dt);
      //write_to_file(W,t);
    }
  free(W);

}


void init(double * W)
{
  //read_from_file()
  W = (double *)malloc(M *sizeof(double));
}
