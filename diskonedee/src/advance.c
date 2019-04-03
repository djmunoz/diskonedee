#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "disk.h"
#include "global.h"




void advance(double * Q, double * G, double dt)
{

  advance_euler(Q, G, dt);
  
  return;

}


void advance_euler(double * Q, double * G, double dt)
{
  Q[0] = 0;
  Q[M-1] = 0;
  for (int k = 1; k < M-1; k++)
    {
       
      Q[k] += -dt * (G->vals[j+1].A_coeff_val * Q[k-1] + G->vals[j+1].B_coeff_val * Q[k] + G->vals[j+1].C_coeff_val * Q[k+1];

    }
  
  return;

}
