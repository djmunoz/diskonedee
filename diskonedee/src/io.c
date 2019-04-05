#include <stdio.h>
#include <time.h>
#include <math.h>

#include "disk.h"
#include "global.h"




void write_header(struct grid *G, FILE *output)
{
  for (int j=1; j < M; j++)
    fprintf(output, "%3g ", G->vals[j].center_val);

  fprintf(output,"\n");

}


void write_to_file(double * Q, struct grid *G, FILE *output)
{
  for (int j=1; j < M; j++)
    fprintf(output, "%3g ", Q[j]/G->vals[j].center_val);

  fprintf(output,"\n");
}

