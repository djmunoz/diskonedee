#include <stdio.h>
#include <time.h>
#include <math.h>

#include "disk.h"
#include "global.h"




void write_header(struct grid *G, FILE *output)
{
  fprintf(output,"000\t");
  for (int j=1; j < M; j++)
    fprintf(output, "%3g ", G->vals[j].center_val);
  fprintf(output,"\n");  

  fprintf(output,"000\t");
  for (int j=1; j < M; j++)
    fprintf(output, "%3g ", eval_omega_func(G->vals[j].center_val));
  fprintf(output,"\n");
  
  fprintf(output,"000\t");
  for (int j=1; j < M; j++)
    fprintf(output, "%3g ", eval_omegaprime_func(G->vals[j].center_val));
  fprintf(output,"\n");

  fprintf(output,"000\t");
  for (int j=1; j < M; j++)
    fprintf(output, "%3g ", eval_lprime_func(G->vals[j].center_val));
  fprintf(output,"\n");
    
  fprintf(output,"000\t");
  for (int j=1; j < M; j++)
    fprintf(output, "%3g ", eval_nu_func(G->vals[j].center_val));
  fprintf(output,"\n");

  fprintf(output,"000\t");
  for (int j=1; j < M; j++)
    fprintf(output, "%3g ", eval_h_func(G->vals[j].center_val));
  fprintf(output,"\n");

}


void write_to_file(double * Q, struct grid *G, FILE *output, double time)
{
  fprintf(output, "%3g\t", time);
  for (int j=1; j < M; j++)
    fprintf(output, "%3g ", Q[j]/G->vals[j].center_val);

  fprintf(output,"\n");
}

