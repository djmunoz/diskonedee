#include <stdlib.h>
#include <gsl/gsl_matrix_double.h>


//void init_quant(double * Q, double * G);

//void init_grid(double * Q);
//void init_grid(struct grid *G);

struct grid *initGrid(struct grid *G, double R1, double R2, int Npoints);
double *init_quant(double *Q, struct grid *G);

double eval_omegaprime_func(double R);
double eval_nu_func(double R);
double eval_lprime_func(double R);
double eval_inner_boundary(double * Q, struct grid *G);
double eval_outer_boundary(double * Q, struct grid *G);


void advance(double * Q,  struct grid *G, double dt);
void advance_euler(double * Q, struct grid *G, double dt);
void advance_cranknicolson(double * Q, struct grid *G, double dt);

void write_header(struct grid *G, FILE *output);
void write_to_file(double * Q, struct grid *G, FILE *output);


//double *invert_tridiagonal_problem(double *xx,  double * dd, double * du, double * dl, double * b);
void invert_tridiagonal_problem(double *Q, double * dd, double * du, double * dl, double * b);
void invert_tridiagonal_CN(double alpha, double * b, double * x, int M);
gsl_matrix *invert_a_matrix(gsl_matrix *matrix);

