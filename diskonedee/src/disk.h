#include <stdlib.h>
#include <gsl/gsl_matrix_double.h>


//void init_quant(double * Q, double * G);

//void init_grid(double * Q);
//void init_grid(struct grid *G);

struct grid *init_grid(struct grid *G);
double *init_quant(double *Q, struct grid *G);
void init_grid_terms(double *Q, struct grid *G);
void init_boundaries(double *Q, struct grid *G);

double eval_omega_func(double R);
double eval_omegaprimeprime_func(double R);
double eval_omegaprime_func(double R);
double eval_h_func(double R);
double eval_nu_func(double R);
double eval_nuprime_func(double R);
double eval_l_func(double R);
double eval_lprime_func(double R);
double eval_g_func(double R);
double eval_beta_func(double R);
double eval_gamma_func(double R);

double eval_inner_boundary(double * Q, struct grid *G);
double eval_outer_boundary(double * Q, struct grid *G);


void advance(double * Q,  struct grid *G, double dt);
void advance_euler(double * Q, struct grid *G, double dt);
void advance_cranknicolson(double * Q, struct grid *G, double dt);

void read_params(void);
void write_header(struct grid *G, FILE *output);
void write_to_file(double * Q, struct grid *G, FILE *output, double time);


//double *invert_tridiagonal_problem(double *xx,  double * dd, double * du, double * dl, double * b);
void invert_tridiagonal_problem(double *Q, double * dd, double * du, double * dl, double * b);
void invert_tridiagonal_CN(double alpha, double * b, double * x, int M);
gsl_matrix *invert_a_matrix(gsl_matrix *matrix);


double get_profile(double R);
double get_profile_powerlaw_truncated(double R);
double get_profile_powerlaw(double R);
double get_profile_deltafunc(double R);
