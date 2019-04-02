#include <stdlib.h>


//void init_quant(double * Q, double * G);

//void init_grid(double * Q);
//void init_grid(struct grid *G);
void init_quant(double *Q);//, struct grid *G);
struct grid *initGrid(struct grid *G, double R1, double R2, int Npoints);
void advance(double * Q, double * G, double dt);
