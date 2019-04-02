#include <stddef.h>
#include <stdio.h>



extern int M;
extern int N;
extern double timemax;
extern double Rmin;
extern double Rmax;


extern struct parameters
{
  double AlphaCofficient;
  double VerticalAspectRatio;
  double TempProfileIndex;
}
  params ;

//struct grid;


struct gridvals
{
  double radius_val;
  double alpha_val;
  double beta_val;
};

struct grid
{
  int Npoints;
  struct gridvals vals[];
  //double Alpha[];
  //double Beta[];
};
