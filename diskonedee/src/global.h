#include <stddef.h>
#include <stdio.h>



extern int M;
extern int N;
extern double timemax;
extern double Rmin;
extern double Rmax;


extern struct parameters
{
  double AlphaCoefficient;
  double VerticalAspectRatio;
  double TempProfileIndex;
}
  params ;

//struct grid;


struct gridvals
{
  double center_val;
  double side_val;
  double A_coeff_val;
  double B_coeff_val;
  double C_coeff_val;
};

struct grid
{
  int Npoints;
  struct gridvals vals[];
  //double Alpha[];
  //double Beta[];
};
