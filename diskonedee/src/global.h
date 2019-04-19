#include <stddef.h>
#include <stdio.h>

extern int M;
extern int N;
extern double Rmin;
extern double Rmax;
extern double dt;

extern struct parameters
{
  double TimeMax;
  double AlphaCoefficient;
  double VerticalAspectRatio;
  double TempProfileIndex;
  int InnerBoundaryConditionType;
  int OuterBoundaryConditionType;
  int ExternalSources;
  double Softening;
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
  double cn_upper_diag;
  double cn_middle_diag;
  double cn_lower_diag;
  double cn_rhs;
};

struct grid
{
  int Npoints;
  double val_in;
  double val_out;
  struct gridvals vals[];
};
