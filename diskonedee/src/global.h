#include <stddef.h>
#include <stdio.h>

extern int M;
extern int N;
extern double dt;
extern char ParameterFile[40];

extern double Rmin;
extern double Rmax;

extern struct parameters
{
  char OutputFile[50];
  double TimeMax;
  double dt;
  double Rmin;
  double Rmax;
  int Ngrid;
  double AlphaCoefficient;
  double VerticalAspectRatio;
  double TempProfileIndex;
  int InnerBoundaryConditionType;
  int OuterBoundaryConditionType;
  int InitialProfile;
  int ExternalSources;
  double ExternalSourcesRadius;
  double ExternalSourcesWidth;
  double ExternalAccretionStrength;
  double ExternalTorqueStrength;
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
