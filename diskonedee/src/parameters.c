#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "disk.h"
#include "global.h"

#define MAXLEN_PARAM_TAG 30
#define MAXLEN_PARAM_VALUE 30
#define MAX_PARAMETERS 20
#define REAL 1
#define STRING 2
#define INT 3

void read_params(void)
{
  char *in_file = "param.txt";


  char paramtag[MAX_PARAMETERS][MAXLEN_PARAM_TAG];
  int id[MAX_PARAMETERS];
  void *addr[MAX_PARAMETERS];
  char buf[MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 200];
  char buf1[MAXLEN_PARAM_TAG + 200], buf2[MAXLEN_PARAM_VALUE + 200], buf3[MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 400];
    
  int param_handled[MAX_PARAMETERS];
  for(int i = 0; i < MAX_PARAMETERS; i++)
      param_handled[i] = 0;

  /* Allocate string tags to parameters */
  int nt=0;
  strcpy(paramtag[nt], "TimeMax");
  addr[nt] = &params.TimeMax;
  id[nt++] = REAL;

  strcpy(paramtag[nt], "dt");
  addr[nt] = &params.dt;
  id[nt++] = REAL;

  strcpy(paramtag[nt], "Ngrid");
  addr[nt] = &params.Ngrid;
  id[nt++] = INT;

  strcpy(paramtag[nt], "AlphaCoefficient");
  addr[nt] = &params.AlphaCoefficient;
  id[nt++] = REAL;
  
  strcpy(paramtag[nt], "VerticalAspectRatio");
  addr[nt] = &params.VerticalAspectRatio;
  id[nt++] = REAL;

  strcpy(paramtag[nt], "TempProfileIndex");
  addr[nt] = &params.TempProfileIndex;
  id[nt++] = REAL;

  strcpy(paramtag[nt], "InnerBoundaryConditionType");
  addr[nt] = &params.InnerBoundaryConditionType;
  id[nt++] = INT;

  strcpy(paramtag[nt], "OuterBoundaryConditionType");
  addr[nt] = &params.OuterBoundaryConditionType;
  id[nt++] = INT;  

  strcpy(paramtag[nt], "ExternalSources");
  addr[nt] = &params.ExternalSources;
  id[nt++] = INT;

  strcpy(paramtag[nt], "ExternalAccretionStrength");
  addr[nt] = &params.ExternalAccretionStrength;
  id[nt++] = REAL;
  
  strcpy(paramtag[nt], "ExternalTorqueStrength");
  addr[nt] = &params.ExternalTorqueStrength;
  id[nt++] = REAL;

  strcpy(paramtag[nt], "Softening");
  addr[nt] = &params.Softening;
  id[nt++] = REAL;
  
  /* Read from file */
  FILE *fd;
  int j;

  if((fd = fopen(in_file, "r")))
    {
      printf("Obtaining parameters from file '%s'.\n", in_file);
      while(!feof(fd))
	{
	  *buf = 0;
	  fgets(buf, MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 200, fd);

	  if(sscanf(buf, "%s%s", buf1, buf2, buf3) < 2)
	    continue;
	  
	  if(buf1[0] == '%')
	    continue;

	  for(int k = 0; k < nt; k++)
	    if(strcmp(buf1, paramtag[k]) == 0)
	      {
		if(param_handled[k] == 0)
		  {
		    j = k;
		    param_handled[k] = 1;
		    break;
		  }
		else
		  {
		    j = -2;
		    break;
		  }
	      }

	      
	  if(j >= 0)
	    {
	      switch (id[j])
		{
		case REAL:
		  *((double *) addr[j]) = atof(buf2);
		  break;
		case STRING:
		  strcpy((char *) addr[j], buf2);
		  break;
		case INT:
		  *((int *) addr[j]) = atoi(buf2);
		  break;
		}
	    }

	  
	}
      fclose(fd);
      
      printf("Your parameters are:\n");
      for(int k = 0; k < nt; k++)
	if(param_handled[k] == 1)
	  {
	    if (id[k] == REAL)
	      printf("%s=%g\n",paramtag[k],*((double *) addr[k]));
	    if (id[k] == INT)
	      printf("%s=%d\n",paramtag[k],*((int *) addr[k]));
	    if (id[k] == STRING)
	      printf("%s=%s\n",paramtag[k],(char *) addr[k]);
	       

	  }
      
    }

}
