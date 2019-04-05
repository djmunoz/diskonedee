#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

#include "disk.h"
#include "global.h"


gsl_matrix *invert_a_matrix(gsl_matrix *matrix)
{
  
    gsl_permutation *p = gsl_permutation_alloc(M);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(M, M);
    gsl_linalg_LU_invert(matrix, p, inv);
    gsl_permutation_free(p);

    return inv;
}
