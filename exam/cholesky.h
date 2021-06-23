#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

void ch_decomp(gsl_matrix* A, gsl_matrix* L);
void matrix_print(char s[],gsl_matrix* A);

