#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

void ch_decomp(gsl_matrix* A, gsl_matrix* L);
void matrix_print(char s[],gsl_matrix* A);
void ch_solve(gsl_matrix* L, gsl_vector* b, gsl_vector* x);
void vector_print(char s[],gsl_vector* v);
double ch_det(gsl_matrix* L);
void ch_inv(gsl_matrix* L, gsl_matrix* Ainv);
