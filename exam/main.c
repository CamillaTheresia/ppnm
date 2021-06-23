#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "cholesky.h"
#define RND (double)rand()/RAND_MAX

int main(){
	int n = 4;
	gsl_matrix* B = gsl_matrix_calloc(n,n);
	gsl_matrix* A = gsl_matrix_calloc(n,n);
	gsl_matrix* L = gsl_matrix_calloc(n,n);
	gsl_matrix* LLt = gsl_matrix_calloc(n,n);
	for(int i=0; i<B->size1; i++){
		for(int j=0; j<B->size2; j++){
			double num = RND;
			gsl_matrix_set(B,i,j,num);
		}
	}
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,B,B,0,A);
	matrix_print("A is a symmetric positive-definite matrix:",A);
	ch_decomp(A,L);
	matrix_print("L is a lower triangular matrix:",L);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,L,L,0,LLt);
	matrix_print("LL^T calculated to check that A=LL^T:",LLt);

gsl_matrix_free(A);
gsl_matrix_free(L);
return 0;
}
