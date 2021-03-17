#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>

int main(){
	int m=100000;
	gsl_matrix* A=gsl_matrix_alloc(m,m);
	gsl_vector* v=gsl_vector_alloc(m);
	gsl_linalg_QR_decomp(A,v);
return 0;
}
