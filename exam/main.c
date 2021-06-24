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
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* b = gsl_vector_alloc(n);
	gsl_vector* Ax = gsl_vector_alloc(n);
	gsl_matrix* Ainv = gsl_matrix_calloc(n,n);
	gsl_matrix* AAinv = gsl_matrix_alloc(n,n);
	gsl_matrix* AinvA = gsl_matrix_alloc(n,n);
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

	for(int i=0; i<b->size;i++){
		double num = RND;
		gsl_vector_set(b,i,num);
	}
	ch_solve(L,b,x);
	printf("Solving the linear equation Ax=b:\n\n");
	vector_print("b is:",b);
	vector_print("x is found to be:",x);
	gsl_blas_dgemv(CblasNoTrans,1,A,x,0,Ax);
	vector_print("Ax is calculated to check that Ax=b:",Ax);

	double det = ch_det(L);
	printf("The determinant of A is found by det(A)=det(L)*det(L^T)=det(L)^2=\n");
	printf("(sum of the diagonal elements of L)^2:\n\n");
	printf("det(A)=%10g\n\n\n",det);

	ch_inv(L,Ainv);
	matrix_print("The inverse of A is:",Ainv);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,A,Ainv,0,AAinv);
	matrix_print("A*Ainv should be equal to the identity matrix:",AAinv);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Ainv,A,0,AinvA);
	matrix_print("Ainv*A should be equal to the identity matrix:",AinvA);

gsl_matrix_free(A);
gsl_matrix_free(L);
gsl_matrix_free(B);
gsl_matrix_free(LLt);
gsl_vector_free(x);
gsl_vector_free(b);
gsl_vector_free(Ax);
gsl_matrix_free(Ainv);
gsl_matrix_free(AAinv);
gsl_matrix_free(AinvA);
return 0;
}

