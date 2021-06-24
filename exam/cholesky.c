#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "cholesky.h"

void ch_decomp(gsl_matrix* A, gsl_matrix* L){
	for(int i=0; i<A->size1;i++){
		for(int j=0; j<=i; j++){
			double sum=0;
			for(int k=0; k<j; k++){
				sum += gsl_matrix_get(L,i,k)*gsl_matrix_get(L,j,k);
			}
			if (i==j){
				double a = gsl_matrix_get(A,i,i)-sum;
				double b = sqrt(a);
				gsl_matrix_set(L,i,j,b);
			}
			else {
				double c = (1.0/gsl_matrix_get(L,j,j)*(gsl_matrix_get(A,i,j)-sum));
				gsl_matrix_set(L,i,j,c);
			}
		}
	}
}

void matrix_print(char s[],gsl_matrix* A){
	printf("%s\n\n",s);
	for(int i=0; i<A->size1; i++){
    		for(int j=0;j<A->size2;j++){
      			printf("%15g",gsl_matrix_get(A,i,j));
		}
		printf("\n\n");
	}
	printf("\n");
}

void ch_solve(gsl_matrix* L, gsl_vector* b, gsl_vector* x){
	gsl_vector* y = gsl_vector_alloc(b->size);
	gsl_vector_memcpy(y,b);
	for(int i=0; i<y->size; i++){
		double yi=gsl_vector_get(y,i);
		for(int k=0; k<i;k++){
			yi-=gsl_matrix_get(L,i,k)*gsl_vector_get(y,k);
		}
		gsl_vector_set(y,i,yi/gsl_matrix_get(L,i,i));
	}

	gsl_matrix* Lt = gsl_matrix_alloc(L->size2,L->size1);
	gsl_matrix_transpose_memcpy(Lt,L);
	gsl_vector_memcpy(x,y);
	for(int i=x->size-1;i>=0;i--){
		double xi=gsl_vector_get(x,i);
		for(int k=i+1; k<x->size; k++){
			xi-=gsl_matrix_get(Lt,i,k)*gsl_vector_get(x,k);
		}
		gsl_vector_set(x,i,xi/gsl_matrix_get(Lt,i,i));
	}

gsl_vector_free(y);
gsl_matrix_free(Lt);
}

void vector_print(char s[],gsl_vector* v){
	printf("%s\n\n",s);
	for(int i=0; i<v->size; i++){
		printf("%10g\n\n",gsl_vector_get(v,i));
	}
	printf("\n");
}

double ch_det(gsl_matrix* L){
	double d = 1;
	for(int i=0; i<L->size1;i++){
		d*=gsl_matrix_get(L,i,i);
	}
	return d*d;
}

void ch_inv(gsl_matrix* L, gsl_matrix* Ainv){
	gsl_vector* b = gsl_vector_calloc(L->size2);
	gsl_vector* x = gsl_vector_calloc(L->size2);
	for(int i=0; i<L->size2; i++){
		gsl_vector_set(b,i,1.0);
		ch_solve(L,b,x);
		gsl_vector_set(b,i,0.0);
		gsl_matrix_set_col(Ainv,i,x);
	}

gsl_vector_free(b);
gsl_vector_free(x);
}







