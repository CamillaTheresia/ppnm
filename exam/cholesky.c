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
	printf("%s\n",s);
	for(int i=0; i<A->size1; i++){
    		for(int j=0;j<A->size2;j++){
      			printf("index(%d,%d) = %g\n",i,j,gsl_matrix_get(A,i,j));
		}
	printf("\n");
	}
}

