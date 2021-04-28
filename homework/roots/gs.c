#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

void GS_decomp(gsl_matrix* A, gsl_matrix* R){
        for(int i=0;i<A->size2;i++){
                double s=0;
                for(int j=0;j<A->size1;j++){
                        s+=gsl_matrix_get(A,j,i)*gsl_matrix_get(A,j,i);
                }
                double Rii=sqrt(s);
                gsl_matrix_set(R,i,i,Rii);
                for(int j=0;j<A->size1;j++){
                        double qji=gsl_matrix_get(A,j,i)/Rii;
                        gsl_matrix_set(A,j,i,qji);
                }
                for(int j=i+1;j<A->size2;j++){
                        double d=0;
                        for(int k=0;k<A->size1;k++){
                                d+=gsl_matrix_get(A,k,i)*gsl_matrix_get(A,k,j);
                        }
                        double Rij=d;
                        gsl_matrix_set(R,i,j,Rij);
                        for(int k=0;k<A->size1;k++){
                                double akj=gsl_matrix_get(A,k,j)-gsl_matrix_get(A,k,i)*Rij;
                                gsl_matrix_set(A,k,j,akj);
                        }
                }
        }
}

void GS_solve(gsl_matrix* Q, gsl_matrix* Qt, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
        for(int i=0;i<Qt->size2;i++){
                for(int j=0;j<Qt->size1;j++){
                        double Qtij=gsl_matrix_get(Q,j,i);
                        gsl_matrix_set(Qt,i,j,Qtij);
                }
        }
        for(int i=0;i<Qt->size1;i++){
                double xi=0;
                for(int k=0;k<b->size;k++){
                        xi+=gsl_matrix_get(Qt,i,k)*gsl_vector_get(b,k);
                }
                gsl_vector_set(x,i,xi);
        }
        for(int i=x->size-1;i>=0;i--){
                double s=gsl_vector_get(x,i);
                for(int k=i+1;k<x->size;k++){
                        s-=gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
                }
                gsl_vector_set(x,i,s/gsl_matrix_get(R,i,i));
        }
}

