#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#define RND (double)rand()/RAND_MAX

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

int main(){
	int m=100000;
	gsl_matrix* A=gsl_matrix_alloc(m,m);
	gsl_matrix* R=gsl_matrix_calloc(m,m);
	for(int i=0; i< A->size1; i++){
		for(int j=0; j<A->size2; j++){
			double Aij=RND;
			gsl_matrix_set(A,i,j,Aij);
		}
	}
	GS_decomp(A,R);
gsl_matrix_free(A);
gsl_matrix_free(R);
return 0;
}
