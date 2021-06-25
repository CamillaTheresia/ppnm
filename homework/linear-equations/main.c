#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#define RND (double)(rand()%10)
#define FMT "%7.3f"

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

void matrix_print(char s[],gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0; i<A->size1; i++){
    		for(int j=0;j<A->size2;j++){
      			printf("index(%d,%d) = %g\n",i,j,gsl_matrix_get(A,i,j));
		}
	printf("\n");
	}
}

void vector_print(char s[],gsl_vector* v){
	printf("%s\n",s);
	for(int i=0; i<v->size; i++){
		printf("index(%d) = %g\n",i,gsl_vector_get(v,i));
	}
	printf("\n");
}

void matrix_product(gsl_matrix* A, gsl_matrix* R, gsl_matrix* P){
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<R->size2;j++){
			double Pij=0;
			for(int k=0;k<A->size2;k++){
				Pij+=gsl_matrix_get(A,i,k)*gsl_matrix_get(R,k,j);
			gsl_matrix_set(P,i,j,Pij);
			}
		}
	}
}

void mv_product(gsl_matrix* A, gsl_vector* x, gsl_vector* P){
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<x->size;j++){
			double xi=0;
			for(int k=0;k<A->size2;k++){
				xi+=gsl_matrix_get(A,i,k)*gsl_vector_get(x,k);
			gsl_vector_set(P,i,xi);
			}
		}
	}
}

void matrix_transpose(gsl_matrix* A, gsl_matrix* At){
	for(int i=0;i<A->size2;i++){
		for(int j=0;j<A->size1;j++){
			double Atij=gsl_matrix_get(A,j,i);
			gsl_matrix_set(At,i,j,Atij);
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

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
	gsl_matrix* Qt=gsl_matrix_alloc(Q->size2,Q->size1);
	gsl_vector* e1=gsl_vector_calloc(B->size2);
	gsl_vector* e2=gsl_vector_calloc(B->size2);
	gsl_vector* e3=gsl_vector_calloc(B->size2);
	gsl_vector* x1=gsl_vector_alloc(B->size2);
	gsl_vector* x2=gsl_vector_alloc(B->size2);
	gsl_vector* x3=gsl_vector_alloc(B->size2);

	gsl_vector_set(e1,0,1);
	gsl_vector_set(e2,1,1);
	gsl_vector_set(e3,2,1);

	GS_solve(Q, Qt, R, e1, x1);
	GS_solve(Q, Qt, R, e2, x2);
	GS_solve(Q, Qt, R, e3, x3);

	for(int i=0;i<B->size1;i++){
	gsl_matrix_set(B,i,0,gsl_vector_get(x1,i));
	}
	for(int i=0;i<B->size1;i++){
	gsl_matrix_set(B,i,1,gsl_vector_get(x2,i));
	}
	for(int i=0;i<B->size1;i++){
	gsl_matrix_set(B,i,2,gsl_vector_get(x3,i));
	}

gsl_matrix_free(Qt);
gsl_vector_free(e1);
gsl_vector_free(e2);
gsl_vector_free(e3);
gsl_vector_free(x1);
gsl_vector_free(x2);
gsl_vector_free(x3);
}

int main(int argc, char** argv){
	int n=4;
	int m=3;
if(argc>1)n=atoi(argv[1]);
if(argc>2)m=atoi(argv[2]);
if(n<7){
	gsl_matrix* A=gsl_matrix_alloc(n,m);
	gsl_matrix* Acopy=gsl_matrix_alloc(n,m);
	gsl_matrix* R=gsl_matrix_calloc(m,m);
	gsl_matrix* qr=gsl_matrix_alloc(n,m);
	gsl_matrix* At=gsl_matrix_alloc(m,n);
	gsl_matrix* i=gsl_matrix_alloc(m,m);
	gsl_matrix* a1=gsl_matrix_alloc(m,m);
	gsl_matrix* acopy1=gsl_matrix_alloc(m,m);
	gsl_matrix* r1=gsl_matrix_calloc(m,m);
	gsl_matrix* at1=gsl_matrix_alloc(m,m);
	gsl_vector* b=gsl_vector_alloc(m);
	gsl_vector* x=gsl_vector_alloc(m);
	gsl_vector* a1x=gsl_vector_alloc(m);
	gsl_matrix* a2=gsl_matrix_alloc(m,m);
	gsl_matrix* r2=gsl_matrix_calloc(m,m);
	gsl_matrix* acopy2=gsl_matrix_alloc(m,m);
	gsl_matrix* B=gsl_matrix_alloc(m,m);
	gsl_matrix* ab=gsl_matrix_alloc(m,m);
	gsl_matrix* ba=gsl_matrix_alloc(m,m);
	for(int i=0; i< A->size1; i++){
		for(int j=0; j<A->size2; j++){
			double Aij=RND;
			gsl_matrix_set(A,i,j,Aij);
		}
	}
	gsl_matrix_memcpy(Acopy,A);
	GS_decomp(A,R);
	matrix_print("Checks that R is upper triangular:",R);
	matrix_product(A,R,qr);
	matrix_transpose(A,At);
	matrix_product(At,A,i);
	matrix_print("Checks that Q^T*Q is equal to 1:",i);
	matrix_print("Q*R is equal to:",qr);
	matrix_print("A is equal to:",Acopy);
	matrix_print("Q:",A);
	for(int i=0; i< a1->size1; i++){
		for(int j=0; j<a1->size2; j++){
			double aij=RND;
			gsl_matrix_set(a1,i,j,aij);
		}
	}
	gsl_matrix_memcpy(acopy1,a1);
	GS_decomp(a1,r1);
	for(int i=0;i<b->size;i++){
		double bi=RND;
		gsl_vector_set(b,i,bi);
	}
	GS_solve(a1,at1,r1,b,x);
	printf("Solving the equation QRx=b:");
	vector_print("x is found to be:",x);
	mv_product(acopy1,x,a1x);
	vector_print("Ax is:",a1x);
	vector_print("b is:",b);
	for(int i=0; i< a2->size1; i++){
		for(int j=0; j<a2->size2; j++){
			double aij=RND;
			gsl_matrix_set(a2,i,j,aij);
		}
	}
	gsl_matrix_memcpy(acopy2,a2);
	GS_decomp(a2,r2);
	GS_inverse(a2,r2,B);
	matrix_print("The inverse matrix B is:",B);
	matrix_product(acopy2,B,ab);
	matrix_product(B,acopy2,ba);
	matrix_print("AB is:",ab);
	matrix_print("BA is:",ba);
gsl_matrix_free(A);
gsl_matrix_free(Acopy);
gsl_matrix_free(R);
gsl_matrix_free(qr);
gsl_matrix_free(At);
gsl_matrix_free(i);
gsl_matrix_free(a1);
gsl_matrix_free(acopy1);
gsl_matrix_free(r1);
gsl_matrix_free(at1);
gsl_vector_free(b);
gsl_vector_free(x);
gsl_vector_free(a1x);
gsl_matrix_free(a2);
gsl_matrix_free(r2);
gsl_matrix_free(acopy2);
gsl_matrix_free(B);
gsl_matrix_free(ab);
gsl_matrix_free(ba);
}
else{
	printf("n=%i\n",n);
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* R=gsl_matrix_alloc(n,n);
	for(int i=0;i<m;i++)for(int j=0;j<m;j++)gsl_matrix_set(A,i,j,RND);
	GS_decomp(A,R);
	gsl_matrix_free(A);
	gsl_matrix_free(R);
}
}
