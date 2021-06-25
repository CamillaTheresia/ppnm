#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<stdlib.h>
#define RND (double)rand()/RAND_MAX

double funs(int i, double x){
	switch(i){
		case 0: return 1; break;
		case 1: return x; break;
		case 2: return x*x; break;
		default: return NAN;
	}
}

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

void GS_solve(gsl_vector* x, gsl_vector* y, gsl_vector* dy, gsl_vector* c,int m){
	gsl_matrix* Q=gsl_matrix_alloc(x->size,m);
	gsl_matrix* R=gsl_matrix_calloc(m,m);
	gsl_matrix* Qt=gsl_matrix_alloc(m,x->size);
	gsl_vector* b=gsl_vector_alloc(x->size);

	for(int i=0;i<x->size;i++){
		gsl_vector_set(b,i,gsl_vector_get(y,i)/gsl_vector_get(dy,i));
		for(int n=0;n<m;n++){
			gsl_matrix_set(Q,i,n,funs(n,gsl_vector_get(x,i))/gsl_vector_get(dy,i));
		}
	}

	GS_decomp(Q,R);

	for(int j=0;j<Qt->size2;j++){
		for(int i=0;i<Qt->size1;i++){
			double Qtij=gsl_matrix_get(Q,j,i);
			gsl_matrix_set(Qt,i,j,Qtij);
		}
	}
	for(int i=0;i<Qt->size1;i++){
		double ci=0;
		for(int k=0;k<b->size;k++){
			ci+=gsl_matrix_get(Qt,i,k)*gsl_vector_get(b,k);
		}
		gsl_vector_set(c,i,ci);
	}
	for(int i=c->size-1;i>=0;i--){
		double s=gsl_vector_get(c,i);
		for(int k=i+1;k<c->size;k++){
			s-=gsl_matrix_get(R,i,k)*gsl_vector_get(c,k);
		}
		gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));
	}
gsl_matrix_free(Q);
gsl_matrix_free(R);
gsl_matrix_free(Qt);
gsl_vector_free(b);
}

void GS_solvunc(gsl_vector* x, gsl_vector* y, gsl_vector* dy, gsl_vector* c,int m ,gsl_matrix* AtAinv){
	gsl_matrix* Q=gsl_matrix_alloc(x->size,m);
	gsl_matrix* R=gsl_matrix_calloc(m,m);
	gsl_matrix* Qt=gsl_matrix_alloc(m,x->size);
	gsl_vector* b=gsl_vector_alloc(x->size);
	gsl_matrix* A=gsl_matrix_alloc(x->size,m);
	gsl_matrix* At=gsl_matrix_alloc(m,x->size);
	gsl_matrix* AtA=gsl_matrix_alloc(m,m);
	gsl_matrix* Ra=gsl_matrix_calloc(m,m);
	gsl_matrix* AtAT=gsl_matrix_alloc(m,m);
	gsl_vector* dc1=gsl_vector_alloc(m);
	gsl_vector* dc2=gsl_vector_alloc(m);
	gsl_vector* e1=gsl_vector_calloc(m);
	gsl_vector* e2=gsl_vector_calloc(m);

	gsl_vector_set(e1,0,1);
	gsl_vector_set(e2,1,1);

	for(int i=0;i<x->size;i++){
		gsl_vector_set(b,i,gsl_vector_get(y,i)/gsl_vector_get(dy,i));
		for(int n=0;n<m;n++){
			gsl_matrix_set(Q,i,n,funs(n,gsl_vector_get(x,i))/gsl_vector_get(dy,i));
		}
	}
	gsl_matrix_memcpy(A,Q);
	GS_decomp(Q,R);

	for(int j=0;j<Qt->size2;j++){
		for(int i=0;i<Qt->size1;i++){
			double Qtij=gsl_matrix_get(Q,j,i);
			gsl_matrix_set(Qt,i,j,Qtij);
		}
	}
	for(int j=0;j<At->size2;j++){
		for(int i=0;i<At->size1;i++){
			double Atij=gsl_matrix_get(A,j,i);
			gsl_matrix_set(At,i,j,Atij);
		}
	}
	for(int i=0;i<Qt->size1;i++){
		double ci=0;
		for(int k=0;k<b->size;k++){
			ci+=gsl_matrix_get(Qt,i,k)*gsl_vector_get(b,k);
		}
		gsl_vector_set(c,i,ci);
	}
	for(int i=c->size-1;i>=0;i--){
		double s=gsl_vector_get(c,i);
		for(int k=i+1;k<c->size;k++){
			s-=gsl_matrix_get(R,i,k)*gsl_vector_get(c,k);
		}
		gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));
	}

	matrix_product(At,A,AtA);
	GS_decomp(AtA,Ra);
	for(int j=0;j<AtAT->size2;j++){
		for(int i=0;i<AtAT->size1;i++){
			double Atij=gsl_matrix_get(AtA,j,i);
			gsl_matrix_set(AtAT,i,j,Atij);
		}
	}
	for(int i=0;i<AtAT->size1;i++){
		double ci=0;
		for(int k=0;k<e1->size;k++){
			ci+=gsl_matrix_get(AtAT,i,k)*gsl_vector_get(e1,k);
		}
		gsl_vector_set(dc1,i,ci);
	}
	for(int i=dc1->size-1;i>=0;i--){
		double s=gsl_vector_get(dc1,i);
		for(int k=i+1;k<dc1->size;k++){
			s-=gsl_matrix_get(Ra,i,k)*gsl_vector_get(dc2,k);
		}
		gsl_vector_set(dc1,i,s/gsl_matrix_get(Ra,i,i));
	}
	gsl_matrix_set(AtAinv,0,0,gsl_vector_get(dc1,0));
	gsl_matrix_set(AtAinv,1,0,gsl_vector_get(dc1,1));
	for(int i=0;i<AtAT->size1;i++){
		double ci=0;
		for(int k=0;k<e1->size;k++){
			ci+=gsl_matrix_get(AtAT,i,k)*gsl_vector_get(e2,k);
		}
		gsl_vector_set(dc2,i,ci);
	}
	for(int i=dc1->size-1;i>=0;i--){
		double s=gsl_vector_get(dc2,i);
		for(int k=i+1;k<dc1->size;k++){
			s-=gsl_matrix_get(Ra,i,k)*gsl_vector_get(dc2,k);
		}
		gsl_vector_set(dc2,i,s/gsl_matrix_get(Ra,i,i));
	}
	gsl_matrix_set(AtAinv,0,1,gsl_vector_get(dc2,0));
	gsl_matrix_set(AtAinv,1,1,gsl_vector_get(dc2,1));

gsl_matrix_free(Q);
gsl_matrix_free(R);
gsl_matrix_free(Qt);
gsl_vector_free(b);
gsl_matrix_free(A);
gsl_matrix_free(At);
gsl_matrix_free(AtA);
gsl_vector_free(dc1);
gsl_vector_free(dc2);
gsl_vector_free(e1);
gsl_vector_free(e2);
gsl_matrix_free(Ra);
gsl_matrix_free(AtAT);
}

int main(){
	int n=10;
	int m=3;
	int nt=9;
	int mt=2;
	gsl_matrix* A=gsl_matrix_alloc(n,m);
	gsl_matrix* Acopy=gsl_matrix_alloc(n,m);
	gsl_matrix* R=gsl_matrix_calloc(m,m);
	gsl_matrix* At=gsl_matrix_alloc(m,n);
	gsl_matrix* ar=gsl_matrix_alloc(n,m);
	gsl_matrix* AtA=gsl_matrix_alloc(m,m);
	gsl_vector* x=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_alloc(n);
	gsl_vector* dy=gsl_vector_alloc(n);
	gsl_vector* c=gsl_vector_alloc(m);
	gsl_vector* t=gsl_vector_alloc(nt);
	gsl_vector* yt=gsl_vector_alloc(nt);
	gsl_vector* dyt=gsl_vector_alloc(nt);
	gsl_vector* ct=gsl_vector_alloc(mt);
	gsl_matrix* AtAinv=gsl_matrix_alloc(mt,mt);
	gsl_vector* ctunc=gsl_vector_alloc(mt);

	for(int i=0; i< A->size1; i++)
		for(int j=0; j<A->size2; j++)
		{
		double Aij=RND;
		gsl_matrix_set(A,i,j,Aij);
		}
	gsl_matrix_memcpy(Acopy,A);
	GS_decomp(A,R);
	printf("# index 0: Checks that QR-decomposition works\n");
	matrix_transpose(A,At);
	matrix_product(At,A,AtA);
	matrix_print("QtQ is:",AtA);
	matrix_product(A,R,ar);
	matrix_print("QR is:",ar);
	matrix_print("A is:",Acopy);
	printf("\n\n");

	for(int i=0;i<n;i++){
		gsl_vector_set(x,i,i+1);
		gsl_vector_set(y,i,2*(i+1)+1);
		gsl_vector_set(dy,i,1./2);
	}
	GS_solve(x,y,dy,c,m);
	printf("# index 1: ck's\n");
	vector_print("ck's:",c);
	printf("\n\n");
	printf("# index 2: least square\n");
	for(int i=0;i<n;i++){
		double f=gsl_vector_get(c,0)*funs(0,gsl_vector_get(x,i))+gsl_vector_get(c,1)*funs(1,gsl_vector_get(x,i))+gsl_vector_get(c,2)*funs(2,gsl_vector_get(x,i));
		printf("%g %g %g %g\n",gsl_vector_get(x,i),gsl_vector_get(y,i),gsl_vector_get(dy,i),f);
	}
	printf("\n\n");

	gsl_vector_set(t,0,1);
	gsl_vector_set(t,1,2);
	gsl_vector_set(t,2,3);
	gsl_vector_set(t,3,4);
	gsl_vector_set(t,4,6);
	gsl_vector_set(t,5,9);
	gsl_vector_set(t,6,10);
	gsl_vector_set(t,7,13);
	gsl_vector_set(t,8,15);

	gsl_vector_set(yt,0,log(117));
	gsl_vector_set(yt,1,log(100));
	gsl_vector_set(yt,2,log(88));
	gsl_vector_set(yt,3,log(72));
	gsl_vector_set(yt,4,log(53));
	gsl_vector_set(yt,5,log(29.5));
	gsl_vector_set(yt,6,log(25.2));
	gsl_vector_set(yt,7,log(15.2));
	gsl_vector_set(yt,8,log(11.1));

	for(int i=0;i<nt;i++){
		gsl_vector_set(dyt,i,1./20);
	}
	GS_solve(t,yt,dyt,ct,mt);
	printf("# index 3: ck's for exponential function\n");
	printf("Coefficients: a=%g, gamma=%g, gamma_real=0.278\n",exp(gsl_vector_get(ct,0)),-1*gsl_vector_get(ct,1));
	printf("The calculated value does not correspond with the modern value.\n");
	printf("\n\n");
	printf("# index 4: least square\n");
	for(int i=0;i<nt;i++){
		double ft=gsl_vector_get(ct,0)*funs(0,gsl_vector_get(t,i))+gsl_vector_get(ct,1)*funs(1,gsl_vector_get(t,i));
		printf("%g %g %g %g\n",gsl_vector_get(t,i),exp(gsl_vector_get(yt,i)),exp(gsl_vector_get(yt,i))/20,exp(ft));
	}
	printf("\n\n");
	GS_solvunc(t,yt,dyt,ctunc,mt,AtAinv);
	printf("# index 5: covariance matrix\n");
	matrix_print("Covariance matrix:",AtAinv);
	printf("\n\n");
	printf("# index 6: ck's for exponential function with uncertainties\n");
	printf("Coefficients: a=%g, a_unc=%g, gamma=%g, gamma_unc=%g, gamma_real=0.278\n",exp(gsl_vector_get(ctunc,0)),exp(sqrt(gsl_matrix_get(AtAinv,0,0))),-1*gsl_vector_get(ctunc,1),sqrt(gsl_matrix_get(AtAinv,1,1)));
	printf("Even with uncertainties, the calculated value does not correspond to the modern value.\n");
	printf("\n\n");
	printf("# index 7: least square\n");
	for(int i=0;i<nt;i++){
		double ft=gsl_vector_get(ctunc,0)*funs(0,gsl_vector_get(t,i))+gsl_vector_get(ctunc,1)*funs(1,gsl_vector_get(t,i));
		double ftp=(gsl_vector_get(ctunc,0)+sqrt(gsl_matrix_get(AtAinv,0,0)))*funs(0,gsl_vector_get(t,i))+(gsl_vector_get(ctunc,1)+sqrt(gsl_matrix_get(AtAinv,1,1)))*funs(1,gsl_vector_get(t,i));
		double ftm=(gsl_vector_get(ctunc,0)-sqrt(gsl_matrix_get(AtAinv,0,0)))*funs(0,gsl_vector_get(t,i))+(gsl_vector_get(ctunc,1)-sqrt(gsl_matrix_get(AtAinv,1,1)))*funs(1,gsl_vector_get(t,i));
		printf("%g %g %g %g %g %g\n",gsl_vector_get(t,i),exp(gsl_vector_get(yt,i)),exp(gsl_vector_get(yt,i))/20,exp(ft),exp(ftp),exp(ftm));
	}
	printf("\n\n");
gsl_matrix_free(A);
gsl_matrix_free(Acopy);
gsl_matrix_free(R);
gsl_matrix_free(At);
gsl_matrix_free(ar);
gsl_matrix_free(AtA);
gsl_vector_free(x);
gsl_vector_free(y);
gsl_vector_free(dy);
gsl_vector_free(c);
gsl_vector_free(t);
gsl_vector_free(yt);
gsl_vector_free(dyt);
gsl_vector_free(ct);
gsl_matrix_free(AtAinv);
gsl_vector_free(ctunc);
return 0;
}
