#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<gsl/gsl_blas.h>
#define RND (double)rand()/RAND_MAX

void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}

void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int j=0;j<A->size2;j++){
		double new_apj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		double new_aqj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_apj);
		gsl_matrix_set(A,q,j,new_aqj);
		}
}

void jacobi_diag(gsl_matrix* A, gsl_matrix* V, int n){
	int changed;
	do{
		changed=0;
		for(int p=0;p<n-1;p++)
		for(int q=p+1;q<n;q++){
			double apq=gsl_matrix_get(A,p,q);
			double app=gsl_matrix_get(A,p,p);
			double aqq=gsl_matrix_get(A,q,q);
			double theta=0.5*atan2(2*apq,aqq-app);
			double c=cos(theta),s=sin(theta);
			double new_app=c*c*app-2*s*c*apq+s*s*aqq;
			double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
			if(new_app!=app || new_aqq!=aqq) // do rotation
				{
				changed=1;
				timesJ(A,p,q, theta);
				Jtimes(A,p,q,-theta); // A←J^T*A*J 
				timesJ(V,p,q, theta); // V←V*J
				}
		}
	}while(changed!=0);
}

void matrix_transpose(gsl_matrix* A,gsl_matrix* B){
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++){
			double Bji=gsl_matrix_get(A,i,j);
			gsl_matrix_set(B,j,i,Bji);
		}
	}
}

void matrix_print(char s[], gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++){
			printf("matrix(%d,%d)=%g\n",i,j,gsl_matrix_get(A,i,j));
		}
	}
	printf("\n");
}

int main(){
	int n=5;
	int m=20;
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
	gsl_matrix* V=gsl_matrix_calloc(n,n);
	gsl_matrix* VtA=gsl_matrix_alloc(n,n);
	gsl_matrix* VtAV=gsl_matrix_alloc(n,n);
	gsl_matrix* Vd=gsl_matrix_alloc(n,n);
	gsl_matrix* VdVt=gsl_matrix_alloc(n,n);
	gsl_matrix* VtV=gsl_matrix_alloc(n,n);
	gsl_matrix* H=gsl_matrix_alloc(m,m);
	gsl_matrix* E=gsl_matrix_calloc(m,m);

	for(int i=0; i< A->size1; i++)
		for(int j=i; j<A->size2; j++)
		{
		double Aij=RND;
		gsl_matrix_set(A,i,j,Aij);
		gsl_matrix_set(A,j,i,Aij);
		}
	gsl_matrix_memcpy(Acopy,A);

	for(int i=0;i<V->size1;i++){
		gsl_matrix_set(V,i,i,1);
	}

	jacobi_diag(A,V,n);

	double s=1.0/(m+1);
	for(int i=0;i<m-1;i++){
		gsl_matrix_set(H,i,i,-2);
		gsl_matrix_set(H,i,i+1,1);
		gsl_matrix_set(H,i+1,i,1);
	}
	gsl_matrix_set(H,m-1,m-1,-2);
	gsl_matrix_scale(H,-1/s/s);

	for(int i=0;i<E->size1;i++){
		gsl_matrix_set(E,i,i,1);
	}

	jacobi_diag(H,E,m);

	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, Acopy, 0.0, VtA);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, VtA, V, 0.0, VtAV);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, A, 0.0, Vd);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Vd, V, 0.0, VdVt);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, V, 0.0, VtV);

	printf("# index 0: part A\n");
	matrix_print("VTAV:",VtAV);
	matrix_print("D:",A);
	matrix_print("VDVT:",VdVt);
	matrix_print("A:",Acopy);
	matrix_print("VTV:",VtV);
	printf("\n\n");

	printf("# index 1: check that the energies are correct\n");
	for (int k=0;k<m/3; k++){
		double exact = M_PI*M_PI*(k+1)*(k+1);
		double calculated = gsl_matrix_get(H,k,k);
		printf("%i %g %g\n",k,calculated,exact);
	}
	printf("\n\n");

	printf("# index 2: eigenfunctions\n");
	for(int k=0;k<3;k++){
		printf("0 0\n");
		for(int i=0;i<m;i++)
			printf("%g %g\n",(i+1.0)/(m+1), gsl_matrix_get(E,i,k));
		printf("1 0\n");
	}
	printf("\n\n");

	printf("# index 3: analytical\n");
	for(double i=0;i<=1;i+=1./16){
		printf("%g %g %g %g\n",i,0.32*sin(M_PI*i),-0.32*sin(2*M_PI*i),0.32*sin(3*M_PI*i));
	}
	printf("\n\n");

gsl_matrix_free(A);
gsl_matrix_free(Acopy);
gsl_matrix_free(V);
gsl_matrix_free(VtA);
gsl_matrix_free(VtAV);
gsl_matrix_free(Vd);
gsl_matrix_free(VdVt);
gsl_matrix_free(VtV);
gsl_matrix_free(H);
gsl_matrix_free(E);
return 0;
}





