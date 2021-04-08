#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

void rkstep12(
	void f(double x, gsl_vector* yx, gsl_vector* dydx),
	double x, gsl_vector* yx, double h, gsl_vector* yh, gsl_vector* err,gsl_vector* k0, gsl_vector* k1
){
	f(x,yx,k0);
	gsl_vector_memcpy(yh,yx);
	gsl_blas_daxpy(h/2,k0,yh);
	f(x+h/2,yh,k1);
	gsl_vector_memcpy(yh,yx);
	gsl_blas_daxpy(h,k1,yh);
	gsl_vector_memcpy(err,k1);
	gsl_blas_daxpy(-1,k0,err);
}

void driver(
	void f(double x,gsl_vector* yx,gsl_vector* dydx),
	double a,
	gsl_vector* yx,
	double b,
	double h,
	double acc,
	double eps
){
	double x=a;
	fprintf(stderr,"%g ",x);
	for(int i=0;i<yx->size;i++) fprintf(stderr,"%g ",gsl_vector_get(yx,i));
	fprintf(stderr,"\n");
	gsl_vector* yh =gsl_vector_alloc(yx->size);
	gsl_vector* err=gsl_vector_alloc(yx->size);
	gsl_vector* k0 =gsl_vector_alloc(yx->size);
	gsl_vector* k1 =gsl_vector_alloc(yx->size);
	while(x<b){
		if(x+h>b)h=b-x;
		rkstep12(f,x,yx,h,yh,err,k0,k1);
		double norm_err=gsl_blas_dnrm2(err);
		double norm_y  =gsl_blas_dnrm2(yh);
		double tol=(norm_y*eps+acc)*sqrt(h/(b-a));
		if(norm_err<tol){
			x=x+h;
			gsl_vector_memcpy(yx,yh);
			fprintf(stderr,"%g ",x);
			for(int i=0;i<yx->size;i++)fprintf(stderr,"%g ",gsl_vector_get(yx,i));
			fprintf(stderr,"\n");
		}
		if(norm_err>0) h*=0.95*pow(tol/norm_err,0.25);
		else h*=2;
	}
gsl_vector_free(yh);
gsl_vector_free(err);
gsl_vector_free(k0);
gsl_vector_free(k1);
}

void f1(double x, gsl_vector* yx, gsl_vector* dydx){
	gsl_vector_set(dydx,0,gsl_vector_get(yx,1));
	gsl_vector_set(dydx,1,-gsl_vector_get(yx,0));
}

void f2(double x, gsl_vector* yx, gsl_vector* dydx){
	double N=6*1e6,Tc=2,Tr=14;
	gsl_vector_set(dydx,0,-gsl_vector_get(yx,1)*gsl_vector_get(yx,0)/(N*Tc));
	gsl_vector_set(dydx,1,gsl_vector_get(yx,1)*gsl_vector_get(yx,0)/(N*Tc)-gsl_vector_get(yx,1)/Tr);
	gsl_vector_set(dydx,2,gsl_vector_get(yx,1)/Tr);
}

void f4(double x, gsl_vector* yx, gsl_vector* dydx){
	double N=6*1e6,Tc=4,Tr=14;
	gsl_vector_set(dydx,0,-gsl_vector_get(yx,1)*gsl_vector_get(yx,0)/(N*Tc));
	gsl_vector_set(dydx,1,gsl_vector_get(yx,1)*gsl_vector_get(yx,0)/(N*Tc)-gsl_vector_get(yx,1)/Tr);
	gsl_vector_set(dydx,2,gsl_vector_get(yx,1)/Tr);
}
void f8(double x, gsl_vector* yx, gsl_vector* dydx){
	double N=6*1e6,Tc=8,Tr=14;
	gsl_vector_set(dydx,0,-gsl_vector_get(yx,1)*gsl_vector_get(yx,0)/(N*Tc));
	gsl_vector_set(dydx,1,gsl_vector_get(yx,1)*gsl_vector_get(yx,0)/(N*Tc)-gsl_vector_get(yx,1)/Tr);
	gsl_vector_set(dydx,2,gsl_vector_get(yx,1)/Tr);
}

int main(){
	double a=0,b=7,h=0.1,acc=1e-2,eps=1e-2;
	int n=2;
	gsl_vector* yx=gsl_vector_alloc(n);
	gsl_vector_set(yx,0,0);
	gsl_vector_set(yx,1,1);
	fprintf(stderr,"# index 0: u''=-u\n");
	driver(f1,a,yx,b,h,acc,eps);
	fprintf(stderr,"\n\n");

	double b1=14;
	int m=3;
	gsl_vector* yx2=gsl_vector_alloc(m);
	gsl_vector_set(yx2,0,5998900);
	gsl_vector_set(yx2,1,100000);
	gsl_vector_set(yx2,2,100);
	fprintf(stderr,"# index 1: Tc=2\n");
	driver(f2,a,yx2,b1,h,acc,eps);
	fprintf(stderr,"\n\n");

	gsl_vector* yx4=gsl_vector_alloc(m);
	gsl_vector_set(yx4,0,5998900);
	gsl_vector_set(yx4,1,100000);
	gsl_vector_set(yx4,2,100);
	fprintf(stderr,"# index 2: Tc=4\n");
	driver(f4,a,yx4,b1,h,acc,eps);
	fprintf(stderr,"\n\n");

	gsl_vector* yx8=gsl_vector_alloc(m);
	gsl_vector_set(yx8,0,5998900);
	gsl_vector_set(yx8,1,100000);
	gsl_vector_set(yx8,2,100);
	fprintf(stderr,"# index 3: Tc=8\n");
	driver(f8,a,yx8,b1,h,acc,eps);
	fprintf(stderr,"\n\n");

gsl_vector_free(yx);
gsl_vector_free(yx4);
gsl_vector_free(yx2);
gsl_vector_free(yx8);
return 0;
}
