#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>

complex plainmc(int dim,double f(int dim,double*x),double*a,double*b,int N);
complex quasimc(int dim,double f(int dim,double*x),double*a,double*b,int N);

#define R 0.9
double f(int dim,double* p){
	assert(dim==2);
	double x=p[0],y=p[1];
	if(x*x+y*y<R*R)return 1;
	else return 0;
	}

double f1(int dim,double* p){
	assert(dim==3);
	double x=p[0],y=p[1],z=p[2];
	double a = 1/pow(M_PI,3)*1/(1-cos(x)*cos(y)*cos(z));
	return a;
	}

double f2(int dim, double* p){
	assert(dim==1);
	double x=p[0];
	double a=sin(x)*sin(x)*sin(x)/(sin(x)*sin(x)*sin(x)+cos(x)*cos(x)*cos(x));
	return a;
}

int main(int argc,char** argv){
double a[]={0,0},b[]={1,1};
int dim=sizeof(a)/sizeof(a[0]);
if(argc>1){
	srand(42);
	int N=(int)atof(argv[1]);
	complex result_p=plainmc(dim,f,a,b,N);
	complex result_q=quasimc(dim,f,a,b,N);
	double integ_p=creal(result_p);
	double integ_q=creal(result_q);
	double exact=M_PI/4*R*R;
	double error_p=fabs(integ_p-exact);
	double error_q=fabs(integ_q-exact);
	printf("%i %g %g\n",N,error_p,error_q);
	}

else{
	int N=(int)1e4;
	complex result_p=plainmc(dim,f,a,b,N);
	complex result_q=quasimc(dim,f,a,b,N);
	double integ_p=creal(result_p), esterr_p=cimag(result_p);
	double integ_q=creal(result_q), esterr_q=cimag(result_q);
	double exact=M_PI/4*R*R;
	double error_p=fabs(integ_p-exact);
	double error_q=fabs(integ_q-exact);
	printf("integrating x*x+y*y<%g*%g?1:0 with N=%i\n",R,R,N);
	printf("pseudo-random:\n");
	printf("integral = %f error estimate = %f\n",integ_p,esterr_p);
	printf("exact    = %f actual error   = %f\n",exact,error_p);
	printf("quasi-random:\n");
	printf("integral = %f error estimate = %f\n",integ_q,esterr_q);
	printf("exact    = %f actual error   = %f\n",exact,error_q);

	double a1[]={0,0,0},b1[]={M_PI,M_PI,M_PI};
	int dim1=sizeof(a1)/sizeof(a1[0]);
	int M=(int)1e6;
	complex result_p1=plainmc(dim1,f1,a1,b1,M);
	double integ_p1=creal(result_p1), esterr_p1=cimag(result_p1);
	double exact1=1.3932039296856768591842462603255;
	double error_p1=fabs(integ_p1-exact1);
	printf("integrating 1/pow(M_PI,3)*1/(1-cos(x)*cos(y)*cos(z))dxdydz from 0 to pi with N=%i\n",M);
	printf("pseudo-random:\n");
	printf("integral = %f error estimate = %f\n",integ_p1,esterr_p1);
	printf("exact    = %f actual error   = %f\n",exact1,error_p1);

	double a2[]={0},b2[]={M_PI/2};
	int dim2=sizeof(a2)/sizeof(a2[0]);
	int O=(int)1e4;
	complex result_p2=plainmc(dim2,f2,a2,b2,O);
	double integ_p2=creal(result_p2), esterr_p2=cimag(result_p2);
	double exact2=M_PI/4;
	double error_p2=fabs(integ_p2-exact2);
	printf("integrating sin(x)*sin(x)*sin(x)/(sin(x)*sin(x)*sin(x)+cos(x)*cos(x)*cos(x)) from 0 to pi/2 with N=%i\n",O);
	printf("pseudo-random:\n");
	printf("integral = %f error estimate = %f\n",integ_p2,esterr_p2);
	printf("exact    = %f actual error   = %f\n",exact2,error_p2);
	}

return 0;
}
