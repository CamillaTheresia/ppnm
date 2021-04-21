#include<stdio.h>
#include<math.h>
#define R 0.8

double strata(
        int dim,
        double f(int dim,double*x),
        double*a,double*b,
        double acc,double eps,int nreuse,double mreuse);

double f(int dim,double*p){
	double x=p[0],y=p[1];
fprintf(stderr,"%g %g\n",x,y);
	double r=x*x+y*y<R*R?1:0;
	return r;
	}

int main(){
	double a[]={0,0},b[]={1,1},acc=1e-3,eps=0;
	int dim=sizeof(a)/sizeof(a[0]);
	double integ=strata(dim,f,a,b,acc,eps,0,0);
	double exact=M_PI*R*R/4;
	printf("strata integration f(x,y)=x*x+y*y<R*R?1:0, R=%g\n",R);
	printf("a = {%g,%g} b = {%g,%g}\n",a[0],a[1],b[0],b[1]);
	printf("acc = %g eps = %g\n\n",acc,eps);
	printf("integ = %g error-estimate = %g\n",integ,acc+fabs(integ)*eps);
	printf("exact = %g actual-error   = %g\n",exact,fabs(integ-exact));
}
