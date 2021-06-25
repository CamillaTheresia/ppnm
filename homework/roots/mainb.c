#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>

double Fe(double e,double r);
int newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps);

static double rmax;
void master(gsl_vector*x,gsl_vector*fx){
	double e=gsl_vector_get(x,0);
	assert(e<0);
	double frmax=Fe(e,rmax);
	gsl_vector_set(fx,0,frmax);
	}

int main(int argc, char** argv){
	rmax = argc>1 ? atof(argv[1]):4;
	int dim=1;
	gsl_vector* x=gsl_vector_alloc(dim);
	gsl_vector_set(x,0,-1);
	newton(master,x,1e-3);
	double e=gsl_vector_get(x,0);

	printf("# rmax, e=\n");
	printf("%g %g\n",rmax,e);
	printf("\n\n");

	printf("# r, Fe(e,r), exact\n");
	for(double r=0;r<=rmax;r+=rmax/64)
		printf("%g %g %g\n",r,Fe(e,r),r*exp(-r));
	printf("\n\n");
return 0;
}
