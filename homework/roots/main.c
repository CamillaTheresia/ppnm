#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<float.h>
#include<math.h>

int newton (
	void f(gsl_vector* x,gsl_vector* fx),
	gsl_vector* x, double eps);
void rkstep12(
	void f(double x, gsl_vector* yx, gsl_vector* dydx),
	double x, gsl_vector* yx, double h, gsl_vector* yh, gsl_vector* err,gsl_vector* k0, gsl_vector* k1
);

void vector_print(char* s,gsl_vector* v){
	printf("%s",s);
	for(int i=0;i<v->size;i++)printf("%10.4g ",gsl_vector_get(v,i));
	printf("\n");
	}

static int ncalls;
void ft(gsl_vector* p, gsl_vector* fx){ /* quadratic test function */
	ncalls++;
	double x=gsl_vector_get(p,0);
	gsl_vector_set(fx,0,3*2*x+2);
	}

void f(gsl_vector* p,gsl_vector* fx){
	ncalls++;
	double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
	gsl_vector_set(fx,0, 2*(1-x)*(-1)+100*2*(y-x*x)*(-1)*2*x);
	gsl_vector_set(fx,1, 100*2*(y-x*x));
	}

int main() {
	gsl_vector* xt=gsl_vector_alloc(1);
	gsl_vector_set(xt,0,-3);
	gsl_vector* fxt=gsl_vector_alloc(1);
	gsl_vector* x=gsl_vector_alloc(2);
	gsl_vector_set(x,0,-2);
	gsl_vector_set(x,1,8);
	gsl_vector* fx=gsl_vector_alloc(2);

	printf("Root finding:\n");
	printf("Test function: Find minimum of 3*x^2+2*x+1:\n");
	vector_print("initial vector x: ",xt);
	ft(xt,fxt);
	vector_print("            f(x): ",fxt);
	ncalls=0;
	int stepst=newton(ft,xt,1e-3);
	printf("steps = %i,	ncalls = %i\n",stepst,ncalls);
	vector_print("      solution x: ",xt);
	ft(xt,fxt);
	vector_print("            f(x): ",fxt);
	printf("Minimum of 3*x^2+2x*+1 should be:\n");
	printf("x=-0.3333\n");

	printf("Extremum of the Rosenbrock's function:\n");
	vector_print("initial vector x: ",x);
	f(x,fx);
	vector_print("            f(x): ",fx);
	ncalls=0;
	int steps=newton(f,x,1e-3);
	printf("steps = %i,	ncalls = %i\n",steps,ncalls);
	vector_print("      solution x: ",x);
	f(x,fx);
	vector_print("            f(x): ",fx);
	printf("Minimum of Rosenbrock's function should be:\n");
	printf("x=1	y=1\n");

return 0;
}
