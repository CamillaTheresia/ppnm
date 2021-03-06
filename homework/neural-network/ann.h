#include<gsl/gsl_vector.h>
#ifndef HAVE_NEURONS
#define HAVE_NEURONS
typedef struct {
	int n;
	double(*f)(double);
	double(*d)(double);
	double(*i)(double);
	gsl_vector* params;
	} ann;
ann* ann_alloc(int n,double(*f)(double),double(*d)(double),double(*i)(double));
void ann_free(ann* network);
double ann_response(ann* network,double x);
double ann_derivative(ann* network,double x);
double ann_integral(ann* network,double x);
void ann_train(double cost_function(gsl_vector* p),ann* network,gsl_vector* xlist,gsl_vector* ylist);
#endif
