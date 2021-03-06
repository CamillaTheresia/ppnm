#ifndef HAVE_ADAPT_H
#define HAVE_ADAPT_H
double adapt24 (
	double f(double),double a, double b,
	double acc, double eps, double f2, double f3, int nrec
	);
double adapt (
	double f(double),double a,double b,
	double acc,double eps
	);
double clenshaw_curtis(
	double f(double),double a,double b,double acc,double eps
	);
#endif
