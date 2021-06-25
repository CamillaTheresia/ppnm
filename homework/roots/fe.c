#include<assert.h>
#include<stdio.h>
#include"ode.h"

static double E;
void equation(int n,double r,double* y,double* dydr){
	assert(n==2);
	dydr[0]=y[1];
	dydr[1]=2*(-1/r-E)*y[0]; /* -(1/2)f'' - (1/r)f = e f */
	}

double Fe(double e,double r){
	assert(r>=0);
	E=e;
	const double rmin = 1e-3;
	if(r<rmin) return r-r*r;
	double y[2]={rmin-rmin*rmin,1-2*rmin};
	int n=2;
	double h=0.1,acc=1e-3,eps=1e-3;
	odedriver(n,equation,rmin,y,r,h,acc,eps);
	return y[0];
}
