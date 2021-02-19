#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_gamma.h>
#include"gamma.h"

int main(){
	double xmin=1,xmax=8;
	for(double x=xmin;x<=xmax;x+=1.0/8){
		printf("%10g %10g %10g %10g\n",x,tgammaf(x),gsl_sf_gamma(x),Gamma(x));
		}
return 0;
}
