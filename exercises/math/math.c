#include<math.h>
#include<complex.h>
#include<stdio.h>

int main(){
	double x=tgamma(5);
	double y=jn(1,0.5);
	complex z=csqrt(-2);
	complex a=cexp(I*M_PI);
	complex b=cexp(I);
	complex c=cpow(I,M_E);
	complex d=cpow(I,I);
	printf("tgamma(5)=%g+I%g\n",creal(x),cimag(x));
	printf("j1(0.5)=%g+I%g\n",creal(y),cimag(y));
	printf("sqrt(-2)=%g+I%g\n",creal(z),cimag(z));
	printf("exp(I*pi)=%g+I%g\n",creal(a),cimag(a));
	printf("exp(I)=%g+I%g\n",creal(b),cimag(b));
	printf("pow(I,e)=%g+I%g\n",creal(c),cimag(c));
	printf("pow(I,I)=%g+I%g\n",creal(d),cimag(d));
return 0;
}
