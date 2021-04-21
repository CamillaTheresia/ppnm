#include<math.h>
#include<assert.h>
#include<stdio.h>
#include"adapt.h"
#define SQR2 1.41421356237309504880

int calls;
double f(double x){calls++; return 1/sqrt(x);};
double f2(double x){calls++; return log(x)/sqrt(x);};
double fp(double x){calls++; return 4*sqrt(1-x*x)/M_PI;};

int main()
{
	double a=0,b=1,acc=0.001,eps=0.001;
	calls=0;
	double Q=adapt(f,a,b,acc,eps);
	double exact=2;
	printf("Integration of 1/sqrt(x) from %g to %g : OPEN4\n",a,b);
	printf("              Q = %g\n",Q);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n",fabs(Q-exact));

	a=0,b=1,acc=0.001,eps=0.001;

	calls=0;
	Q=adapt(f2,a,b,acc,eps);
	exact=-4;
	printf("Integration of log(x)/sqrt(x) from %g to %g : OPEN4\n",a,b);
	printf("              Q = %g\n",Q);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n",fabs(Q-exact));

	a=0,b=1,acc=1e-9,eps=1e-9;

	calls=0; Q=adapt(fp,a,b,acc,eps);
	exact=1;
	printf("Integration of 4*sqrt(1-x*x)/M_PI from %g to %g : OPEN4\n",a,b);
	printf("              Q = %25.20f\n",Q);
	printf("          exact = %25.20f\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n",fabs(Q-exact));

return 0 ;
}
