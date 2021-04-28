#include<math.h>
#include<assert.h>
#include<stdio.h>
#include<gsl/gsl_integration.h>
#include"adapt.h"
#define SQR2 1.41421356237309504880

int calls;
double f(double x){calls++; return 1/sqrt(x);};
double f2(double x){calls++; return log(x)/sqrt(x);};
double fp(double x){calls++; return 4*sqrt(1-x*x)/M_PI;};

double fg(double x, void * params){
	calls++;
	double alpha = *(double *) params;
	double f = alpha/sqrt(x);
	return f;
}

double f2g(double x, void * params){
	calls++;
	double alpha = *(double *) params;
	double f = log(alpha*x)/sqrt(x);
	return f;
}

double fpg(double x, void * params){
	calls++;
	double alpha = *(double *) params;
	double f = alpha*4*sqrt(1-x*x)/M_PI;
	return f;
}

int main()
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);
	double result, error;
 	double alpha = 1.0;
	gsl_function F;
	F.function = &fg;
	F.params = &alpha;

	double a=0,b=1,acc=0.001,eps=0.001;
	calls=0;
	double Q=adapt(f,a,b,acc,eps);
	double exact=2;
	printf("Integration of 1/sqrt(x) from %g to %g : OPEN4\n",a,b);
	printf("              Q = %g\n",Q);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n\n",fabs(Q-exact));

	calls=0;
	Q=clenshaw_curtis(f,a,b,acc,eps);
	exact=2;
	printf("Integration of 1/sqrt(x) from %g to %g : CLENSHAW_CURTIS\n",a,b);
	printf("              Q = %g\n",Q);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n\n",fabs(Q-exact));

	calls=0;
	gsl_integration_qags(&F, 0, 1, acc, eps, 10000,w, &result, &error);
	exact=2;
	printf("Integration of 1/sqrt(x) from %g to %g : GSL\n",a,b);
	printf("              Q = %g\n",result);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",error);
	printf("   actual error = %g\n\n",result-exact);

	a=0,b=1,acc=0.001,eps=0.001;

	double result1, error1;
 	double alpha1 = 1.0;
	gsl_function Fone;
	Fone.function = &f2g;
	Fone.params = &alpha1;

	calls=0;
	Q=adapt(f2,a,b,acc,eps);
	exact=-4;
	printf("Integration of log(x)/sqrt(x) from %g to %g : OPEN4\n",a,b);
	printf("              Q = %g\n",Q);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n\n",fabs(Q-exact));

	calls=0;
	Q=clenshaw_curtis(f2,a,b,acc,eps);
	exact=-4;
	printf("Integration of log(x)/sqrt(x) from %g to %g : CLENSHAW_CURTIS\n",a,b);
	printf("              Q = %g\n",Q);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n\n",fabs(Q-exact));

	calls=0;
	gsl_integration_qags(&Fone, 0, 1, acc, eps, 10000,w, &result1, &error1);
	exact=-4;
	printf("Integration of log(x)/sqrt(x) from %g to %g : GSL\n",a,b);
	printf("              Q = %g\n",result1);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",error1);
	printf("   actual error = %g\n\n",result1-exact);

	a=0,b=1,acc=1e-9,eps=1e-9;

	double result2, error2;
 	double alpha2 = 1.0;
	gsl_function Ftwo;
	Ftwo.function = &fpg;
	Ftwo.params = &alpha2;

	calls=0; Q=adapt(fp,a,b,acc,eps);
	exact=1;
	printf("Integration of 4*sqrt(1-x*x)/M_PI from %g to %g : OPEN4\n",a,b);
	printf("              Q = %25.20f\n",Q);
	printf("          exact = %25.20f\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n\n",fabs(Q-exact));

	calls=0; Q=clenshaw_curtis(fp,a,b,acc,eps);
	exact=1;
	printf("Integration of 4*sqrt(1-x*x)/M_PI from %g to %g : CLENSHAW_CURTIS\n",a,b);
	printf("              Q = %25.20f\n",Q);
	printf("          exact = %25.20f\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n\n",fabs(Q-exact));

	calls=0;
	gsl_integration_qags(&Ftwo, 0, 1, acc, eps, 10000,w, &result2, &error2);
	exact=1;
	printf("Integration of 4*sqrt(1-x*x)/M_PI from %g to %g : GSL\n",a,b);
	printf("              Q = %g\n",result2);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",error2);
	printf("   actual error = %g\n\n",result2-exact);

gsl_integration_workspace_free(w);
return 0 ;
}
