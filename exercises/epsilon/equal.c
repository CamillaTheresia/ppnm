#include<stdio.h>
#include<math.h>

int equal(double a, double b, double tau, double epsilon){
	if (fabs(a-b)<tau || fabs(a-b)/(fabs(a)+fabs(b))<epsilon/2) {
		return 1;
	}
	else {
		return 0;
	}
}

int main(){
	double a=1;
	double b=1.5;
	double tau=1;
	double epsilon=0.5;
	int result=equal(a,b,tau,epsilon);
	printf("a=%g\n",a);
	printf("b=%g\n",b);
	printf("tau=%g\n",tau);
	printf("epsilon=%g\n",epsilon);
	printf("result=%d\n",result);
return 0;
}
