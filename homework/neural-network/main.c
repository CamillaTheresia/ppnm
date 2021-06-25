#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include"ann.h"
double activation_function(double x){return x*exp(-x*x);}
double function_to_fit(double x){return cos(5*x-1)*exp(-x*x);}
double derivative(double x){return exp(-x*x)-2*x*x*exp(-x*x);}
double integral(double x){return -0.5*exp(-x*x);}

gsl_vector* vx;
gsl_vector* vy;
ann* network;
double cost_function(gsl_vector* p){
	gsl_vector_memcpy(network->params,p);
	double sum=0;
	for(int i=0;i<vx->size;i++){
		double xi=gsl_vector_get(vx,i);
		double yi=gsl_vector_get(vy,i);
		double fi=ann_response(network,xi);
		sum+=fabs(fi-yi);
	}
	return sum/vx->size;
}

int main(){
	int n=5;
	network=ann_alloc(n,activation_function,derivative,integral);
	double a=-1,b=1;
	int nx=20;
	vx=gsl_vector_alloc(nx);
	vy=gsl_vector_alloc(nx);

	printf("# Function the neural-network should fit to:\n");
	for(int i=0;i<nx;i++){
		double x=a+(b-a)*i/(nx-1);
		double f=function_to_fit(x);
		gsl_vector_set(vx,i,x);
		gsl_vector_set(vy,i,f);
		printf("%g %g \n",x,f);
	}
	printf("\n \n");

	for(int i=0;i<network->n;i++){
		gsl_vector_set(network->params,3*i+0,a+(b-a)*i/(network->n-1));
		gsl_vector_set(network->params,3*i+1,1);
		gsl_vector_set(network->params,3*i+2,1);
//		printf("%g %g %g \n",gsl_vector_get(network->params,3*i+0),gsl_vector_get(network->params,3*i+1),gsl_vector_get(network->params,3*i+2));
	}
//	printf("\n\n");
	ann_train(cost_function,network,vx,vy);
//	for(int i=0;i<network->n;i++){
//		printf("%g %g %g \n",gsl_vector_get(network->params,3*i+0),gsl_vector_get(network->params,3*i+1),gsl_vector_get(network->params,3*i+2));
//	}
//	printf("\n\n");

	printf("# The response of the neural-network:\n");
	double dz=1.0/64;
	for(double z=a;z<=b;z+=dz){
		printf("%g %g %g %g\n",z,ann_response(network,z),ann_derivative(network,z),ann_integral(network,z));
	}

gsl_vector_free(vx);
gsl_vector_free(vy);
ann_free(network);
return 0;
}
