#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

int qnewton(double beta(gsl_vector*),gsl_vector* x,double acc);


void vector_print(char* s,gsl_vector* v){
	printf("%s",s);
	for(int i=0;i<v->size;i++)printf("%10.4g ",gsl_vector_get(v,i));
	printf("\n");
	}

double fr(gsl_vector* p){
	double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
	double fx = (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
	return fx;
	}

double fh(gsl_vector* p){
	double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
	double fx = (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
	return fx;
	}

gsl_vector* E;
gsl_vector* sigma;
gsl_vector* dsigma;
double chi2(gsl_vector* p){
	double A=gsl_vector_get(p,0), m=gsl_vector_get(p,1), L=gsl_vector_get(p,2), sum=0, F;
	for(int i=0;i<30;i++){
		F=A/((gsl_vector_get(E,i)-m)*(gsl_vector_get(E,i)-m)+L*L/4);
		sum+=(F-gsl_vector_get(sigma,i))*(F-gsl_vector_get(sigma,i))/(gsl_vector_get(dsigma,i)*gsl_vector_get(dsigma,i));
	}
	return sum;
}

int main(){
	gsl_vector* xr=gsl_vector_alloc(2);
	gsl_vector_set(xr,0,-2);
	gsl_vector_set(xr,1,8);
	double accr=1e-6;
	fprintf(stderr,"minimum of (1-x)^2+100*(y-x)^2\n");
	qnewton(fr,xr,accr);
	fprintf(stderr,"calculated: x=%g, y=%g\n",gsl_vector_get(xr,0),gsl_vector_get(xr,1));
	fprintf(stderr,"should be: x=1, y=1\n");
	fprintf(stderr,"\n");

	gsl_vector* xh=gsl_vector_alloc(2);
	gsl_vector_set(xh,0,-2);
	gsl_vector_set(xh,1,8);
	double acch=1e-6;
	fprintf(stderr,"minimum of (x^2+y-11)^2+(x+y^2-7)^2\n");
	qnewton(fh,xh,acch);
	fprintf(stderr,"calculated: x=%g, y=%g\n",gsl_vector_get(xh,0),gsl_vector_get(xh,1));
	fprintf(stderr,"should be: x=-3.779310, y=-3.283186\n");
	fprintf(stderr,"\n");

	E=gsl_vector_alloc(30);
	sigma=gsl_vector_alloc(30);
	dsigma=gsl_vector_alloc(30);
	double En, sigman, dsigman;
	FILE * data;
	data = fopen("data.txt","r");
	for(int i=0;i<30;i++){
		fscanf(data,"%lg %lg %lg",&En,&sigman,&dsigman);
/*		fprintf(stderr,"%lg %lg %lg\n",En,sigman,dsigman); */
		gsl_vector_set(E,i,En);
		gsl_vector_set(sigma,i,sigman);
		gsl_vector_set(dsigma,i,dsigman);
		}
	fclose(data);
	gsl_vector* H=gsl_vector_alloc(3);
	gsl_vector_set(H,0,2);
	gsl_vector_set(H,1,130);
	gsl_vector_set(H,2,1);
	double accH=1e-6;
	qnewton(chi2,H,accH);
	fprintf(stderr,"calculated: A=%g, m=%g, gamma=%g\n",gsl_vector_get(H,0),gsl_vector_get(H,1),gsl_vector_get(H,2));
	fprintf(stderr,"should be: m=125.3\n");
	fprintf(stderr,"\n");

	for(double j=101;j<=159;j+=1.0/8){
	printf("%g %g\n",j,gsl_vector_get(H,0)/((j-gsl_vector_get(H,1))*(j-gsl_vector_get(H,1))+gsl_vector_get(H,2)*gsl_vector_get(H,2)/4));
	}

gsl_vector_free(xr);
gsl_vector_free(xh);
gsl_vector_free(E);
gsl_vector_free(sigma);
gsl_vector_free(dsigma);
gsl_vector_free(H);
return 0;
}
