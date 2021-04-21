#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<complex.h>

double corput(int n, int base){
	double q=0, bk=(double)1/base;
	while(n>0) { q += (n % base)*bk; n /= base; bk /= base; }
	return q;
}

void halton1(int n, int d, double *a, double *b, double *x){
	int base[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79};
	int dmax=sizeof(base)/sizeof(int); assert(d <= dmax);
	for(int i=0;i<d;i++)x[i]=a[i]+corput(n+1,base[i])*(b[i]-a[i]);
}

void halton2(int n, int d, double *a, double *b, double *x){
	int base[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83};
	int dmax=sizeof(base)/sizeof(int); assert(d <= dmax);
	for(int i=0;i<d;i++)x[i]=a[i]+corput(n+1,base[i])*(b[i]-a[i]);
}

complex quasimc(int d, double f(int,double*), double a[], double b[], int N){
	double x[d];
	double vol=1; for(int i=0;i<d;i++) vol*=b[i]-a[i];
	double sum1=0, sum2=0;
	for(int i=0;i<N/2;i++){
		halton1(i,d,a,b,x); sum1+=f(d,x);
		}
	for(int i=0;i<N/2;i++){
		halton2(i,d,a,b,x); sum2+=f(d,x);
		}
	double integ=(sum1+sum2)/N*vol;
	double error=fabs(sum1-sum2)/N*vol;
	return integ+I*error;
	}
