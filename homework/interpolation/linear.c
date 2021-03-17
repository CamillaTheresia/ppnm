#include<stdio.h>
#include<assert.h>
#include <gsl/gsl_interp.h>

double linerp(int n, double x[], double y[], double z){
	assert(x[0]<=z && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	assert(x[i+1]>x[i]);
	return y[i]+(y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
}

double linerp_integ(int n, double x[], double y[], double z){
	assert(x[0]<=z && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	assert(x[i+1]>x[i]);
	return 0.5*(y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i])*(z-x[i])+y[i]*(z-x[i])+0.5*3*(x[i]-x[0])*(x[i]-x[0])+y[0]*(x[i]-x[0])+y[0]*(x[i]-x[0]);
}

int main(){
	int n=10;
	double x[n];
	double y[n];
	double y_integ[n];
	FILE* xyout=fopen("out.linear.xy.txt","w");
	for(int i=0;i<n;i++){
		x[i]=i;
		y[i]=3*i;
		y_integ[i]=0.5*3*i*i;
		fprintf(xyout,"%10g %10g %10g\n",x[i],y[i],y_integ[i]);
	}
	fclose(xyout);

	gsl_interp* linear = gsl_interp_alloc(gsl_interp_linear,n);
	gsl_interp_init(linear,x,y,n);
	double zmin=1,zmax=8;
	for(double z=zmin;z<zmax;z+=1./8){
		double interp_l=gsl_interp_eval(linear ,x,y,z,NULL);
		double integ_l=gsl_interp_eval_integ(linear,x,y,x[0],z,NULL);
		printf("%10g %10g %10g %10g %10g\n",z,linerp(n,x,y,z),linerp_integ(n,x,y,z),interp_l,integ_l);
	}
gsl_interp_free(linear);
return 0;
}
