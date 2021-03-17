#include<stdlib.h>
#include<assert.h>
#include<stdio.h>

typedef struct {int n; double *x, *y, *b, *c;} qspline;

qspline* qspline_alloc(int n, double *x, double *y){
	qspline *s=(qspline*)malloc(sizeof(qspline));
	s->b = (double*)malloc((n-1)*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));
	s->n = n;
	for(int i=0;i<n;i++){s->x[i]=x[i]; s->y[i]=y[i];}
	int i; double p[n-1], h[n-1];
	for(i=0;i<n-1;i++){h[i]=x[i+1]-x[i]; p[i]=(y[i+1]-y[i])/h[i];}
	s->c[0]=0;
	for(i=0;i<n-2;i++)s->c[i+1]=(p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
	s->c[n-2]/=2;
	for(i=n-3;i>=0;i--)s->c[i]=(p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
	for(i=0;i<n-1;i++)s->b[i]=p[i]-s->c[i]*h[i];
	return s;
}

double qspline_eval(qspline *s,double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i=0, j=s->n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>s->x[mid]) i=mid; else j=mid;
		}
	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*s->c[i]);
}

double qspline_integ_eval(qspline *s,double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i=0, j=s->n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>s->x[mid]) i=mid; else j=mid;
		}
	double h=z-s->x[i];
	double w=s->y[i]*h+0.5*h*h*s->b[i]+1./3*h*h*h*s->c[i];
	for(int k=0;k<i;k++){
		w+=s->y[k]*(s->x[k+1]-s->x[k])+0.5*(s->x[k+1]-s->x[k])*(s->x[k+1]-s->x[k])*s->b[k]+1./3*(s->x[k+1]-s->x[k])*(s->x[k+1]-s->x[k])*(s->x[k+1]-s->x[k])*s->c[k];
	}
	return w;
}

double qspline_diff_eval(qspline *s,double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i=0, j=s->n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>s->x[mid]) i=mid; else j=mid;
		}
	double h=z-s->x[i];
	return h*s->b[i]+2*h*s->c[i];
}

void qspline_free(qspline *s){
	free(s->x); free(s->y); free(s->b); free(s->c); free(s);
}

int main(){
	int n=5;
	double x[n],y1[n],y2[n],y3[n],y1_diff[n],y2_diff[n],y3_diff[n],y1_integ[n],y2_integ[n],y3_integ[n];
	printf("x: y1: y2: y3:\n");
	for(int i=0;i<n;i++){
		x[i]=i+1;
		y1[i]=1;
		y2[i]=i+1;
		y3[i]=(i+1)*(i+1);
		printf("%10g %10g %10g %10g\n\n",x[i],y1[i],y2[i],y3[i]);
	}
	qspline* s1=qspline_alloc(n,x,y1);
	qspline* s2=qspline_alloc(n,x,y2);
	qspline* s3=qspline_alloc(n,x,y3);
	printf("z: qspline(y1): qspline(y2): qspline(y3):\n");
	double zmin=2,zmax=4;
	for(double z=zmin;z<=zmax;z+=1./2){
		printf("%10g %10g %10g %10g\n\n",z,qspline_eval(s1,z),qspline_eval(s2,z),qspline_eval(s3,z));
	}
	printf("x: y1_diff: y2_diff: y3_diff:\n");
	for(int i=0;i<n;i++){
		x[i]=i+1;
		y1_diff[i]=0;
		y2_diff[i]=1;
		y3_diff[i]=2*(i+1);
		printf("%10g %10g %10g %10g\n\n",x[i],y1_diff[i],y2_diff[i],y3_diff[i]);
	}
	printf("z: qspline_diff(y1): qspline_diff(y2): qspline_diff(y3):\n");
	for(double z=zmin;z<=zmax;z+=1./2){
		printf("%10g %10g %10g %10g\n\n",z,qspline_diff_eval(s1,z),qspline_diff_eval(s2,z),qspline_diff_eval(s3,z));
	}
	printf("x: integ(y1): integ(y2): integ(y3):\n");
	for(int i=0;i<n;i++){
		x[i]=i+1;
		y1_integ[i]=(i+1)-1;
		y2_integ[i]=0.5*(i+1)*(i+1)-0.5;
		y3_integ[i]=1./3*(i+1)*(i+1)*(i+1)-1./3;
		printf("%10g %10g %10g %10g\n\n",x[i],y1_integ[i],y2_integ[i],y3_integ[i]);
	}
	printf("z: qspline_integ(y1): qspline_integ(y2): qspline_integ(y3):\n");
	for(double z=zmin;z<=zmax;z+=1./2){
		printf("%10g %10g %10g %10g\n\n",z,qspline_integ_eval(s1,z),qspline_integ_eval(s2,z),qspline_integ_eval(s3,z));
	}
	printf("b1: b2: b3:\n");
	for(int i=0;i<n-1;i++){
		printf("%10g %10g %10g\n\n",(double)s1->b[i],(double)s2->b[i],(double)s3->b[i]);
	}
	printf("c1: c2: c3:\n");
	for(int i=0;i<n-1;i++){
		printf("%10g %10g %10g\n\n",(double)s1->c[i],(double)s2->c[i],(double)s3->c[i]);
	}
	qspline_free(s1);
	qspline_free(s2);
	qspline_free(s3);
return 0;
}


