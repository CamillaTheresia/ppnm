#include<stdio.h>
#include<limits.h>
#include<float.h>

int main(){
	int i=1; while(i+1>i) {i++;}
	int v;
	for(int j=1;j+1>j;j++) {v=j;};
	int k=1;do{
		k++;
		}while(k+1>k);
	printf("my max int with while =%i\n",i);
	printf("my max int with for =%i\n",v);
	printf("my max int with do while =%i\n",k);
	printf("value of INT_MAX =%i\n",INT_MAX);

	int a=1; while(a-1<a) {a=a-1;}
	int w;
	for(int b=1;b-1<b;b=b-1) {w=b;};
	int c=1;do{
		c=c-1;
		}while(c-1<c);
	printf("my min int with while =%i\n",a);
	printf("my min int with for =%i\n",w);
	printf("my min int with do while =%i\n",c);
	printf("value of INT_MAX =%i\n",INT_MIN);

	float x=1; while(1+x!=1) {x/=2;} x*=2;
	double  x_d=1; while(1+x_d!=1) {x_d/=2;} x_d*=2;
	long double x_ld=1; while(1+x_ld!=1) {x_ld/=2;} x_ld*=2;
	float y; for(y=1;1+y!=1;y/=2){} y*=2;
	double y_d; for(y_d=1;1+y_d!=1;y_d/=2){} y_d*=2;
	long double y_ld; for(y_ld=1;1+y_ld!=1;y_ld/=2){} y_ld*=2;
	float z=1; do{z/=2;}while(1+z!=1); z*=2;
	double z_d=1; do{z_d/=2;}while(1+z_d!=1); z_d*=2;
	long double z_ld=1; do{z_ld/=2;}while(1+z_ld!=1); z_ld*=2;
	printf("float epsilon with while =%g\n",x);
	printf("double epsilon with while =%lg\n",x_d);
	printf("long double epsilon with while =%Lg\n",x_ld);
	printf("float epsilon with for =%g\n",y);
	printf("double epsilon with for =%lg\n",y_d);
	printf("long double epsilon with for =%Lg\n",y_ld);
	printf("float epsilon with do while=%g\n",z);
	printf("double epsilon with do while =%lg\n",z_d);
	printf("long double epsilon with do while =%Lg\n",z_ld);
	printf("FLT_EPSILON =%g\n",FLT_EPSILON);
	printf("DBL_EPSILON =%g\n",DBL_EPSILON);
	printf("LDBL_EPSILON =%Lg\n",LDBL_EPSILON);

	int max=INT_MAX/2;
	float sum_up_float=0;
	double sum_up_double=0;
	for(int d=1;d<=max;d++){
		sum_up_float+=1.0/d;
		sum_up_double+=1.0/d;
	}
	printf("sum_up_float=%g\n",sum_up_float);
	printf("sum_up_double=%g\n",sum_up_double);
	float sum_down_float=0;
	double sum_down_double=0;
	for(int e=max;e>=1;e--){
		sum_down_float+=1.0/e;
		sum_down_double+=1.0/e;
	}
	printf("sum_down_float=%g\n",sum_down_float);
	printf("sum_down_double=%g\n",sum_down_double);
return 0;
}
