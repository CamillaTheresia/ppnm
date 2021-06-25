#include<math.h>
#include<stdio.h>
#define FOR(i) for(int i=0;i<n;i++)

void midpoint(
	int n,
	void f(int n,double t,double y[],double dydt[]),
	double t, double y[], double h,
	double yh[], double dy[]
){
	double k0[n],kmid[n],ymid[n];
	f(n,t,y,k0);
	FOR(i) ymid[i]=y[i]+(0.5 *h)*k0[i];
	f(n,t+0.5*h,ymid,kmid);
	FOR(i) yh[i]=y[i]+h*kmid[i];
	FOR(i) dy[i]=(kmid[i]-k0[i])*h;
}

#define PRINT(x) fprintf(stderr,"%9.3g ",x)
//#define TRACE(t,y) PRINT(t);FOR(i)PRINT(y[i]);fprintf(stderr,"\n")
#define TRACE(t,y)

int odedriver(
	int  n, /* y[n] */
	void f(int n,double t,double*y,double*dydt), /* dy/dt=f(t,y) */
	double a,              /* the start-point a */
	double*y,                    /* y(a) -> y(b) */
	double b,              /* the end-point of the integration */
	double h,                    /* initial step-size */
	double acc,            /* absolute accuracy goal */
	double eps             /* relative accuracy goal */
){
	int steps=0;
	double t=a;
	TRACE(t,y);
	while(t<b){
		if(t+h>b)h=b-t;
		double yh[n],dy[n];
		midpoint(n,f,t,y,h,yh,dy);
		double sum=0; FOR(i)sum+=y[i]*y[i];
		double norm_y=sqrt(sum);
		sum=0; FOR(i)sum+=dy[i]*dy[i];
		double err=sqrt(sum);
		double tol=(acc+eps*norm_y)*sqrt(h/(b-a));
		if(err<tol){
			steps++;
			t=t+h;
			FOR(i) y[i]=yh[i];
			TRACE(t,y);
		}
		if(err>0) h*=0.95*pow(tol/err,0.25);
		else h*=2;
	}//while
return steps;
}
