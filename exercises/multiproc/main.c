#include<stdio.h>
#include<pthread.h>
#include<math.h>
#include<stdlib.h>

struct params {int n, seed; int count;};

void* bar(void* arg){
	struct params * p = (struct params *)arg;
	int N = (*p).n;
	unsigned int seed = (unsigned int)(*p).seed;
	int c=0;
	double x,y,z;
	for(int i=0;i<N;i++) {
		x = (double)rand_r(&seed)/RAND_MAX;
		y = (double)rand_r(&seed)/RAND_MAX;
		z = x*x+y*y;
		if(z<=1) c++;
	}
	(*p).count=c;
return NULL;
}

int main(){
	int N=(int)1e6;
	int count1=0,count2=0,count3=0;
	pthread_t t1,t2,t3;
	struct params p1 = {.n=N/3, .seed=1, .count=count1};
	struct params p2 = {.n=N/3, .seed=2, .count=count2};
	struct params p3 = {.n=N/3, .seed=3, .count=count3};
	pthread_create(&t1,NULL,bar,(void*)&p1);
	pthread_create(&t2,NULL,bar,(void*)&p2);
	pthread_create(&t3,NULL,bar,(void*)&p3);
	pthread_join(t1,NULL);
	pthread_join(t2,NULL);
	pthread_join(t3,NULL);
	int count=p1.count+p2.count+p3.count;
	double pi=(double)count/N*4;
	printf("Calculated value of pi = %g\n",pi);
	printf("Actual value of pi =  3.14159265\n");
return 0;
}
