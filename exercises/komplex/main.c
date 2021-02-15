#include"komplex.h"
#include<stdio.h>
#define TINY 1e-6

int main() {
	komplex a = {1,2}, b = {3,4};
	komplex_print("a=",a);
	komplex_print("b=",b);
	printf("testing komplex_add...\n");
	komplex r = komplex_add(a,b);
	komplex R = {4,6};
	komplex_print("a+b should = ", R);
	komplex_print("a+b actually = ", r);
	printf("testing komplex_sub...\n");
	komplex d = komplex_sub(a,b);
	komplex D = {-2,-2};
	komplex_print("a-b should = ", D);
	komplex_print("a-b actually = ", d);
return 0;
}
