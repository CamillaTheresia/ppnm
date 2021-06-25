#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#define RND (double)(rand()%10)
#define FMT "%7.3f"


int main(int argc, char** argv){
        int n=4;
        int m=3;
if(argc>1)n=atoi(argv[1]);
if(argc>2)m=atoi(argv[2]);
        printf("n=%i\n",n);
        gsl_matrix* A=gsl_matrix_alloc(n,n);
        gsl_vector* R=gsl_vector_alloc(n);
        for(int i=0;i<m;i++)for(int j=0;j<m;j++)gsl_matrix_set(A,i,j,RND);
        gsl_linalg_QR_decomp(A,R);
        gsl_matrix_free(A);
	gsl_vector_free(R);
}
