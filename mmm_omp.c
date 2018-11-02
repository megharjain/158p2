#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>



unsigned int matrix_checksum(int N, double *M);


void ikj_order(double *A, double *B, double *C, int N) 
{
	struct timespec *start = malloc(sizeof(struct timespec));
	struct timespec *end = malloc(sizeof(struct timespec));
	clock_gettime(CLOCK_MONOTONIC, start);
    
    
    #pragma omp parallel for schedule(dynamic)
    
    for (int i = 0; i < N; i++) 
    {

        for (int k = 0; k < N; k++) 
        {
            double r = A[i*N+k];
            
            for (int j = 0; j < N; j++)
				C[i*N+j] += r*B[k*N+j];
        }
    }
    
    
    clock_gettime(CLOCK_MONOTONIC, end);
	double elapsed = 	(end->tv_sec+end->tv_nsec/1e9)-(start->tv_sec+start->tv_nsec/1e9);
	printf("Running time: %f secs\n", elapsed);
}




int main(int argc, char **argv)
{
	double *A, *B, *C;

	if (argc != 2) 
	{
		fprintf(stderr, "Usage: %s N\n", argv[0]);
		return 1;
	}

	else if (argc == 2) 
	{
		int N = atoi(argv[1]);

		if (N > 2000 || N <= 0) 
        {
			fprintf(stderr, "%s", "Error: wrong matrix order (0 < N <= 2000)\n");         
			return 1;
		}

		A = malloc((N*sizeof(double))*(N*sizeof(double)));
		B = malloc((N*sizeof(double))*(N*sizeof(double)));
		C = malloc((N*sizeof(double))*(N*sizeof(double)));

		for (int i = 0; i < N; i++) 
        {
			for (int j = 0; j < N; j++) 
            {
				A[i*N+j] = i+j;
				B[i*N+j] = i+j*2;
			}
		}
        
        ikj_order(A, B, C, N);
    	printf("A: %u\n", matrix_checksum(N, A));
    	printf("B: %u\n", matrix_checksum(N, B));
    	printf("C: %u\n", matrix_checksum(N, C));
    
        free(A);
        free(B);
        free(C);	

		
	}

	return 0;
}
