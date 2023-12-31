#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include<stdlib.h>
#define MIN(a,b) ((a)<(b) ? (a) : (b))

int main(int argc, char* argv[])
{
	int    count;        /* Local prime count */
	double elapsed_time; /* Parallel execution time */
	int    first;        /* Index of first multiple */
	int    global_count; /* Global prime count */
	int    high_value;   /* Highest value on this proc */
	int    i;
	int    id;           /* Process ID number */
	int    index;        /* Index of current prime */
	int    low_value;    /* Lowest value on this proc */
	char* marked;       /* Portion of 2,...,'n' */
	int    n;            /* Sieving from 2, ..., 'n' */
	int    p;            /* Number of processes */
	int    proc0_size;   /* Size of proc 0's subarray */
	int    prime;        /* Current prime */
	int    size;         /* Elements in 'marked' */
	int    odd;
	int    loc;
	char* primes;
	int    primes_size;
	int    lv;
	int    sec;
	int    chunk;

	MPI_Init(&argc, &argv);

	/* Start the timer */

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = -MPI_Wtime();

	n = atoi(argv[1]);
	odd = (n - 1) / 2 ;
	chunk = 65536;

	/* Figure out this process's share of the array, as
	  well as the integers represented by the first and
	  last array elements */

	low_value = 2 * (id * (long long)odd / p) + 3;
	high_value = 2 * ((id + 1) * (long long)odd / p - 1) + 3;
	size = ((id + 1) * (long long)odd / p - 1) - (id * (long long)odd / p) + 1;
	/* Allocate this process' share of the array */

	marked = (char*)calloc(size,sizeof(char));

	if (marked == NULL) {
		printf("Cannot allocate enough memory\n");
		MPI_Finalize();
		exit(1);
	}
	primes_size = (sqrt(n) - 1) / 2 ;
	primes = (char*)malloc(primes_size);
	if (primes == NULL) {
		printf("Cannot allocate enough memory\n");
		free(marked);
		MPI_Finalize();
		exit(1);
	}
	for (i = 0; i < primes_size; i++) primes[i] = 0;

	index = 0;
	prime = 3;
	do {
		for (i = (prime * prime - 3) / 2; i < primes_size; i += prime)
			primes[i] = 1;
		while (primes[++index]);
		prime = 2 * index + 3;
	} while (prime * prime <= sqrt(n));

	for (sec = 0; sec < size; sec += chunk) {
		index = 0;
		prime = 3;
		lv = 2 * ((low_value - 3) / 2 + sec) + 3;
		do {
			if (prime * prime > lv)
				first = (prime * prime - 3) / 2 - (lv - 3) / 2;
			else {
				if (!(lv % prime)) first = 0;
				else {
					first = prime - lv % prime;
					if (!((lv + first) % 2))
						first = (first + prime) / 2;
					else first /= 2;
				}
			}
			for (i = first + sec; i < first + sec + chunk && i < size; i += prime)
				marked[i] = 1;
			while (primes[++index]);
			prime = 2 * index + 3;
		} while (prime * prime <= n);
	}
	count = 0;
	for (i = 0; i < size; i++)
		if (!marked[i]) count++;
	MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	/* Stop the timer */

	elapsed_time += MPI_Wtime();

	/* Print the results */

	if (!id) {
		printf("There are %d primes less than or equal to %d\n",
			global_count + 1, n);
		printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
	}
	MPI_Finalize();
	return 0;
}