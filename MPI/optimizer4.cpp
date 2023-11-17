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
	register int    i;           
	int    id;           /* Process ID number */
	int    index;        /* Index of current prime */
	int    low_value;    /* Lowest value on this proc */
	char* marked;       /* Portion of 2,...,'n' */
	int    n;            /* Sieving from 2, ..., 'n' */
	int    p;            /* Number of processes */
	int    proc0_size;   /* Size of proc 0's subarray */
	register int    prime;        /* Current prime */
	int    size;         /* Elements in 'marked' */
	int    odd;          /* 奇数总数 */
	char* primes;        /* 进程内素数表 */
	int    primes_size;  /* 进程内内素数表大小 */
	int    lv;           /* 进程最小数据对应的表格位置 */
	int    sec;          /* 各数据块第一个数据对应的表格位置 */
	int    chunk;        /* Cache块大小 */
	//MPI框架初始化
	MPI_Init(&argc, &argv);

	/* Start the timer */

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Barrier(MPI_COMM_WORLD);
	//记录代码开始时间
	elapsed_time = -MPI_Wtime();
	//读入数据量并计算对应奇数数量
	n = atoi(argv[1]);
	odd = (n - 1) / 2;
	//读入Cache大小
	chunk = 65536;
	/* Figure out this process's share of the array, as
	  well as the integers represented by the first and
	  last array elements */
	//计算各进程所分配的数上下限及范围
	low_value = 2 * (id * (long long)odd / p) + 3;
	high_value = 2 * ((id + 1) * (long long)odd / p - 1) + 3;
	size = ((id + 1) * (long long)odd / p - 1) - (id * (long long)odd / p) + 1;

	/* Allocate this process' share of the array */
	//创建表格并判断内存是否足够
	marked = (char*)calloc(size,sizeof(char));
	if (marked == NULL) {
		printf("Cannot allocate enough memory\n");
		MPI_Finalize();
		exit(1);
	}
	//创建进程内素数表并判断内存是否足够
	primes_size = (sqrt(n) - 1) / 2;
	primes = (char*)malloc(primes_size);
	if (primes == NULL) {
		printf("Cannot allocate enough memory\n");
		free(marked);
		MPI_Finalize();
		exit(1);
	}
	//初始化素数表为0
	for (i = 0; i < primes_size; i++) primes[i] = 0;
	//重构循环下的进程内素数查找
	for (sec = 0; sec < size; sec += chunk) {//逐个查找放入Cache的数据块
		index = 0;
		prime = 3;
		do {
			//判断第一个素数倍数
			if (prime * prime > sec)
				first = (prime * prime - 3) / 2-sec;
			else {
				if (!(sec % prime)) first = 0;
				else {
					first = prime - sec % prime;
					if (!((sec + first) % 2))
						first = (first + prime) / 2;
					else first /= 2;
				}
			}
			//将当前素数所有倍数的表格位置置1
			for (i = (prime * prime - 3) / 2; i < primes_size; i += prime)
				primes[i] = 1;
			//查找下一个素数
			while (primes[++index]);
			prime = 2 * index + 3;
		} while (prime * prime <= sqrt(n));
	}
	//重构循环下的全局素数查找
	for (sec = 0; sec < size; sec += chunk) {//外循环为逐个查找Cache中的数据块
		index = 0;
		prime = 3;
		lv = 2 * ((low_value - 3) / 2 + sec) + 3;
		do {
			//各进程判断分配数据中第一个素数倍数
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
			//各进程将分配数据中当前素数所有倍数在表格中置1
			for (i = first + sec; i < first + sec + chunk && i < size; i += prime)
				marked[i] = 1;
			//查找下一个素数
			while (primes[++index]);
			prime = 2 * index + 3;
		} while (prime * prime <= n);
	}
	//统计分配的数中在表格对应位置为0的数
	count = 0;
	for (i = 0; i < size; i++)
		if (!marked[i]) count++;
	//规约所有进程结果
	MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	/* Stop the timer */
	//计算花费时间
	elapsed_time += MPI_Wtime();
	/* Print the results */
	//0号进程返回素数总数
	if (!id) {
		printf("There are %d primes less than or equal to %d\n",
			global_count+1 , n);
		printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
	}
	MPI_Finalize();
	return 0;
}