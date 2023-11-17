#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

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

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();//记录开始时间
    //检查参数是否足够
    if (argc != 2) {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }
    //读取筛选最大数
    n = atoi(argv[1]);
    //为每个进程分配筛选数范围
    low_value = 2 + id * ((long long)n - 1) / p;
    high_value = 1 + (id + 1) * ((long long)n - 1) / p;
    size = high_value - low_value + 1;
    /* Bail out if all the primes used for sieving are
       not all held by process 0 */
    //计算0号进程分配数
    proc0_size = (n - 1) / p;
    
    //判断0号进程是否分配到所有素数
    if ((long long)((2 + proc0_size)) * (long long)((2 + proc0_size)) < n) {
        if (!id) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    /* Allocate this process's share of the array. */
    //创建表并判断是否有足够内存
    marked = (char*)malloc(size);

    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }
    //初始化表格为0
    for (i = 0; i < size; i++) marked[i] = 0;
    //0号进程从表格第一个数开始选素数
    if (!id) index = 0;
    //第一个素数为2
    prime = 2;
    //各进程计算分配的数中的素数第一个倍数
    do {
        if (prime * prime > low_value)
            first = prime * prime - low_value;
        else {
            if (!(low_value % prime)) first = 0;
            else first = prime - (low_value % prime);
        }
        for (i = first; i < size; i += prime) marked[i] = 1;
        if (!id) {
            while (marked[++index]);
            prime = index + 2;
        }
        if (p > 1) MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } while (prime * prime <= n);
    //统计所有表格中不为0的数
    count = 0;
    for (i = 0; i < size; i++)
        if (!marked[i]) count++;
    if (p > 1) MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
        0, MPI_COMM_WORLD);
     else {
         global_count = count;
     }

     /* Stop the timer */

    elapsed_time += MPI_Wtime();


    /* Print the results */
    //0号进程返回素数数量结果
    if (!id) {
        printf("There are %d primes less than or equal to %d\n",
            global_count, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}