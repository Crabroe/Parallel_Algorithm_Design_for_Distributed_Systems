#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"
#include "check.h"
#include <cuda_runtime.h>
#define SOFTENING 1e-9f
#define BLOCK_SIZE 128
#define BLOCK_STEP 32
#define BLOCK_NUM 32
#define MOD(a,b) ((a) - (a) / (b) * (b))
/*
 * Each body contains x, y, and z coordinate positions,
 * as well as velocities in the x, y, and z directions.
 */
typedef struct { float x, y, z, vx, vy, vz; } Body;
/*
 * Do not modify this function. A constraint of this exercise is
 * that it remain a host function.
 */

void randomizeBodies(float* data, int n) {
    for (int i = 0; i < n; i++) {
        data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    }
}
/*
 * This function calculates the gravitational impact of all bodies in the system
 * on all others, but does not update their positions.
 */
__global__ void bodyForce(Body* p, float dt, int n) {
    //每个线程算一个点的一部分
    int i = MOD(threadIdx.x + blockIdx.x * blockDim.x, n);//线程负责的天体序号
    int nn = n / (BLOCK_STEP * BLOCK_SIZE);//考虑线程数增加后的线程内外层循环次数
    float Fx = 0.0f;
    float Fy = 0.0f;
    float Fz = 0.0f;
    __shared__ float3 p_shared[BLOCK_SIZE];//声明线程块内线程块大小的共享内存
    float xi = p[i].x;
    float yi = p[i].y;
    float zi = p[i].z;//储存对应位置的数据
    float dx, dy, dz, distSqr, invDist, invDist3;
    int loc;
 //循环展开
#pragma unroll 32
    for (register int j = (blockIdx.x + blockIdx.x / BLOCK_NUM); nn-- > 0; j += BLOCK_STEP) {
        //从当前块开始向前推进
        j = MOD(j, BLOCK_NUM);
        //获得当前块中的当前相对位置的数据，并写到shared中
        loc = j * BLOCK_SIZE + threadIdx.x;
        p_shared[threadIdx.x] = make_float3(p[loc].x, p[loc].y, p[loc].z);
        //同步锁防止共享内存未修改完成就使用
        __syncthreads();
 //循环展开
#pragma unroll 32
        //线程计算其负责的天体与其余线程块大小数量的天体的作用力
        for (register int k = 0; k < BLOCK_SIZE; k++) {
            dx = p_shared[k].x - xi;
            dy = p_shared[k].y - yi;
            dz = p_shared[k].z - zi;
            distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
            invDist = rsqrtf(distSqr);
            invDist3 = invDist * invDist * invDist;
            Fx += dx * invDist3;
            Fy += dy * invDist3;
            Fz += dz * invDist3;
        }
        //同步锁保证共享内存使用完毕后才能更新
        __syncthreads();
    }
    //原子加防止数据竞争导致出错
    atomicAdd(&p[i].vx, dt * Fx);
    atomicAdd(&p[i].vy, dt * Fy);
    atomicAdd(&p[i].vz, dt * Fz);
    if (blockIdx.x / BLOCK_NUM == BLOCK_STEP - 1) {
        atomicAdd(&p[i].x, p[i].vx * dt);
        atomicAdd(&p[i].y, p[i].vy * dt);
        atomicAdd(&p[i].z, p[i].vz * dt);
    }
}
int main(const int argc, const char** argv) {
    /*
     * Do not change the value for `nBodies` here. If you would like to modify it,
     * pass values into the command line.
     */
    int nBodies = 2 << 11;
    int salt = 1;
    if (argc > 1) nBodies = 2 << atoi(argv[1]);
    /*
     * This salt is for assessment reasons. Tampering with it will result in automatic failure.
     */
    if (argc > 2) salt = atoi(argv[2]);
    const float dt = 0.01f; // time step
    const int nIters = 10;  // simulation iterations
    int bytes = nBodies * sizeof(Body);
    float* buf;
    cudaMallocHost((void**)&buf, bytes);
    /*
     * As a constraint of this exercise, `randomizeBodies` must remain a host function.
     */
    randomizeBodies(buf, 6 * nBodies); // Init pos / vel data
    float* d_buf;
    cudaMalloc((void**)&d_buf, bytes);
    Body* d_p = (Body*)d_buf;
    cudaMemcpy(d_buf, buf, bytes, cudaMemcpyHostToDevice);
    double totalTime = 0.0;
    /*
     * This simulation will run for 10 cycles of time, calculating gravitational
     * interaction amongst bodies, and adjusting their positions to reflect.
     */
     /*******************************************************************/
     // Do not modify these 2 lines of code.
    for (register int iter = 0; iter < nIters; iter++) {
        StartTimer();
        /*******************************************************************/
        /*
         * You will likely wish to refactor the work being done in `bodyForce`,
         * as well as the work to integrate the positions.
         */
        bodyForce << <BLOCK_NUM * BLOCK_STEP, BLOCK_SIZE >> > (d_p, dt, nBodies); // compute interbody forces
      /*
       * This position integration cannot occur until this round of `bodyForce` has completed.
       * Also, the next round of `bodyForce` cannot begin until the integration is complete.
       */
        if (iter == nIters - 1)
            cudaMemcpy(buf, d_buf, bytes, cudaMemcpyDeviceToHost);//最后一个时间片再写入CPU
        /*******************************************************************/
        // Do not modify the code in this section.
        const double tElapsed = GetTimer() / 1000.0;
        totalTime += tElapsed;
    }
    double avgTime = totalTime / (double)(nIters);
    float billionsOfOpsPerSecond = 1e-9 * nBodies * nBodies / avgTime;
#ifdef ASSESS
    checkPerformance(buf, billionsOfOpsPerSecond, salt);
#else
    checkAccuracy(buf, nBodies);
    printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies, billionsOfOpsPerSecond);
    salt += 1;
#endif
    /*******************************************************************/
    /*
     * Feel free to modify code below.
     */
    cudaFree(buf);
}