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
    //ÿ���߳���һ�����һ����
    int i = MOD(threadIdx.x + blockIdx.x * blockDim.x, n);//�̸߳�����������
    int nn = n / (BLOCK_STEP * BLOCK_SIZE);//�����߳������Ӻ���߳������ѭ������
    float Fx = 0.0f;
    float Fy = 0.0f;
    float Fz = 0.0f;
    __shared__ float3 p_shared[BLOCK_SIZE];//�����߳̿����߳̿��С�Ĺ����ڴ�
    float xi = p[i].x;
    float yi = p[i].y;
    float zi = p[i].z;//�����Ӧλ�õ�����
    float dx, dy, dz, distSqr, invDist, invDist3;
    int loc;
 //ѭ��չ��
#pragma unroll 32
    for (register int j = (blockIdx.x + blockIdx.x / BLOCK_NUM); nn-- > 0; j += BLOCK_STEP) {
        //�ӵ�ǰ�鿪ʼ��ǰ�ƽ�
        j = MOD(j, BLOCK_NUM);
        //��õ�ǰ���еĵ�ǰ���λ�õ����ݣ���д��shared��
        loc = j * BLOCK_SIZE + threadIdx.x;
        p_shared[threadIdx.x] = make_float3(p[loc].x, p[loc].y, p[loc].z);
        //ͬ������ֹ�����ڴ�δ�޸���ɾ�ʹ��
        __syncthreads();
 //ѭ��չ��
#pragma unroll 32
        //�̼߳����为��������������߳̿��С�����������������
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
        //ͬ������֤�����ڴ�ʹ����Ϻ���ܸ���
        __syncthreads();
    }
    //ԭ�Ӽӷ�ֹ���ݾ������³���
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
            cudaMemcpy(buf, d_buf, bytes, cudaMemcpyDeviceToHost);//���һ��ʱ��Ƭ��д��CPU
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