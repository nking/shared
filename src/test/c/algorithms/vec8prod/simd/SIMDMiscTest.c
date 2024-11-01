#include "SIMDMisc.h"
#include <assert.h>
#include <stdlib.h>

/* 
context switching is expensive here and in multithread without SIMD.

The algorithmic intensity for this SIMD method is
   time for 3 calcs / time for (4 efficient loads from memory + 3 gather loads from memory)

The better algorithmic intensity for SIMD method compared to multithread is not
easily seen among the much larger times needed for the context switching.
 */
void * test(float* vin, int vLen) {
    
    INIT_TIME();
    START_TOT_TIME();

    float expAns = 1.f;
    for (int i = 0; i < vLen; ++i) {
        expAns *= vin[i];
    }
    
    STOP_TOT_TIME(serial);

    START_TOT_TIME();
    float ans = simd_function(vLen, vin);
    STOP_TOT_TIME(simd);

    //printf("expAns=%f, ans=%f\n", expAns, ans);

    assert(fabsf((expAns/ans) - 1) < 5E-5);

    return NULL;
}

void * test16() {
    INIT_TIME_TITLE(simdtest16);

    float vin[16];
    for (int i = 0; i < 16; ++i) {
        vin[i] = (float)(i + 10);
    }

    test(vin, 16);

    return NULL;
}

void * testRand() {
    // generate  random vector of numbers whose product is <= 3.4E38
    int n = 1<<16;

    // using seconds of time of day:
    unsigned int seed = time(0);
    printf("seed=%d\n", seed);
    srand(seed);

    float factor = 1.0f  - pow(MAXFLOAT, (1.f/(float)n));

    float vin[n];
    for (int i = 0; i < n; ++i) {
        vin[i] = 1.0f + factor * ((float)rand() / (float)RAND_MAX);
    }

    INIT_TIME_TITLE(simdtestrand);

    test(vin, n);

    return NULL;
}

int compare(const void *a, const void *b) {
    return (*(double *)a - *(double *)b);
}

void * estimateCPUStats() {
    int nTests = 10;
    int n = 1000000;
    clock_t start;
    struct timespec tw1, tw2;
    double times[nTests];
    double timesW[nTests];
    double timesWNS[nTests];
    for (int i = 0; i < nTests; ++i) {
	clock_gettime(CLOCK_MONOTONIC, &tw1);
        start = clock();
        for (int j = 0; j < n; ++j) {
            sqrtf((float)j);
        }
        times[i] = (double)(clock() - start);
	clock_gettime(CLOCK_MONOTONIC, &tw2);
	timesW[i] = (1.0 * tw2.tv_sec + 1E-9 * tw2.tv_nsec)
		  - (1.0 * tw1.tv_sec + 1E-9 * tw1.tv_nsec);
	timesWNS[i] = tw2.tv_nsec - tw1.tv_nsec;
    }

    qsort(times, nTests, sizeof(double), compare);
    qsort(timesW, nTests, sizeof(double), compare);
    qsort(timesWNS, nTests, sizeof(double), compare);
    //for (int i = 0; i < nTests; ++i) {
    //    printf("timeW=%f\n", timesW[i]);
    //}

    // number of clocks for n instructions
    double med = times[nTests/2];
    double instrPerClock = (double)n/med;

    printf("number of clocks measured for n=%d instructions=%f\n", n, med);
    printf("1 clock = %f instructions\n", instrPerClock);
    printf("CLOCKS_PER_SEC = %d\n", CLOCKS_PER_SEC);

    // using wall clock time, have instructions time in ns
    double medW = (timesW[nTests/2]); 
    printf("Wall: time for n=%d instructions in sec = %f\n", n, medW);
    printf("Wall: time for n=%d instructions in nsec = %f\n", n, (medW/1E-9));
    printf("Wall: time for n=%d instructions in nsec = %f\n", n, timesWNS[nTests/2]);
     
    double estClockPerSec = ((double)n/medW)*(1./instrPerClock); 
    printf("estimated CLOCKS_PER_SEC = %f\n", estClockPerSec);
      
    return NULL;
}

int main() {
    // uncomment as needed.  no coding at this time for runtime flags to choose tests.
    //test16();
    testRand();
    //estimateCPUStats();
    //cache_explore_function8();
    //cache_explore_function4();
    //cache_explore_simd_more_math();
    //cache_explore_serial_more_math();
    return 0;
}
