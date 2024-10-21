#include "SimISPCMultiThreadMisc.h"

void fisherYatesShuffle(float*x, int xLen) {
    float tmp;
    for (int i = xLen - 1; i >= 0; --i) {
        // random number between 0 and i, inclusive
        int j = rand() % (i + 1);
        if (i != j) {
            tmp = x[i];
            x[i] = x[j];
            x[j] = tmp;
        } 
    }
}

/* 
context switching is expensive here and in SIMD (no thread pooling in either).

The algorithmic intensity for the multInThread method is
   time for calcs / time for data ops

the number of calculations for the partition of 8 
   N/2 calcs, then N/4, then N/8, ...
   which is a geometric series of (N/2) * series with r = 0.5
   so the number of calcs = (N/2) * (1/(1-0.5) = 8 = O(N)

there are 2 data ops for each calculation.
algorithmic intensity for the multInThread method is
    1/2 
*/

void * test(float* vin, int vLen, int isWidth) {

    INIT_TIME();
    START_TOT_TIME();

    float expAns = 1.f;
    for (int i = 0; i < vLen; ++i) {
        expAns *= vin[i];
    }

    STOP_TOT_TIME(serial);

    START_TOT_TIME();
    float ans = multithread_function(vin, vLen, isWidth);
    STOP_TOT_TIME(simsimd);

    //printf("expAns=%f, ans=%f\n", expAns, ans);

    assert(fabsf((expAns/ans) - 1) < 5E-5);

    return NULL;
}

void * test16() {
    INIT_TIME_TITLE(simsimdtest16);

    float vin[16];
    for (int i = 0; i < 16; ++i) {
        vin[i] = (float)(i + 10);
    }

    test(vin, 16, 8);

    return NULL;
}

void * testRand() {
    // generate random vector of numbers whose product is <= 3.4E38

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

    INIT_TIME_TITLE(simsimdtestrand);

    test(vin, n, 8);

    return NULL;
}

int main() {
   test16();
   testRand();
   return 0;
}
