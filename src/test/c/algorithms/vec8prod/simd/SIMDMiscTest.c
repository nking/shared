#include "SIMDMisc.h"
#include <assert.h>

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

int main() {
    test16();
    testRand();
    return 0;
}
