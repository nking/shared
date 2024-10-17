#include <stdio.h>
#include <stdlib.h>
#include <math.h>   //  fabsf
#include <assert.h>
#include "../time_log.h"

// Include the header file that the ispc compiler generates
#include "ispc_function.h"

void * test(float* vin, int vLen) {

    INIT_TIME();
    START_TOT_TIME();

    float expAns = 1.f;
    for (int i = 0; i < vLen; ++i) {
        expAns *= vin[i];
    }

    STOP_TOT_TIME(serial);

    START_TOT_TIME();
    float ans = ispc_function(vLen, vin);
    STOP_TOT_TIME(simd);

    //printf("expAns=%f, ans=%f\n", expAns, ans);

    assert(fabsf((expAns/ans) - 1) < 5E-5);

    return NULL;
}

void * test16() {
    INIT_TIME_TITLE(ispctest16);

    float vin[16];
    for (int i = 0; i < 16; ++i) {
        vin[i] = (float)(i + 10);
    }

    test(vin, 16);

    return NULL;
}

void * testRand() {
    // generate 128 random vector of numbers in range [1,1.65] whose product is <= 3.4E38

    // using seconds of time of day:
    unsigned int seed = 1729184686;//time(0);
    printf("seed=%d\n", seed);
    srand(seed);

    float vin[128];
    for (int i = 0; i < 128; ++i) {
        vin[i] = 0.65f + ((float)rand() / (float)RAND_MAX);
    }

    INIT_TIME_TITLE(ispctestrand);

    test(vin, 128);

    return NULL;
}

int main() {
    test16();
    testRand();
    return 0;
}

