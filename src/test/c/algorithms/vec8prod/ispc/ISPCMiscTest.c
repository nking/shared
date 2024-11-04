#include <stdio.h>
#include <stdlib.h>
#include <math.h>   //  fabsf
#include <assert.h>
#include "../time_log.h"

// Include the header file that the ispc compiler generates
#include "ispc_function.h"

void * test(float* vin, int vLen, float tol) {

    INIT_TIME();
    START_TOT_TIME();

    float expAns = 1.f;
    for (int i = 0; i < vLen; ++i) {
        expAns *= vin[i];
    }
    //printf("expAns=%f\n", expAns);

    STOP_TOT_TIME(serial);

    int nTests = 10;
    for (int i = 0; i < nTests; ++i) {
        START_TOT_TIME();
        //float ans = ispc_function_tasks(vLen, vin);//can have larger running time than ispc_function
        float ans = ispc_function(vLen, vin);
        //float ans = ispc_function2(vLen, vin);
        STOP_TOT_TIME(ispc);
        //printf("expAns=%f, ans=%f\n", expAns, ans);
        assert(fabsf((expAns/ans) - 1) < tol);
    }

    return NULL;
}

void * test16() {
    INIT_TIME_TITLE(ispctest16);

    float vin[16];
    for (int i = 0; i < 16; ++i) {
        vin[i] = (float)(i + 10);
    }

    test(vin, 16, 0.00005f);

    return NULL;
}

void * testRandLarge() {
    int n = 1<<16;
    // generate random vector of numbers whose product is <= 3.4E38

    // using seconds of time of day:
    unsigned int seed = time(0);
    //seed = 1730760287;
    printf("seed=%d\n", seed);
    srand(seed);

    float factor = powf(MAXFLOAT, (1.f/(float)n)) - 1.0f;

    float vin[n];
    for (int i = 0; i < n; ++i) {
        vin[i] = 1.0f + factor * ((float)rand() / (float)RAND_MAX);
    }
    //printf("\n");

    INIT_TIME_TITLE(ispctestrandLarge);

    test(vin, n, 0.0001f);

    return NULL;
}

void * _test(float* vin, int vLen) {
    ispc_higher_arith_func(vLen, vin);
    return NULL;
}
void * testHigherArithInt() {
    int n = 262144;
    // generate random vector of numbers whose product is <= 3.4E38

    // using seconds of time of day:
    unsigned int seed = time(0);
    printf("seed=%d\n", seed);
    srand(seed);

    float factor = 1.0f  - powf(MAXFLOAT, (1.f/(float)n));
    float vin[n];
    for (int i = 0; i < n; ++i) {
        vin[i] = 1.0f + factor * ((float)rand() / (float)RAND_MAX);
    }

    int nTests = 10;

    INIT_TIME_TITLE(ispc262144HigherArithInt);
    INIT_TIME();

    for (int i = 0; i < nTests; ++i) {
        START_TOT_TIME();
        _test(vin, n);
        STOP_TOT_TIME(ispc);
    }    

    return NULL;
}


int main() {
    //test16();
    testRandLarge();
    //testHigherArithInt();
    return 0;
}

