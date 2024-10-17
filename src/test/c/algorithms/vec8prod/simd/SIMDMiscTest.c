#include "SIMDMisc.h"
#include <assert.h>

void * test16() {

    float vin[16];
    for (int i = 0; i < 16; ++i) {
        vin[i] = (float)(i + 10);
    }

    INIT_TIME();
    START_TOT_TIME();

    float expAns = 1.f;
    for (int i = 0; i < 16; ++i) {
        expAns *= vin[i];
    }
    
    STOP_TOT_TIME(serial);

    START_TOT_TIME();
    float ans = simd_function(16, vin);
    STOP_TOT_TIME(simd);

    //printf("expAns=%f, ans=%f\n", expAns, ans);

    assert(fabsf((expAns/ans) - 1) < 5E-5);

    return NULL;
}

void * testRand() {
    return NULL;
}

int main() {
    test16();
    testRand();
    return 0;
}
