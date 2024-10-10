#include <stdio.h>
#include <math.h>   //  fabsf
#include <assert.h>
#include <immintrin.h>  // Include AVX header

#include "SIMDMisc.h"

int main() {

    // TODO: consider making a test with very large number of numbers and compare the
    // runtime executions
    float vin[16];
    for (int i = 0; i < 16; ++i) {
        vin[i] = (float)(i + 10);
    }

    float expAns = 1.f;
    for (int i = 0; i < 16; ++i) {
        expAns *= vin[i];
    }

    float ans = simd_function(16, vin);

    //printf("expAns=%f, ans=%f\n", expAns, ans);

    assert(fabsf((expAns/ans) - 1) < 5E-5);

    printf("Done\n");

    return 0;
}
