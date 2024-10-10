#include <stdio.h>
#include <stdlib.h>
#include <math.h>   //  fabsf
#include <assert.h>

// Include the header file that the ispc compiler generates
#include "ispc_function.h"

int main() {
    float vin[16];

    float expAns = 1.f;

    for (int i = 0; i < 16; ++i) {
        vin[i] = (float)(i + 10);
        expAns *= vin[i];
    }

    float ans = ispc_function(16, vin);

    //printf("expAns=%f, ans=%f\n", expAns, ans);

    assert(fabsf((expAns/ans) - 1) < 5E-5);

    return 0;
}

