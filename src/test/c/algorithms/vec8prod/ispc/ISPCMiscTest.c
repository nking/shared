#include <stdio.h>
#include <stdlib.h>
#include <math.h>   //  fabsf
#include <assert.h>

// Include the header file that the ispc compiler generates
#include "ispc_function.h"

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

    float ansISPC = ispc_function(16, vin);

    // not avail for architecture x86_64
    //float ansISPCTask = ispc_function_tasks(16, vin);

    //printf("expAns=%f, ans=%f\n", expAns, ans);

    assert(fabsf((expAns/ansISPC) - 1) < 5E-5);
    //assert(fabsf((expAns/ansISPCTask) - 1) < 5E-5);

    printf("Done\n");

    return 0;
}

