#include <stdio.h>
#include <stdlib.h>

// Include the header file that the ispc compiler generates
#include "ispc_function.h"

int main() {
    float vin[16], vout[16];

    float expAns = 1.f;

    // Initialize input buffer
    for (int i = 0; i < 16; ++i) {
        vin[i] = (float)(i + 10);
        expAns *= vin[i];
    }

    // Call simple() function from simple.ispc file
    ispc_function(16, vin, vout);

    // Print results
    for (int i = 0; i < 16; ++i)
        printf("%d: simple(%f) = %f\n", i, vin[i], vout[i]);

    printf("expAns=%f\n", expAns);

    return 0;
}

