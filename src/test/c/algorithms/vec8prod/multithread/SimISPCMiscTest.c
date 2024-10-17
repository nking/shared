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
    // generate 128 random vector of numbers in range [1,1.65] whose product is <= 3.4E38

    // using seconds of time of day:
    unsigned int seed = time(0);
    printf("seed=%d\n", seed);
    srand(seed);

    float vin[128];
    for (int i = 0; i < 128; ++i) {
        vin[i] = 0.65f + ((float)rand() / (float)RAND_MAX);
    }

    INIT_TIME_TITLE(simsimdtestrand);

    test(vin, 128, 8);

    return NULL;
}

int main() {
   test16();
   testRand();
   return 0;
}
