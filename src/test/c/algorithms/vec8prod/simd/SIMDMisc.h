#include <stdio.h>
#include <math.h> //ceil()
#include <immintrin.h>  // Include AVX header
#include "../time_log.h"

//public:
extern float simd_function(int N,  float *x);

//private:
void *multInThread(void *arg);