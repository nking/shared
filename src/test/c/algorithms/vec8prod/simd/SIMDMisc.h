#include <stdio.h>
#include <math.h> //ceil()
#include <immintrin.h>  // Include AVX header
#include "../time_log.h"

//public:
extern float simd_function(int N,  float *x);

extern void cache_explore_function8();
extern void cache_explore_function4();

//private:
void *intrinsicsThread(void *arg);
