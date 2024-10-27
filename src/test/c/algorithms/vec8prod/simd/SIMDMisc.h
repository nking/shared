#include <stdio.h>
#include <math.h> //ceil()
#include <immintrin.h>  // Include AVX header
#include "../time_log.h"

//public:
extern float simd_function(int N,  float *x);

extern void cache_explore_function8();
extern void cache_explore_function4();
extern void cache_explore_simd_more_math();
extern void cache_explore_serial_more_math(); 

//private
void *intrinsicsThread(void *arg);
void _cache_explore_serial_more_math(float *x, const int n);
void _cache_explore_simd_more_math(float *x, const int n);
