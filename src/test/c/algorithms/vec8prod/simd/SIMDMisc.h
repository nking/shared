#include <stdio.h>
#include <math.h> //ceil()
#include <pthread.h>
#include <immintrin.h>  // Include AVX header

//public:
extern float simd_function(int N,  float *x);

//private:
void *multInThread(void *arg);