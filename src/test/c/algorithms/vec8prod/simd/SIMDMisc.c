#include "SIMDMisc.h"
#include <pthread.h>

struct thread_data {
   float* x;
   int n;
   int instanceNumber;
   int isWidth;
};

/*timing tests to explore when cache is filled
 cache_explore_function8()
 cache_explore_function4()

    level 3 Cache: 3 MB shared (not used here)
    Cache Level 1:  32k for each core
    Cache Level 2:  256k for each core
    max memory bus: 26 GB/sec
    processor : 1.6 GHz

    at each cycle memory transfer at most:
        max memory bus (26 GB/sec) * (sec/1.6G cycles) = 16 Bytes / cycle
          = 128 bits / cycle

    this algorithm on this computer is memory bandwidth bound.

    looking at the peaks in data load times:  cache evictions?

    8-wide vec of 4 bytes = 32 bytes = 256 bits

    measured avg load is 1 cycle

    measured diff in i = ~1500 to ~7000 for loads that are > 10 stdev above avg.
    (32 * 3000)/(1024) ~ 100 kB. 

    can use tools like Intel® VTune™ Profiler to explore further.
*/
void cache_explore_function8() {
    INIT_TIME();
    START_THR_TIME();

    float x[8]; 
    for (int i = 0; i < 125000; ++i) {
       for (int j = 0; j < 8; ++j) {
           x[j] = i + j;
       }
       START_D_TIME();
       //__m256 avx_x = _mm256_loadu_ps(&x[0]);
       __m256 avx_x = _mm256_loadu_ps(&x[0]);
       STOP_D_TIME2(i);
    }
    STOP_THR_TIME(loopnotthr);
}

void cache_explore_function4() {
    INIT_TIME();
    START_THR_TIME();

    float x[4]; 
    for (int i = 0; i < 250000; ++i) {
       for (int j = 0; j < 4; ++j) {
           x[j] = i + j;
       }
       START_D_TIME();
       //__m128 avx_x = _mm_loadu_ps(&x[0]);
       __m128 avx_x = _mm_load_ps(&x[0]);
       STOP_D_TIME2(i);
    }
    STOP_THR_TIME(loopnotthr);
}

float simd_function(int N,  float * x) {
   // use 8-wide vec avx256
   int isWidth = 8;

   //NOTE: no thread pooling is used
   int nThreads = (int)ceil( (float) N / isWidth);

   for (int i = 0; i < nThreads; ++i) {

      struct thread_data data;
      data.x = x;
      data.instanceNumber = i;
      data.isWidth = isWidth;
      data.n = N;

      pthread_t thread_id;

      pthread_create(&thread_id, NULL, intrinsicsThread, &data);
      pthread_join(thread_id, NULL);
   }

   float res = 1.f;
   for (int i = 0; i < N; i += isWidth) {
      res *= x[i];
   }

   return res;
}

void printVec(__m256 vec) {
   float* f = (float*)&vec;
    //__m128 _mm256_extractf32x4_ps (__m256 a, int imm8)
    //__m128 _mm256_extractf128_ps (__m256 a, const int imm8)
    //float _mm_cvtss_f32 (__m128 a)
    //float tx = _mm_cvtss_f32(_mm256_extractf128_ps(avx_x, 0));
    printf("      ");
    for (int i = 0; i < 8; ++i) {
        printf("%.0f, ", f[i]);
    }
    printf("\n");
}

void *intrinsicsThread(void *arg) {
   INIT_TIME();
   START_THR_TIME();

   struct thread_data *data = (struct thread_data *)arg;

   int idx0 = data->isWidth * data->instanceNumber;
   //int idx1 = data->isWidth * (data->instanceNumber + 1) -1;

   //printf("thread %d, idx0=%d, x[idx0]=%f  (xref=%p)\n", data->instanceNumber,
   //    idx0, data->x[idx0], (void*)&(data->x));

   //ps - vectors contain floats (ps stands for packed single-precision)

   //256-bit vector containing 8 floats
   START_D_TIME();
   __m256 avx_x = _mm256_loadu_ps(&(data->x[idx0]));
   STOP_D_TIME(load);

   // shift right by 1 (intel little endian)
   START_D_TIME();
   //__m256 avx_y = _mm256_permutevar8x32_ps(avx_x, _mm256_set_epi32(0,7,6,5,4,3,2,1));
   //float * ptr = ((float*)&avx_x);
   //__m256 avx_y = _mm256_loadu_ps(&(ptr[1]));
   __m256 avx_y = _mm256_loadu_ps(&(((float*)&avx_x)[1]));
   STOP_D_TIME(loadpermute);
   
   avx_x = _mm256_mul_ps(avx_x, avx_y);

   // shift right by 2 (intel little endian)
   START_D_TIME();
   //avx_y = _mm256_permutevar8x32_ps(avx_x, _mm256_set_epi32(0,0, 7,6,5,4,3,2));
   avx_y = _mm256_loadu_ps(&(((float*)&avx_x)[2]));
   STOP_D_TIME(loadpermute);

   avx_x = _mm256_mul_ps(avx_x, avx_y);

   // shift right by 4 (intel little endian)
   START_D_TIME();
   //avx_y = _mm256_permutevar8x32_ps(avx_x, _mm256_set_epi32(0,0,0,0, 7,6,5,4));
   avx_y = _mm256_loadu_ps(&(((float*)&avx_x)[4]));
   STOP_D_TIME(loadpermute);

   avx_x = _mm256_mul_ps(avx_x, avx_y);

   // store result back into x[idx0]
   START_D_TIME();
   data->x[idx0] = ((float*)&avx_x)[0];//_mm_cvtss_f32(_mm256_extractf128_ps(avx_x, 0));
   STOP_D_TIME(store);

   //printf("   thread %d, result x[idx0]=%f\n", data->instanceNumber, data->x[idx0]);
   
   STOP_THR_TIME(thr);

   return NULL;
}

float * generate_float_array(int n, float min, float max) {
   unsigned int seed = time(0);
   seed = 1234567;
   printf("seed=%d\n", seed);
   srand(seed);

   float x[n];
   
   float factor = max - min;
   for (int i = 0; i < n; ++i) {
      x[i] = min + factor * ((float)rand() / (float)RAND_MAX);
   }

   float *ptr;
   ptr = &x[0];
   return ptr;
}

void cache_explore_simd_more_math() {
   int n = 250000;
   float * x = generate_float_array(n, 1.0f, 10.f);
   _cache_explore_simd_more_math(x, n);
}
void cache_explore_serial_more_math() {
   int n = 250000;
   float * x = generate_float_array(n, 1.0f, 10.f);
   _cache_explore_serial_more_math(x, n);
}

void _cache_explore_simd_more_math(float *x, const int n) {
   // goal is to compare serial function to serial simd-ized function
   // for an algorithm that simply performs alot of math on data,
   // independently.
   
   // intrinsics multiplication is by vectors.
   // e.g. unary operations: reciprocal, sqrt
   
   // will use 4-wide vectors of 32-bit floats to stay within this platforms' memory bandwidth.

   float factors[4] = {3.f, 3.f, 3.f, 3.f};
   __m128 avx_f = _mm_loadu_ps(&factors[0]);
   float s[4] = {1.5f, 1.5f, 1.5f, 1.5f};
   __m128 avx_s = _mm_loadu_ps(&s[0]);

   INIT_TIME();
   START_THR_TIME();
   for (int i = 0; i < n; i += 4) {
      START_D_TIME();
      __m128 avx_x = _mm_load_ps(&x[i]);
      STOP_D_TIME2(i);

      for (int j = 0; j < 4; ++j) {
         avx_x = _mm_fmaddsub_ps(avx_x, avx_f, avx_s);
         avx_x = _mm_rcp_ps(avx_x);
         avx_x = _mm_sqrt_ps(avx_x);
      }

      START_D_TIME();
      _mm_store_ps(&x[i], avx_x);
      STOP_D_TIME(store);
   }
   STOP_THR_TIME(simdmoremath);
}

void _cache_explore_serial_more_math(float *x, const int n) {
   
   float factors[4] = {3.f, 3.f, 3.f, 3.f};
   float s[4] = {+1.5f, -1.5f, +1.5f, -1.5f};

   INIT_TIME();
   START_THR_TIME();
   for (int i = 0; i < n; i += 1) {
      for (int j = 0; j < 4; ++j) {
         x[i] = (x[i] * factors[0]) + s[i%2];
         x[i] = 1./x[i];
         x[i] = sqrtf(x[i]);
      }
   }
   STOP_THR_TIME(serialmoremath);
}

void cache_explore_serial_sterling_gamma() {

   unsigned int seed = time(0);
   //seed = 1234567;
   printf("seed=%d\n", seed);
   srand(seed);

   // notes here are included to explore how to use intrinsics for 
   // sterling gamma function

   const float SQRT2PI    = 2.50662827463100024157f;
   const float SC[] = {1.f, 0.08333333333333333f, 0.003472222222222222f
       -0.0026813271604938273f, -0.00022947209362139917f};
   float tmp[5] = {1.f, 0.f, 0.f, 0.f};

   int n = 1 << 16;
   float x[n];
   float factor = 20.f;

   for (int i = 0; i < n; ++i) {
      x[i] = 0.00001f + factor * ((float)rand() / (float)RAND_MAX);
   }
   
   for (int i = 0; i < n; ++i) {

        // to SIMD-ize, need to use bit masks for branchless programming
        if (x[i] < 0) {
           x[i] = FP_NAN;
           return;
        } else if (x[i] == 0) {
           x[i] = 1;
           return;
        }

        // simd rcp_ps
        tmp[1] = 1./x[i]; // r1
        tmp[2] = tmp[0]*tmp[0]; //r2
        tmp[4] = tmp[1]*tmp[1]; //r4
        tmp[3] = tmp[1] * tmp[4];
 
        // simd dp_ps dot product
        float f2 = 0;
        for (int j = 0; j < 4; ++j) {
             f2 += tmp[j] * SC[j];
        }

        //SIMD sqrt
        float f1 = sqrtf(x[i]);

        //SIMD pow_ps
        float f3 = powf(x[i]/M_E, x[i]);
        
        //SIMD _mul_ps, with these factors placed in a vector... a gather not aligned load unfortunately.
        // and making the factors requires scatter and gather ops
        x[i] = SQRT2PI * f1 * f2 * f3;
    }
}        
