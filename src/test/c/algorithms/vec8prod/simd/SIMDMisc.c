#include "SIMDMisc.h"

//TODO: add timing with RDTSC clock

struct thread_data {
   float* x;
   int n;
   int instanceNumber;
   int isWidth;
};

float simd_function(int N,  float * x) {
   // use 8-wide vec avx256
   int isWidth = 8;

   int nThreads = (int)ceil( (float) N / isWidth);

   for (int i = 0; i < nThreads; ++i) {

      struct thread_data data;
      data.x = x;
      data.instanceNumber = i;
      data.isWidth = isWidth;
      data.n = N;

      pthread_t thread_id;

      pthread_create(&thread_id, NULL, multInThread, &data);
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

void *multInThread(void *arg) {

   struct thread_data *data = (struct thread_data *)arg;

   int idx0 = data->isWidth * data->instanceNumber;
   //int idx1 = data->isWidth * (data->instanceNumber + 1) -1;

   //printf("thread %d, idx0=%d, x[idx0]=%f  (xref=%p)\n", data->instanceNumber,
   //    idx0, data->x[idx0], (void*)&(data->x));

   //ps - vectors contain floats (ps stands for packed single-precision)

   //256-bit vector containing 8 floats
   __m256 avx_x = _mm256_loadu_ps(&(data->x[idx0]));

   int nIter = 0;
   while (nIter < 3) {

      int shift = 1 << nIter;

      __m256 avx_y = _mm256_loadu_ps(&( ((float*)&avx_x)[shift] ));

      //printf("   thread %d, nIter=%d, shift=%d, avx_x, avx_y:\n", data->instanceNumber, nIter, shift);
      //printVec(avx_x);
      //printVec(avx_y);

      avx_x = _mm256_mul_ps(avx_x, avx_y);

      ++nIter;
   }

   // store result back into x[idx0]
   data->x[idx0] = ((float*)&avx_x)[0];//_mm_cvtss_f32(_mm256_extractf128_ps(avx_x, 0));

   //printf("   thread %d, result x[idx0]=%f\n", data->instanceNumber, data->x[idx0]);

   return NULL;
}


