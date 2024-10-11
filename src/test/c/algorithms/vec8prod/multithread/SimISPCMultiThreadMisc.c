#include "SimISPCMultiThreadMisc.h"

/*
emmulating the ISPC example from stanford lecture 4 CS149 parallel computing
with software threads.
doing so just to look at the divide and conquer math and vector elements.
https://gfxcourses.stanford.edu/cs149/fall23/lecture/

   the ISPC code is:
   export void vec8product(uniform float* x, uniform float* result) {
      float val1 = x[programIndex];
      float val2 = shift(val1, -1);
      if (programIndex % 2 == 0)
         val1 = val1 * val2;

      val2 = shift(val1, -2);
      if (programIndex % 4 == 0)
         val1 = val1 * val2;

      val2 = shift(val1, -4);
      if (programIndex % 8 == 0)
         *result = val1 * val2;
   }

   the ISPC example takes a large vector of data, divides it among
   N "gangs" (name is from "gang scheduling")
   and computes the product of each gang's data separately,
   then combines the result.

   The ISPC gang uses single process multiple data (SPMD) model.

   below are notes from Pharr and Mark, "ispc: A SPMD Compiler for High Performance CPU Programming"
*/


struct thread_data {
   float* x; 
   int instanceNumber;
   int isWidth;
};


void *multInThread(void *arg) {

   // x is a uniform shared variable array.

   struct thread_data *data = (struct thread_data *)arg;

   // begin spmd replacement.  
   // emulating the vector lane as offsets in the uniform shared variable x
   int idx0 = data->isWidth * data->instanceNumber;
   int idx1 = data->isWidth * (data->instanceNumber + 1) -1;

   // the span of this program is log_2(isWidth) = 3 for SPMD execution

   for (int prId = 2; prId <= data->isWidth; prId <<= 1) {
      // each prId is 1 round of computation on the vector of data between 
      // idx0 and idx1 .
      //  we pretend the vector operations in each prId level is SPMD 
      // (single process multiple data) and is executed
      // simultaneously for this lane.

      int off0 = prId - 1;
      int off1 = (prId/2) - 1;

      // pretend this is SPMD
      for (int j = idx0; j < idx1; j+= prId) {
         data->x[j + off0] *= data->x[j+off1];
      }
   }
   return NULL;
}

float multithread_function(float* x, int xLen, int isWidth) {
   assert(xLen % isWidth == 0); 

   int nInstances = xLen/isWidth;

   /* emulating ISPC dividing the work among xLen/isWidth workers.
      NOTE that the ISPC gang instances all run in the same hardware 
   thread and context.
   there are no more than 2 program instances running in a gang - no more 
   than twice the SIMD width of the hardware its running on.
      e.g. for a CPU running 4-wide SSE instruction set, there can be 4 or 8 
          program instances running in a gang.
      for a CPU running 8-wide AVX, there can be 8 or 16 program instances 
          running in a gang.
   The SPMD parallelization is across SIMD lanes of a single core.

   Note: a gang is roughly equiv to a CUDA warp.
   */


   for (int i = 0; i < nInstances; ++i) {

      struct thread_data data;
      data.x = x;
      data.instanceNumber = i;
      data.isWidth = isWidth;

      pthread_t thread_id;

      pthread_create(&thread_id, NULL, multInThread, &data);
      pthread_join(thread_id, NULL);
   }

   // when done with all threads, multiply results for each lane
   float res = 1.f;
   for (int i = isWidth-1; i < xLen; i += isWidth) {
      res *= x[i];
   }

   return res;
}