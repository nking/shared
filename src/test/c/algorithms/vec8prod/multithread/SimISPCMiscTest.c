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

void test0() {

   int n = 16;
   int isWidth = 8; // instruction set width
   float x[n]; // loading the data "coherently"
   float origX[n]; // loading the data "coherently"
   float expAns = 1.f;
   // simply using 10:25 for values
   for (int i = 0; i < n; ++i) {
      x[i] = (float)(i + 10);
      origX[i] = x[i];
      expAns *= x[i];
   }

   float ans = multithread_function(x, n, isWidth);
   assert(fabsf((expAns/ans)-1) < 5E-5);

   //permute x to show order doesn't matter
   fisherYatesShuffle(origX, n);

   ans = multithread_function(origX, n, isWidth);

   assert(fabsf((expAns/ans)-1) < 5E-5);

   printf("Done\n");
}

int main() {
   test0();
   return 0;
}
