#include "time_log_ispc.h"

/*
Implicit SPMD Program Compiler (ispc)

https://ispc.github.io/ispc.html#uniform-and-varying-qualifiers
If a variable has a uniform qualifier, then there is only a single instance of that variable shared by all program
instances in a gang. (In other words, it necessarily has the same value across all of the program instances.)
In addition to requiring less storage than varying values, uniform variables lead to a number of performance advantages
when they are applicable (see Uniform Control Flow, for example.) Varying variables may be qualified with varying,
though doing so has no effect, as varying is the default.

ISPC uses an SPMD model.
https://en.wikipedia.org/wiki/Single_program,_multiple_data

to compile:
   ispc ispc_function.ispc -o ispc_function.o
   see CMakeLists.txt
*/

void print_float_array(uniform float * uniform array, uniform int i0, uniform int i1) {
    for (uniform int i = i0; i < i1; ++i) {
        #pragma ignore warning(perf)
        print("%,\n",  array[i]);
    }
    //print("\n");
}

void vecprod_1(uniform float * uniform x, uniform int ishift) {
    INIT_TIME();
    START_THR_TIME();

    START_D_TIME();
    float val1 = x[programIndex];
    START_D_TIME();

    START_D_TIME();
    float val2 = shift(val1, ishift);
    START_D_TIME();

    //print("val1=        %\n", val1);
    //print("shift=% val2=%\n", ishift, val2);
    x[programIndex] = val1 * val2;

    STOP_THR_TIME();
}

void higher_arith_intensity(uniform float * uniform x) {
    
    INIT_TIME();
    START_THR_TIME();

    float val1 = x[programIndex];

    for (int i = 0; i < 4; ++i) {
        val1 *= 4;
        if (programIndex % 2 == 0) {
            val1 += 1.5;
        } else {
            val1 -= 1.5;
        }
        val1 = 1./val1;
        val1 = sqrt(val1);
    }

    x[programIndex] = val1;

    STOP_THR_TIME();
}
void vecprod(uniform float * uniform x, uniform float * uniform y) {
    //print("\nprogramIndex=%\n", programIndex);

    INIT_TIME();
    START_THR_TIME();

    START_D_TIME();
    float val1 = x[programIndex];
    STOP_D_LOAD_TIME();
    START_D_TIME();
    float val2 = shift(val1, 1);
    STOP_D_LOAD_TIME();

    //print("val1=%\n", val1);
    //print("val2=%\n", val2);
    if (programIndex % 2 == 0) {
        val1 = val1 * val2;
    }

    START_D_TIME();
    val2 = shift(val1, 2);
    STOP_D_LOAD_TIME();
    if (programIndex % 4 == 0) {
        val1 = val1 * val2;
    }

    START_D_TIME();
    val2 = shift(val1, 4);
    STOP_D_LOAD_TIME();
    if (programIndex % 8 == 0) {
        val1 = val1 * val2;
        //*result = val1;
        if (programCount == 8) {
            START_D_TIME();
            y[programIndex] = val1;
            STOP_D_STORE_TIME();
        }
    }

    if (programCount > 8) {
        START_D_TIME();
        val2 = shift(val1, 8);
        STOP_D_LOAD_TIME();
        if (programIndex % 16 == 0) {
            val1 = val1 * val2;
            if (programCount == 16) {
                START_D_TIME();
                y[programIndex] = val1;
                STOP_D_STORE_TIME();
            }
        }
    }

    if (programCount > 16) {
        START_D_TIME();
        val2 = shift(val1, 16);
        STOP_D_LOAD_TIME();
        if (programIndex % 32 == 0) {
            val1 = val1 * val2;
            if (programCount == 32) {
                START_D_TIME();
                y[programIndex] = val1;
                STOP_D_STORE_TIME();
            }
        }
    }

    if (programCount > 32) {
        START_D_TIME();
        val2 = shift(val1, 32);
        STOP_D_LOAD_TIME();
        if (programIndex % 64 == 0) {
            val1 = val1 * val2;
            if (programCount == 64) {
                START_D_TIME();
                y[programIndex] = val1;
                STOP_D_STORE_TIME();
            }
        }
    }

    STOP_THR_TIME();
}

/*
this is the code from the cs149 lecture slides:
export void vec8product(uniform float* x, uniform float* result) {

    float val1 = x[programIndex];
    float val2 = shift(val1, -1);

    if (programIndex % 2 == 0)
        val1 = val1 * val2;

    val2 = shift(val1, -2);
    if (programIndex % 4 == 0)
        val1 = val1 * val2;

    val2 = shift(val1, -4);
    if (programIndex % 8 == 0) {
        val1 = val1 * val2;
        *result = val1;
    }
}
*/

//  if compile with --addressing=32 --target=avx2-i32x8
//  this method's foreach is invoked twice to cover the 16 elements with an 8-wide vector index,
//  but each vec8product invocation has programIndex=[0,1,2,3,4,5,6,7] and does not
//  get incremented to [8,9,...], so sending the variable i too, which is specific to the instance.
//  not also, if compile with --addressing=32 --target=avx2-i32x16, all 16 elements are in
//  programIndex.
//  The programIndex is like a cuda thread block lane.

/*
memory allocation notes from ispc manual:

Assuming that a float is 4 bytes in memory.
let the gang size = 8 program instances.

uniform float * uniform p1 = uniform new uniform float[10];  // 4*10 = 40 bytes
float * uniform p2 = uniform new float[10];                  // 8*4*10 = 320 bytes
float * p3 = new float[10];                                  // 8*40 = 320 bytes
varying float * p4 = new varying float[10];                  // 8 * 320 = 2560 bytes
*/
export uniform float ispc_function(uniform int N, uniform float * uniform x) {

    INIT_TIME();

    uniform int n_res = N / programCount;

    //print("programCount=%, n_results=%\n", programCount, n_res);

    uniform float * uniform tmpX = x;
    //if N=65536, and use tmpX = copy x here instead of reference to x, runtime
    //      execution increases from 80 cycles to 270 cycles.
    //uniform float * uniform tmpX = uniform new float[N];
    //memcpy(tmpX, x, N*4);
    uniform float * uniform tmpY = uniform new float[programCount];
    uniform float * uniform res = uniform new float[n_res];

    int nIter = 0;
    // handles to res new arrays that are deleted after no longer used.
    uniform float * uniform prev0 = NULL;
    uniform float * uniform prev1 = NULL;

    //uniform float logN = log((uniform float)N)/log((uniform float)2);
    //int nIter = (int)(floor((logN - 8)/3));
    // print("logN=%, nIter=%\n", logN, nnIter);

    uniform int i_res = 0;

    while (true) {
       foreach (i = 0 ... N) { // parallel execution
          // increment pointer for tmpX
          if (i_res > 0) {
              tmpX += programCount;
          }

          vecprod(tmpX, tmpY);

          //print("(%) programCount=%, after tmpY=\n", i_res, programCount);
          //print_float_array(tmpY, 0, programCount);

          START_D_TIME();
          res[i_res] = tmpY[0];
          STOP_D_STORE_TIME();

          //print("i_res=%, programCount=%\n", i_res, programCount);
          //print_float_array(res, 0, n_res);

          i_res += 1;
      }
      if (n_res <= programCount) {
          break;
      }
      N = n_res;
      i_res = 0;
      n_res = N / programCount; 
    
      if (prev0 == NULL) {
          prev0 = res;
      } else if (prev1 == NULL) {
          prev1 = res;
      } else {
          delete[] prev0;
          prev0 = prev1;
          prev1 = res;    
      }
      tmpX = res;
      res = uniform new float[n_res]; 
      nIter += 1;
    }
   
    uniform float prod = 1.f;
    for (uniform int i = 0; i < n_res; ++i) {
        #pragma ignore warning(perf)
        prod *= res[i];
    }
    delete[] res;
    delete[] tmpY;
    if (prev0 != NULL && prev1 != NULL) {
        delete[] prev0;
    }
    return prod;
}

export uniform float ispc_higher_arith_func(uniform int N, uniform float * uniform x) {

    INIT_TIME();

    uniform float * uniform tmpX = x;
       
    uniform int i_res = 0;

    foreach (i = 0 ... N) { // parallel execution
        // increment pointer for tmpX
        if (i_res > 0) {
            tmpX += programCount;
        }

        higher_arith_intensity(tmpX);

        i_res += 1;
    }
    
    uniform float ret = 0;
    return ret;
}

uniform float ispc_function1(uniform int N, uniform float * uniform x) {

     INIT_TIME();

     uniform int n_res = N / programCount;

     //print("programCount=%, n_results=%\n", programCount, n_res);

     uniform float * uniform tmpX = x;
     uniform float * uniform tmpY = uniform new float[programCount];
     uniform float * uniform res = uniform new float[n_res];

     uniform int i0 = 0;
     uniform int i_res = 0;
     //for (uniform int i = 0; i < N; i += programCount) {
     foreach (i = 0 ... N) { // parallel execution
         // increment pointer for tmpX
         if (i_res > 0) {
             tmpX += programCount;
         }

         vecprod(tmpX, tmpY);

         //print("(%) programCount=%, after tmpY=\n", i_res, programCount);
         //print_float_array(tmpY, 0, programCount);

         START_D_TIME();
         res[i_res] = tmpY[0];
         STOP_D_STORE_TIME();

         //print("(%) res=\n", i_res);
         //print_float_array(res, 0, n_res);

         i_res += 1;
     }

     uniform float prod = 1.f;

     // trying a recursive approach increased the runtime execution
     if (n_res > 8) {
         prod = ispc_function1(n_res, res);
     } else {
         for (uniform int i = 0; i < n_res; ++i) {
             #pragma ignore warning(perf)
             prod *= res[i];
         }
     }
     delete[] res;
     delete[] tmpY;
     return prod;
}

uniform float _ispc_function(uniform int N, uniform float * uniform x) {

     INIT_TIME();

     uniform int n_res = N / programCount;

     //print("programCount=%, n_results=%\n", programCount, n_res);

     uniform float * uniform tmpX = x;
     uniform float * uniform tmpY = uniform new float[programCount];
     uniform float * uniform res = uniform new float[n_res];

     uniform int i0 = 0;
     uniform int i_res = 0;
     //for (uniform int i = 0; i < N; i += programCount) {
     foreach (i = 0 ... N) { // parallel execution
         // increment pointer for tmpX
         if (i_res > 0) {
             tmpX += programCount;
         }

         vecprod(tmpX, tmpY);

         //print("(%) programCount=%, after tmpY=\n", i_res, programCount);
         //print_float_array(tmpY, 0, programCount);

         START_D_TIME();
         res[i_res] = tmpY[0];
         STOP_D_STORE_TIME();

         //print("(%) res=\n", i_res);
         //print_float_array(res, 0, n_res);

         i_res += 1;
     }
     //res array could now be tmpX with N/=8 and repeating
     // above loop until N==8.
     uniform float prod = 1.f;
     for (uniform int i = 0; i < n_res; ++i) {
         #pragma ignore warning(perf)
         prod *= res[i];
     }
     delete[] res;
     delete[] tmpY;
     return prod;
}

/**
used to compare to ispc_function(...).
also useful to look at effect of programIndex filtering on what vecprod_1 sees.
*/
export uniform float ispc_function10(uniform int N, uniform float * uniform x) {

     uniform int n_res = N / programCount;

     //print("programCount=%, n_results=%\n", programCount, n_res);

     uniform float * uniform tmpX = x;
     uniform float * uniform res = uniform new float[n_res];

     uniform int i0 = 0;
     uniform int i_res = 0;
     
     foreach (i = 0 ... N) { // parallel execution
         // increment pointer for tmpX
         if (i_res > 0) {
             tmpX += programCount;
         }

         uniform int nIter = 0;
         uniform int nIterMax = 3; // for programCount==8
         if (programCount == 32) {
            nIterMax = 5;
         } else if (programCount == 16) {
            nIterMax = 4;
         }
         
         while (nIter < nIterMax) {
            uniform int ishift = (uniform int)1 << nIter;
            
            //if (programIndex % ishift == 0) { // Modulus operator with varying types is very inefficient.
            if ((nIter==0 && programIndex % 1 == 0) || 
            (nIter == 1 && programIndex % 2 == 0) ||
            (nIter == 2 && programIndex % 4 == 0) ||
            (nIter == 3 && programIndex % 8 == 0) ||
            (nIter == 4 && programIndex % 16 == 0)
            ) {
                //print("nIter=%, ishift=%\n", nIter, ishift);
                //print("(%)  before: ", nIter);
                //print_float_array(tmpX, 0, programCount);
                
                vecprod_1(tmpX, ishift);
                
                //print("(%)  after: ", nIter);
                //print_float_array(tmpX, 0, programCount);
            }
            ++nIter;
         }

         res[i_res] = tmpX[0];

         //print("(%) res=\n", i_res);
         //print_float_array(res, 0, n_res);

         i_res += 1;
     }
     uniform float prod = 1.f;
     for (uniform int i = 0; i < n_res; ++i) {
         #pragma ignore warning(perf)
         prod *= res[i];
     }
     delete[] res;
     return prod;
}


task void ispc_function_task(uniform int32 N, uniform float * uniform x, 
uniform int32 idx, uniform float * uniform prod_res) {

    INIT_TIME();

    uniform int n_res = N / programCount;

    //print("ispc_function_task: idx=%, programCount=%, n_results=%\n", idx, programCount, n_res);

    uniform float * uniform tmpX = x;
    //if N=65536, and use tmpX = copy x here instead of reference to x, runtime
    //      execution increases from 80 cycles to 270 cycles.
    //uniform float * uniform tmpX = uniform new float[N];
    //memcpy(tmpX, x, N*4);
    uniform float * uniform tmpY = uniform new float[programCount];
    uniform float * uniform res = uniform new float[n_res];

    int nIter = 0;
    // handles to res new arrays that are deleted after no longer used.
    uniform float * uniform prev0 = NULL;
    uniform float * uniform prev1 = NULL;

    //uniform float logN = log((uniform float)N)/log((uniform float)2);
    //int nIter = (int)(floor((logN - 8)/3));
    // print("logN=%, nIter=%\n", logN, nnIter);

    uniform int i_res = 0;

    while (true) {
       foreach (i = 0 ... N) { // parallel execution
          // increment pointer for tmpX
          if (i_res > 0) {
              tmpX += programCount;
          }

          vecprod(tmpX, tmpY);

          START_D_TIME();
          res[i_res] = tmpY[0];
          STOP_D_STORE_TIME();

          i_res += 1;
      }
      if (n_res <= programCount) {
          break;
      }
      N = n_res;
      i_res = 0;
      n_res = N / programCount; 
    
      if (prev0 == NULL) {
          prev0 = res;
      } else if (prev1 == NULL) {
          prev1 = res;
      } else {
          delete[] prev0;
          prev0 = prev1;
          prev1 = res;    
      }
      tmpX = res;
      res = uniform new float[n_res]; 
      nIter += 1;
    }

    #pragma ignore warning(perf)
    prod_res[idx] = res[0];
    for (uniform int i = 1; i < n_res; ++i) {
         #pragma ignore warning(perf)
         prod_res[idx] *= res[i];
    }

    delete[] res;
    delete[] tmpY;
    if (prev0 != NULL && prev1 != NULL) {
        delete[] prev0;
    }
}

export uniform float ispc_function_tasks(uniform int N, uniform float * uniform x) {
    
    uniform int N2 = N / 2;
    uniform int idx0 = 0;
    uniform int idx1 = 1;
    uniform float * uniform tmpX0 = x;
    uniform float * uniform tmpX1 = x += N2;

    uniform float * uniform prod_results = uniform new float[2];
    
    launch ispc_function_task(N2, tmpX0, idx0, prod_results);
    launch ispc_function_task(N2, tmpX1, idx1, prod_results);
    sync;

    //print("return from sync\n");

    uniform float prod = 1.f;
    prod = prod_results[0] * prod_results[1];
    delete[] prod_results;
    return prod;
}

/*
memory allocation notes from ispc manual:

Assuming that a float is 4 bytes in memory.
let the gang size = 8 program instances.

uniform float * uniform p1 = uniform new uniform float[10];  // 4*10 = 40 bytes
float * uniform p2 = uniform new float[10];                  // 8*4*10 = 320 bytes
float * p3 = new float[10];                                  // 8*40 = 320 bytes
varying float * p4 = new varying float[10];                  // 8 * 320 = 2560 bytes

*/
