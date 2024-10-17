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

void print_float_array(uniform float *array, uniform int i0, uniform int i1) {
    for (uniform int i = i0; i < i1; ++i) {
        #pragma ignore warning(perf)
        print("   %,\n", array[i]);
    }
    print("\n");
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
     uniform float * uniform res = uniform new float[n_res];

     uniform int i0 = 0;
     uniform int i_res = 0;
     for (uniform int i = 0; i < N; i += programCount) {

         // increment pointer for tmpX
         if (i > 0) {
             tmpX += programCount;
         }
         uniform float * uniform tmpY = uniform new float[programCount];

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
     for (uniform int i = 0; i < n_res; ++i) {
         #pragma ignore warning(perf)
         prod *= res[i];
     }
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
