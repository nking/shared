/*
https://ispc.github.io/ispc.html#uniform-and-varying-qualifiers
If a variable has a uniform qualifier, then there is only a single instance of that variable shared by all program
instances in a gang. (In other words, it necessarily has the same value across all of the program instances.)
In addition to requiring less storage than varying values, uniform variables lead to a number of performance advantages
when they are applicable (see Uniform Control Flow, for example.) Varying variables may be qualified with varying,
though doing so has no effect, as varying is the default.

to compile:
   ispc ispc_function.ispc -o ispc_function.o
   see CMakeLists.txt
*/

void vec8product(float * uniform x, float * uniform y, int i /*instance i*/, int iShift) {
    // TODO: add time spans to see if the "gather" is a single inefficient fetch from memory to cache.
    // also browse what this is: packed_load_active
    float val1 = x[i];
    float val2 = x[i+iShift];

    //print("\nprogramIndex=%\n", programIndex);
    //print("iShift=%, i=%\n", iShift, i);

    //print("val1=%\n", val1);
    //print("val2=%\n", val2);

    y[i] = val1 * val2;

    //float val5 = y[i];
    //print("after y        =%\n", val5);
    //print("x address=%\n", x);
    //print("y address=%\n", y);
}

void vec16product(float * uniform x, float * uniform y, int iShift) {

    float val1 = x[programIndex];
    float val2 = x[programIndex + iShift];//shift(val1, iShift);
    y[programIndex] = val1 * val2;

    //print("programCount=%\n", programCount);
    //print("programIndex=%\n", programIndex);
    //print("val1=%\n", val1);
    //print("val2=%\n", val2);
}

/*
this is the code from the cs149 lecture slides:
export void vec8product(uniform float* x, uniform float* result) {

    float val1 = x[i + programIndex];
    float val2 = shift(val1, -1);

    if (programIndex % 2 == 0)
        val1 = val1 * val2;

    val2 = shift(val1, -2);
    if (programIndex % 4 == 0)
        val1 = val1 * val2;

    val2 = shift(val1, -4);
    if (programIndex % 8 == 0)
        val1 = val1 * val2;
        *result = val1;
}
*/

//  if compile with --addressing=32 --target=avx2-i32x8
//  this method's foreach is invoked twice to cover the 16 elements with an 8-wide vector lane.
//  but each vec8product invocation has programIndex=[0,1,2,3,4,5,6,7] and does not
//  get incremented to [8,9,...], so sending the variable i too, which is specific to the instance.
//  not also, if compile with --addressing=32 --target=avx2-i32x16, all 16 elements are in
//  programIndex
//export float ispc_function(uniform int N, uniform float* x, uniform float* y) {
export uniform float ispc_function(uniform int N,  float * uniform x) {
//export float ispc_function(uniform int N,  uniform float x[]) {

   // space consumed by x is N*4*8? because it is not uniform float * uniform?
   //uniform float * uniform tmpX = x;
   float * uniform tmpX = x;
   float * uniform tmpY = uniform new float[N];

   //print("programCount=%, N=%, expecting nInstances=%\n", programCount, N, N / programCount);

   uniform int nIter = 0;
   while (nIter < 4) {
       foreach (i=0...N) { // parallel execution
           // each gang process span is 1
           if (programCount==8) {
               vec8product(tmpX, tmpY, i, 1<<nIter);
           } else if (programCount==16 && programCount == N) {
               vec16product(tmpX, tmpY, 1<<nIter);
           }
       }
       float * uniform tmp = tmpX;
       tmpX = tmpY;
       tmpY = tmp;
       ++nIter;
       //print("nIter=%\n", nIter);
   }
   return tmpX[0];
}

/*
memory allocation notes from ispc manual:

Assuming that a float is 4 bytes in memory.
let the gang size = 8 program instances.

uniform float * uniform p1 = uniform new uniform float[10];  // 4*10 = 40 bytes
float * uniform p2 = uniform new float[10];                  // 8*4*10 = 320 bytes
float * p3 = new float[10];                                  // 8*40 = 30 bytes
varying float * p4 = new varying float[10];                  // 8 * 320 = 2560 bytes

*/
