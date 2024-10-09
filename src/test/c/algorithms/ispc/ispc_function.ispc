/*
https://ispc.github.io/ispc.html#uniform-and-varying-qualifiers
If a variable has a uniform qualifier, then there is only a single instance of that variable shared by all program
instances in a gang. (In other words, it necessarily has the same value across all of the program instances.)
In addition to requiring less storage than varying values, uniform variables lead to a number of performance advantages
when they are applicable (see Uniform Control Flow, for example.) Varying variables may be qualified with varying,
though doing so has no effect, as varying is the default.

to compile:
   ispc ispc_function.ispc -o ispc_function.o
*/

//export void vec8product(uniform float* x, uniform float* result) {
void vec8product(float * uniform x, float * uniform result) {
    // x reference and result reference should be incremented by ispc handler
    // and programIndex and programCount are available to this method for use.

    //NLK: error in my invocation of this possibly, or in compile target
    // programIndex is the same array for both invocations here.

    float val1 = x[programIndex];
    float val2 = shift(val1, -1);

    print("programIndex=%\n", programIndex);
    print("val1=%\n", val1);
    print("val2=%\n", val2);

    if (programIndex % 2 == 0) {
        //print("   mod 2: val1 * val2=%\n", val1 * val2);
        val1 = val1 * val2;
    }

    val2 = shift(val1, -2);
    if (programIndex % 4 == 0) {
        //print("   mod 4: val1 * val2=%\n", val1 * val2);
        val1 = val1 * val2;
    }

    val2 = shift(val1, -4);
    if (programIndex % 8 == 0) {
        //*result = val1 * val2;
        val1 = val1 * val2;
    }
    //print("=> x=%\n", x[programIndex]);
    //print("=> val1=%\n", val1);
    result[programIndex] = val1;
}

//export void ispc_function(uniform int N, uniform float* x, uniform float* y) {
export void ispc_function(uniform int N,  float * uniform x, float * uniform y) {
   //int programCount = N/isWidth;
   for (uniform int i = 0; i < N; i += programCount) {
       vec8product(x, y);
   }
   /*
   foreach (i=0...N) { // parallel execution
       // why isn't this i == programIndex in vec8product?
       print("invoke vec8product i=%  programCount=%\n", i, programCount);
       vec8product(x, y);
   }*/
}

