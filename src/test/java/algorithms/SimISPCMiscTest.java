package algorithms;

import algorithms.misc.Shuffle;
import algorithms.util.FormatArray;
import junit.framework.TestCase;

import java.util.Arrays;

public class SimISPCMiscTest extends TestCase {
    /*
    implementing the ISPC example from stanford lecture 4 CS149 parallel computing.
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
    public void test0() throws Exception {

        int n = 16;
        int isWidth = 8; // instruction set width.  usually same as programCount
        float[] x = new float[n]; // loading the data "coherently"
        float expAns = 1.f;
        // simply using 10:25 for values
        for (int i = 0; i < n; ++i) {
            x[i] = (i + 10);
            expAns *= x[i];
        }

        float[] origX = Arrays.copyOf(x, x.length);

        float ans = mult(x, isWidth);
        assertTrue(Math.abs((expAns/ans)-1) <1E-6);

        //permute x to show order doesn't matter
        x =  Arrays.copyOf(origX, x.length);
        Shuffle.fisherYates(x);
        ans = mult(x, isWidth);
        assertTrue(Math.abs((expAns/ans)-1) <1E-6);
    }

    public float mult(float[] x, int isWidth) throws InterruptedException {
        if (x.length % isWidth != 0) throw new IllegalArgumentException("expecting x.length to be a multiple of " + isWidth);

        // isWidth == programCount
        int nInstances = x.length/isWidth;
        // emulating ISPC dividing the work among x.length/isWidth workers.
        /* NOTE that the ISPC gang instances all run in the same hardware thread and context.
        multiple instances of the function run in parallel vi multiple logical threads of control.
        there are no more than 2 program instances running in a gang - no more than twice the SIMD width of the hardware
        its running on.
            e.g. for a CPU running 4-wide SSE instruction set, there can be 4 or 8 program instances running in a gang.
            for a CPU running 8-wide AVX, there can be 8 or 16 program instances running in a gang.
        The SPMD parallelization is across SIMD lanes of a single core.

        Note: a gang is roughly equiv to a CUDA warp.
         */
        System.out.printf("nInstances=%d\n", nInstances);

        Thread[] threads = new Thread[nInstances];
        for (int i = 0; i < nInstances; ++i) {
            final int ii = i;
            threads[i] = new Thread(() -> mult(x, ii, isWidth));
        }

        // start a thread for each.  it helps to show that the order of execution doesn't matter.  there are no
        // dependencies between the work in each there here
        for (Thread thread : threads) {
            thread.start();
            thread.join();
        }

        //System.out.printf("res=%s\n", FormatArray.toString(x, "%.0f"));

        // when done with all threads, multiply results for each lane
        float res = 1.f;
        for (int i = isWidth-1; i < x.length; i += isWidth) {
            res *= x[i];
        }

        return res;
    }

    protected void mult(float[] x, int instanceNumber, int isWidth) {
        // x is a uniform shared variable array.

        // begin spmd replacement.  estimating the lane as offsets in the uniform shared variable x
        // idx0 is programIndex
        int idx0 = isWidth * instanceNumber;
        int idx1 = idx0 + isWidth - 1;//isWidth * (instanceNumber + 1) -1;

        // the span of this program is log_2(isWidth) = 3 for SPMD execution

        for (int prId = 2; prId <= isWidth; prId <<= 1) {
            // each prId is 1 round of computation on the vector of data between idx0 and idx1
            // we pretend the vector operations in each prId level is SPMD (single process multiple data) and is executed
            // simultaneously for this lane.

            int off0 = prId - 1;
            int off1 = (prId/2) - 1;

            //System.out.printf("instanceNumber=%d, prId=%d, idx0=%d, idx1=%d, off0=%d, off1=%d\n",
            //        instanceNumber, prId, idx0, idx1, off0, off1);

            // pretend this is SPMD
            for (int j = idx0; j < idx1; j+= prId) {
                //System.out.printf("   j+off0=%d, j+off1=%d\n", j+off0, j+off1);
                x[j + off0] *= x[j+off1];
            }
        }
    }


}
