package algorithms.misc;

import junit.framework.TestCase;

import java.util.Arrays;
import java.util.Random;

public class ReservoirSamplingTest extends TestCase {

    public void test0() {

        // generate a Gumbel distribution and sample from it
        long seed = System.nanoTime();
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        int n = 10000;
        int k = 100;
        double[] a = new double[n];

        double loc = 4;
        double scale = 1;
        for (int i = 0; i < n; ++i) {
            a[i] = loc - scale * Math.log(Math.log(1./rand.nextDouble()));
        }

        double[] sample = ReservoirSampling.sampleOptimally(a, k);

        // compare the moments of a to sample
        double kMAD = 1.4826;

        double[] mADMinMax0 = MiscMath0.calculateMedianOfAbsoluteDeviation(a);
        double s0 = kMAD*mADMinMax0[0];
        double sigma0 = s0 * Math.sqrt(6.)/Math.PI;
        double[] q0 = MiscMath0.calculateQuartiles(a);
        double skew0 = (q0[2] + q0[0] - 2*q0[1])/(q0[2] - q0[0]);

        double[] mADMinMax1 = MiscMath0.calculateMedianOfAbsoluteDeviation(sample);
        double s1 = kMAD*mADMinMax0[0];
        double sigma1 = s0 * Math.sqrt(6.)/Math.PI;
        double[] q1 = MiscMath0.calculateQuartiles(sample);
        double skew1 = (q0[2] + q0[0] - 2*q0[1])/(q0[2] - q0[0]);

        //System.out.printf("mADMinMax0=%s\n", Arrays.toString(mADMinMax0));
        //System.out.printf("mADMinMax1=%s\n", Arrays.toString(mADMinMax1));
        //System.out.printf("skew0=%f\n", skew0);
        //System.out.printf("skew1=%f\n", skew1);

        assertTrue(Math.abs(sigma0 - sigma1) < 0.01);
        assertTrue(Math.abs(mADMinMax0[1] - mADMinMax1[1]) < sigma0);
        assertTrue(Math.abs(skew1 - skew0) < 1E-5);

        //TODO: consider adding a K-S test
    }
}
