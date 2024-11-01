package algorithms;

import java.util.Arrays;
import java.util.Random;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SerialCompareTest extends TestCase {
    
    public void testRand65536SerialProduct() {
	System.out.println("testRand65536SerialProduct");
	int nTests = 10;
	double start;
	double[] times = new double[nTests];
	for (int i = 0; i < nTests; ++i) {
	    float[] x = generateRandomArray(1 << 16);
            start = System.nanoTime();
            serialProduct(x);
            times[i] = System.nanoTime() - start;
	}

        double med = MedianOfMediansSelect.selectCLRS(times, 0, times.length-1, times.length/2);

	System.out.printf("median: %f\n", med);
    }
    public void testEstimateProcessorSpeed() {
	System.out.println("testEstimateProcessorSpeed");
	int n = 1_000_000;
	int nTests = 20;
	double[] a = generateRandomArrayD(n);
	double[] times = new double[nTests];
	double start;
	for (int i = 0; i < nTests; ++i) {
	    start = System.nanoTime();
            for (int j = 0; j < n; ++j) {
    	        Math.sqrt(a[j]);
            }
	    times[i] = System.nanoTime() - start;
        }
	for (int i = 0; i < nTests; ++i) {
		System.out.printf("time=%f\n", times[i]);
	}
	double med = MedianOfMediansSelect.selectCLRS(times, 0, times.length-1, times.length/2);
	System.out.printf("n=%d instructions in %f nsec\n", n, med);

    }

    public void testRand250000SerialHighArith() {
	System.out.println("testRand250000SerialHighArith");
	int nTests = 10;
	double start;
	double[] times = new double[nTests];
	for (int i = 0; i < nTests; ++i) {
	    float[] x = generateRandomArray(250000);
            start = System.nanoTime();
            serialHigherArithmeticIntensity(x);
            times[i] = System.nanoTime() - start;
	}

        double med = MedianOfMediansSelect.selectCLRS(times, 0, times.length-1, times.length/2);

	System.out.printf("median: %f\n", med);
    }

    protected float[] generateRandomArray(int n) {
         // scaled so that product fits into Float.MAX_VALUE;
	 long seed = System.nanoTime();
	 //System.out.printf("seed=%d\n", seed);
	 Random rand = new Random(seed);

	 float factor = (float)Math.pow(Float.MAX_VALUE, 1.f/n);
	 float[] a = new float[n];
	 for (int i = 0; i < n; ++i) {
	     a[i] = 1.0f + factor * rand.nextFloat();
	 }
	 return a;
    }

    protected double[] generateRandomArrayD(int n) {
         // scaled so that product fits into Float.MAX_VALUE;
	 long seed = System.nanoTime();
	 //System.out.printf("seed=%d\n", seed);
	 Random rand = new Random(seed);

	 double factor = Math.pow(Float.MAX_VALUE, 1./n);
	 double[] a = new double[n];
	 for (int i = 0; i < n; ++i) {
	     a[i] = 1.0 + factor * rand.nextFloat();
	 }
	 return a;
    }

    protected void serialProduct(float[] a) {
	float prod = 1.f;
	for (float aI : a) {
	    prod *= aI; 
        }
    }

    protected void serialHigherArithmeticIntensity(float[] a) {
	// keeping in arrays for comparison with intrinsics methods
        float[] factors = {3.f, 3.f, 3.f, 3.f};
        float[] s = {1.5f, 1.5f, 1.5f, 1.5f};
        int n = a.length;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < 4; ++j) {
		a[i] = (a[i] * factors[0]) + s[i%2];
		a[i] = 1.f/a[i];
		a[i] = (float)Math.sqrt(a[i]);
            }
        }
    }
}
