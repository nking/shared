package algorithms;

import junit.framework.TestCase;

import java.util.Arrays;
import java.util.Random;

public class RandomizedQuickSortTest extends TestCase {

    public void test0() {

        long seed = System.nanoTime();
        //seed = 151408698946400L;
        System.out.printf("seed=%d\n", seed);

        Random rand = new Random(seed);
        int n;
        int nTests = 100;
        double[] a;
        for (int nTest = 0; nTest < nTests; ++nTest) {
            n = 1 + rand.nextInt(1024);
            a = new double[n];
            for (int i = 0; i < n; ++i) {
                a[i] = rand.nextInt();//* rand.nextDouble();
            }

            double[] b = Arrays.copyOf(a, a.length);
            Arrays.sort(b);

            RandomizedQuickSort.sort(a, 0, n-1, rand);

            for (int i = 0; i < n; ++i) {
                assertTrue(Math.abs(a[i] - b[i]) < 1E-11);
            }
        }
    }
}
