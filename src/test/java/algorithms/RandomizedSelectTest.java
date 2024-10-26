package algorithms;

import algorithms.util.FormatArray;
import junit.framework.TestCase;

import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.Random;

public class RandomizedSelectTest extends TestCase {

    // algorithm implementation is not correct
    public void testRandMedian() {
        long seed = System.nanoTime();
        seed = 200341479853666L;
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);
        int n;
        int nTests = 100;
        double[] a;
        for (int nTest = 0; nTest < nTests; ++nTest) {
            n = 1 + rand.nextInt(1024); //set lower to make debugging easier
            a = new double[n];
            for (int i = 0; i < n; ++i) {
                a[i] = rand.nextInt() ;//* rand.nextDouble();
            }

            double[] b = Arrays.copyOf(a, a.length);
            Arrays.sort(b);
            double expAns = b[n/2];

            System.out.println("\n\nnTest=" + nTest);
            System.out.printf("a=%s\n", FormatArray.toString(a, "%.0f"));
            System.out.printf("b=%s\n", FormatArray.toString(b, "%.0f"));
            System.out.printf("n=%d, expected median = %.0f <===\n", n, expAns);

            double ans = RandomizedSelect.select(a, 0, a.length-1,
                    a.length/2, rand);

            System.out.printf("    result a=%s\n", FormatArray.toString(a, "%.0f"));
            System.out.printf("n=%d, expected median = %f <===\n", n, expAns);
            System.out.printf("    result median = %f\n", ans);
            {//DEBUG
                int idx = -1;
                for (int i = 0; i < n; ++i) {
                    if (Math.abs(ans - a[i]) < 1E-11) {
                        idx = i;
                        break;
                    }
                }
                System.out.printf("    result idx = %d, medianidx=%d\n", idx, n/2);
            }

            assertTrue(Math.abs(expAns - ans) < 1E-11);
        }
    }

    public void testQuickSort5() throws NoSuchAlgorithmException {
        //  0  1  2  3  4  5  6  7  8  9  10  11  12
        //  2        2        1        4          6
        //double[] a = new double[]{};
        double[] a = new double[]{0, 2, 0,0, 2, 0,0, 1, 0,0, 4, 0,0, 6, 0, 0};
        MedianOfMedians.quickSort5(a, 1, 13, 3);
        assertEquals(1., a[0+1]);
        assertEquals(2., a[3+1]);
        assertEquals(2., a[6+1]);
        assertEquals(4., a[9+1]);
        assertEquals(6., a[12+1]);
    }


}
