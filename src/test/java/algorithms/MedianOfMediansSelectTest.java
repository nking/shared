package algorithms;

import junit.framework.TestCase;

import java.security.NoSuchAlgorithmException;
import java.util.Random;
import java.util.Arrays;
import algorithms.util.FormatArray;

public class MedianOfMediansSelectTest extends TestCase {

    // TODO: add random test for ranks other than the median

    public void testPartition() {

        double[] a2 = new double[]{4,2,8,1,0, 3, 3};
        int av2 = MedianOfMediansSelect.partitionAround(a2, 0, a2.length-1, 3);
        assert(a2[av2] == 3);

        double[] b2 = new double[]{3, 3, 2,4,8,1,0};
        int bv2 = MedianOfMediansSelect.partitionAround(b2, 0, b2.length-1, 0);
        assert(b2[bv2] == 0);

        double[] b3 = new double[]{3, 3, 2,4,8,1,0};
        int bv3 = MedianOfMediansSelect.partitionAround(b3, 0, b3.length-1, 1);
        assert(b3[bv3] == 1);

        double[] b4 = new double[]{3, 3, 2,4,8,1,0};
        int bv4 = MedianOfMediansSelect.partitionAround(b4, 0, b4.length-1, 8);
        assert(b4[bv4] == 8);

        double[] b5 = new double[]{3, 3, 2,4,8,1,0};
        int bv5 = MedianOfMediansSelect.partitionAround(b5, 0, b5.length-1, 4);
        assert(b5[bv5] == 4);

        double[] b6 = new double[]{3, 3, 2,4,8,1,0};
        int bv6 = MedianOfMediansSelect.partitionAround(b6, 0, b6.length-1, 2);
        assert(b6[bv6] == 2);

        double[] b7 = new double[]{3, 3, 2,4,8,1,0};
        int bv7 = MedianOfMediansSelect.partitionAround(b7, 0, b7.length-1, 3);
        assert(b7[bv7] == 3);
    }

    // algorithm implementation is not correct
    public void testRand() {
        long seed = System.nanoTime();
        //seed = 201396339917283L;
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

            /*System.out.println("\n\nnTest=" + nTest);
            System.out.printf("a=%s\n", FormatArray.toString(a, "%.0f"));
            System.out.printf("b=%s\n", FormatArray.toString(b, "%.0f"));
            System.out.printf("n=%d, expected median = %.0f <===\n", n, expAns);
            */

            double ans = MedianOfMediansSelect.selectCLRS(a, 0, a.length-1, a.length/2);

            //System.out.printf("    result a=%s\n", FormatArray.toString(a, "%.0f"));
            //System.out.printf("    result median = %f\n", ans);

            assertTrue(Math.abs(expAns - ans) < 1E-11);
        }
    }

    public void testQuickSort5() throws NoSuchAlgorithmException {
        //  0  1  2  3  4  5  6  7  8  9  10  11  12
        //  2        2        1        4          6
        //double[] a = new double[]{};
        double[] a = new double[]{0, 2, 0,0, 2, 0,0, 1, 0,0, 4, 0,0, 6, 0, 0};
        MedianOfMediansSelect.quickSort5(a, 1, 13, 3);
        assertEquals(1., a[0+1]);
        assertEquals(2., a[3+1]);
        assertEquals(2., a[6+1]);
        assertEquals(4., a[9+1]);
        assertEquals(6., a[12+1]);
    }


}
