package algorithms;

import algorithms.misc.Shuffle;
import junit.framework.TestCase;

import java.security.NoSuchAlgorithmException;

public class MedianOfMediansTest extends TestCase {

    public void testPartition() {

        double[] a2 = new double[]{4,2,8,1,0, 3, 3};
        int av2 = MedianOfMedians.partitionAround(a2, 0, a2.length-1, 3);
        assert(a2[av2] == 3);

        double[] b2 = new double[]{3, 3, 2,4,8,1,0};
        int bv2 = MedianOfMedians.partitionAround(b2, 0, b2.length-1, 0);
        assert(b2[bv2] == 0);

        double[] b3 = new double[]{3, 3, 2,4,8,1,0};
        int bv3 = MedianOfMedians.partitionAround(b3, 0, b3.length-1, 1);
        assert(b3[bv3] == 1);

        double[] b4 = new double[]{3, 3, 2,4,8,1,0};
        int bv4 = MedianOfMedians.partitionAround(b4, 0, b4.length-1, 8);
        assert(b4[bv4] == 8);

        double[] b5 = new double[]{3, 3, 2,4,8,1,0};
        int bv5 = MedianOfMedians.partitionAround(b5, 0, b5.length-1, 4);
        assert(b5[bv5] == 4);

        double[] b6 = new double[]{3, 3, 2,4,8,1,0};
        int bv6 = MedianOfMedians.partitionAround(b6, 0, b6.length-1, 2);
        assert(b6[bv6] == 2);

        double[] b7 = new double[]{3, 3, 2,4,8,1,0};
        int bv7 = MedianOfMedians.partitionAround(b7, 0, b7.length-1, 3);
        assert(b7[bv7] == 3);
    }

    // TODO: use this and random when have transformed algorithm to 0-based indexes
    public void _test0() throws NoSuchAlgorithmException {
        double[] a = new double[]{1,2,3, 4, 5,6,7};
        Shuffle.fisherYates(a);

        double expAns = 4;
        double ans = MedianOfMedians.select(a,0,a.length-1, 3);
        assertTrue(Math.abs(expAns - ans) < 1E-11);
    }

    public void test1() throws NoSuchAlgorithmException {
        //double[] a = new double[]{0,1,2,2,2,2,2,3,4,5,6,7};
        // temporarily buffer with 0s on ends until transform to 0-based indexes
        double[] a = new double[]{0,1,2,2,2,2,2,3,4,5,6,7, 1,2,2,2,2,2,3,4,5,6,7, 1,2,2,2,2,2,3,4,5,6,7,0};
        Shuffle.fisherYates(a);

        double expAns = 2;
        double ans = MedianOfMedians.select(a,1,a.length-1, a.length/2);
        assertTrue(Math.abs(expAns - ans) < 1E-11);
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
