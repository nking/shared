package algorithms.misc;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.logging.Logger;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MiscMath0Test extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    public void testGetAvgAndStDev() throws Exception {
        
        int[] x = new int[]{2, 2, 2, 4, 4, 4};
        
        double[] avgAndStDev = MiscMath0.getAvgAndStDev(x, x.length);
        assertEquals(avgAndStDev[0], 3.);
        assertTrue(Math.abs(avgAndStDev[1] - 1.1) < 0.01);
        
        avgAndStDev = MiscMath0.getAvgAndStDev(x);
        assertEquals(avgAndStDev[0], 3.);
        assertTrue(Math.abs(avgAndStDev[1] - 1.1) < 0.01);
    }
    
    public void testGetAvgAndStDev2() throws Exception {
        
        double[] x = new double[]{2, 2, 2, 4, 4, 4};
        
        double[] avgAndStDev = MiscMath0.getAvgAndStDev(x);
        assertEquals(avgAndStDev[0], 3.);
        assertTrue(Math.abs(avgAndStDev[1] - 1.1) < 0.01);
        
    }
    
    public void testGet20NeighborOffsets() throws Exception {
                
        PairIntArray offsets = MiscMath0.get20NeighborOffsets();
        
        assertTrue(offsets.getN() == 20);
        
        Set<PairInt> expected = new HashSet<PairInt>();
        for (int dx = -2; dx <= 2; ++dx) {
            for (int dy = -2; dy <= 2; ++dy) {
                if (dx == 0 && dy == 0) {
                    continue;
                }
                expected.add(new PairInt(dx, dy));
            }
        }
        //remove the 4 corners        
        expected.remove(new PairInt(-2, -2));
        expected.remove(new PairInt(-2, 2));
        expected.remove(new PairInt(2, -2));
        expected.remove(new PairInt(2, 2));
        
        assertTrue(expected.size() == 20);
        
        for (int i = 0; i < offsets.getN(); ++i) {
            PairInt p = new PairInt(offsets.getX(i), offsets.getY(i));
            assertTrue(expected.remove(p));
        }
        
        assertTrue(expected.isEmpty());
    }

    /**
     * Test of findPowerOf10 method, of class MiscMath0.
     */
    public void testFindPowerOf10() {

        assertTrue(MiscMath0.findPowerOf10(0.1f) == -1);

        assertTrue(MiscMath0.findPowerOf10(0.11f) == -1);

        assertTrue(MiscMath0.findPowerOf10(0.55f) == -1);

        int p0;
        int p1;
        p0 = MiscMath0.findPowerOf10(0.099999999f);
        p1 = MiscMath0.findPowerOf10_2(0.099999999f);
        //assertTrue(MiscMath0.findPowerOf10(0.099999999f) == -2); <=== precision error before argument reaches method.

        p0 = MiscMath0.findPowerOf10(0.01f);
        p1 = MiscMath0.findPowerOf10_2(0.01f);
        assertTrue(p0 == -2);

        assertTrue(MiscMath0.findPowerOf10(0.011f) == -2);

        assertTrue(MiscMath0.findPowerOf10(0.001f) == -3);

        //assertTrue(MiscMath0.findPowerOf10(0.0099999999f) == -3); <=== precision error before argument reaches method.

        assertTrue(MiscMath0.findPowerOf10(0.f) == 0);

        assertTrue(MiscMath0.findPowerOf10(1.f) == 0);

        assertTrue(MiscMath0.findPowerOf10(10.f) == 1);

        assertTrue(MiscMath0.findPowerOf10(11.f) == 1);

        assertTrue(MiscMath0.findPowerOf10(100.f) == 2);

        assertTrue(MiscMath0.findPowerOf10(-3.1f) == 0);

        assertTrue(MiscMath0.findPowerOf10(-31.1f) == 1);

        assertTrue(MiscMath0.findPowerOf10(-310.1f) == 2);

        assertTrue(MiscMath0.findPowerOf10(-0.1f) == -1);
    }

    /**
     * Test of roundDownByLargestPower method, of class MiscMath0.
     */
    public void testRoundDownByLargestPower() {

        assertTrue(MiscMath0.roundDownByLargestPower(0.1f) == 0.1f);

        assertTrue(MiscMath0.roundDownByLargestPower(0.11f) == 0.1f);

        assertTrue(MiscMath0.roundDownByLargestPower(3.1f) == 3.0f);

        assertTrue(MiscMath0.roundDownByLargestPower(31.1f) == 30.0f);

        assertTrue(MiscMath0.roundDownByLargestPower(31.0f) == 30.0f);

        assertTrue(MiscMath0.roundDownByLargestPower(-3.1f) == -4.0f);

        assertTrue(MiscMath0.roundDownByLargestPower(-31.1f) == -40.0f);

        assertTrue(MiscMath0.roundDownByLargestPower(-31.0f) == -40.0f);
    }

    /**
     * Test of roundUpByExponent method, of class MiscMath0.
     */
    public void testRoundUpByLargestPower() {

        assertTrue(MiscMath0.roundUpByLargestPower(31.1f) == 40.0f);

        assertTrue(MiscMath0.roundUpByLargestPower(0.11f) == 0.2f);

        assertTrue(MiscMath0.roundUpByLargestPower(-0.011f) == -0.02f);

        assertTrue(MiscMath0.roundUpByLargestPower(-0.11f) == -0.2f);

        assertTrue(MiscMath0.roundUpByLargestPower(310.1f) == 400.0f);

        assertTrue(MiscMath0.roundUpByLargestPower(-3.1f) == -4.0f);

        assertTrue(MiscMath0.roundUpByLargestPower(3.1f) == 4.0f);

        assertTrue(MiscMath0.roundUpByLargestPower(10.0f) == 10.0f);

        assertTrue(MiscMath0.roundUpByLargestPower(5.0f) == 5.0f);
    }

    public void test0() throws Exception {
        // constructor isn't used, this is to complete the coverage.
        MiscMath0 mm = new MiscMath0();
    }

    public void testFindStrongestPeakIndexes() throws Exception {

        float[] xHist = new float[]{
            -279.96667f, -239.90002f, -199.83334f, -159.76666f,
            -119.700005f, -79.63334f, -39.566673f, 0.5000076f,
            40.566658f, 80.63331f, 120.69999f, 160.76666f,
            200.83331f, 240.89996f, 280.96667f};

        int[] yHist = new int[] {0, 24, 141, 259, 482, 894, 1110, 1239,
            890, 608, 302, 151, 56, 0, 0};

        HistogramHolder h = new HistogramHolder();
        h.setYHist(yHist);
        h.setXHist(xHist);

        List<Integer> peakIndexes = MiscMath0.findStrongPeakIndexes(h, 0.1f);

        assertNotNull(peakIndexes);
        assertTrue(peakIndexes.size() == 1);
        assertTrue(peakIndexes.get(0).intValue() == 7);

        List<Integer> sortedPeakIndexes =
            MiscMath0.findStrongPeakIndexesDescSort(h, 0.09f);
        assertNotNull(sortedPeakIndexes);
        assertTrue(sortedPeakIndexes.size() == 1);
        assertTrue(sortedPeakIndexes.get(0).intValue() == 7);

        // this should have 3 peaks
        xHist = new float[]{
            -466.63333f, -399.9f, -333.16666f, -266.43335f, -199.70001f,
            -132.96667f, -66.23337f, 0.49996567f, 67.23331f, 133.96664f,
            200.69998f, 267.43332f, 334.1666f, 400.89993f, 467.63327f
        };

        yHist = new int[] {
            115, 265, 385, 560, 484, 581, 748, 809,
            661, 375, 371, 379, 304, 109, 10
        };

        h = new HistogramHolder();
        h.setYHist(yHist);
        h.setXHist(xHist);

        peakIndexes = MiscMath0.findStrongPeakIndexes(h, 0.09f);

        assertNotNull(peakIndexes);
        assertTrue(peakIndexes.size() == 2);
        assertTrue(peakIndexes.get(0).intValue() == 3);
        assertTrue(peakIndexes.get(1).intValue() == 7);

        sortedPeakIndexes = MiscMath0.findStrongPeakIndexesDescSort(h, 0.09f);

        assertNotNull(sortedPeakIndexes);
        assertTrue(sortedPeakIndexes.size() == 2);
        assertTrue(sortedPeakIndexes.get(0).intValue() == 7);
        assertTrue(sortedPeakIndexes.get(1).intValue() == 3);
    }

    public void testFindMinMaxXY_4() {
        
        /*
         2    #     #
         1    #
         0 #
           0  1  2  3  
        */
        int w = 5;
        TIntSet pixIdxs = new TIntHashSet();
        pixIdxs.add((0 * w) + 0);
        pixIdxs.add((1 * w) + 1);
        pixIdxs.add((2 * w) + 1);
        pixIdxs.add((2 * w) + 3);
        
        int[] minmaxXY = MiscMath0.findMinMaxXY(pixIdxs, w);
        
        assertEquals(0, minmaxXY[0]);
        assertEquals(3, minmaxXY[1]);
        assertEquals(0, minmaxXY[2]);
        assertEquals(2, minmaxXY[3]);
    }

    public void testFindYMaxIndex_list() {
        
        int[] a = new int[]{1,5,7,2,4,8,3};
        TIntList input = new TIntArrayList(a);
        
        int maxIdx = MiscMath0.findYMaxIndex(input);
        assertEquals(5, maxIdx);
        
        double[] a1 = new double[]{1,5,7,2,4,8,3};
        maxIdx = MiscMath0.findYMaxIndex(a1);
        assertEquals(5, maxIdx);
        
        float[] a2 = new float[]{1,5,7,2,4,8,3};
        maxIdx = MiscMath0.findYMaxIndex(a2);
        assertEquals(5, maxIdx);
        
    }
    
    public void testIsAPowerOf2() throws Exception {

        assertTrue(MiscMath0.isAPowerOf2(2));

        assertTrue(MiscMath0.isAPowerOf2(1<<2));

        assertTrue(MiscMath0.isAPowerOf2(1<<7));

        assertFalse(MiscMath0.isAPowerOf2(71));

        assertFalse(MiscMath0.isAPowerOf2(17));
    }

    public void testNumberOfBits() throws Exception {

        for (int i = -256; i >= 0; i--) {
            String bitstring = Integer.toBinaryString(i);
            int expected = bitstring.length();
            int result = MiscMath0.numberOfBits(i);
            assertTrue(expected == result);
        }

        for (int i = 0; i <= 256; i++) {
            String bitstring = Integer.toBinaryString(i);
            int expected = bitstring.length();
            int result = MiscMath0.numberOfBits(i);
   // System.out.println(i + ") " + bitstring 
   // + " exp=" + expected + " r=" + result);
            assertTrue(expected == result);
        }
    }

    public void testBitReverse() throws Exception {

        int max = 1 << 4;
        int nBits = MiscMath0.numberOfBits(max);

        for (int i = 0; i <= max; i++) {

            String bitstring = Integer.toBinaryString(i);
            while (bitstring.length() < nBits) {
                bitstring = "0" + bitstring;
            }
            char[] revBitstring = bitstring.toCharArray();
            for (int ii = 0; ii < (bitstring.length()/2); ii++) {
                int idx2 = bitstring.length() - ii - 1;
                char swap = revBitstring[ii];
                revBitstring[ii] = revBitstring[idx2];
                revBitstring[idx2] = swap;
            }

            String rBitstring = String.valueOf(revBitstring);

            int rev = Integer.parseInt(rBitstring, 2);

            int result = MiscMath0.bitReverse(i, nBits);

            assertTrue(rev == result);
        }

    }
    
    public void testFindMinMax_int_array_2d() {
        
        int[][] a = new int[3][3];
        a[0] = new int[]{1, 2, 3};
        a[1] = new int[]{1, 100, 1};
        a[2] = new int[]{-1, -2, -3};
        
        int max = MiscMath0.findMax(a);
        
        assertEquals(100, max);
        
        int min = MiscMath0.findMin(a);
        
        assertEquals(-3, min);
    }
    
    public void testFindMinMax_float_array_2d() {
        
        float[][] a = new float[3][3];
        a[0] = new float[]{1, 2, 3};
        a[1] = new float[]{1, 100, 1};
        a[2] = new float[]{-1, -2, -3};
        
        float max = MiscMath0.findMax(a);
        
        assertEquals(100.f, max);
        
        float min = MiscMath0.findMin(a);
        
        assertEquals(-3.f, min);
    }
    
    public void testCalcQuartilesTest_float() {
        
        // 0 1 2  4 4 4  4 6 8  12 12 16
        //   3      12     18      40
        //  3/73    ...
        
        float[] a = new float[]{0, 1, 2, 4, 4, 4, 4, 6, 8, 12, 12, 16};
        
        float[] expected = new float[]{
            3.f/73.f, 12.f/73.f, 18.f/73.f, 40.f/73.f};
    
        float[] b = MiscMath0.calcQuartiles(a, true);
        
        assertTrue(Arrays.equals(b, expected));
        
        a = new float[]{4, 4, 4, 4, 0, 1, 2, 6, 8, 12, 12, 16};
        
        b = MiscMath0.calcQuartiles(a, false);
        
        assertTrue(Arrays.equals(b, expected));
    }
    
    public void testCalcQuartilesTest_double() {
        
        // 0 1 2  4 4 4  4 6 8  12 12 16
        //   3      12     18      40
        //  3/73    ...
        
        double[] a = new double[]{0, 1, 2, 4, 4, 4, 4, 6, 8, 12, 12, 16};
        
        double[] expected = new double[]{
            3./73., 12./73., 18./73., 40./73.};
    
        double[] b = MiscMath0.calcQuartiles(a, true);
        
        assertTrue(Arrays.equals(b, expected));
        
        a = new double[]{4, 4, 4, 4, 0, 1, 2, 6, 8, 12, 12, 16};
        
        b = MiscMath0.calcQuartiles(a, false);
        
        assertTrue(Arrays.equals(b, expected));
    }
    
    public void testCalcQuartilesTest_int() {
        
        // 0 1 2  4 4 4  4 6 8  12 12 16
        //   3      12     18      40
        //  3/73    ...
        
        float norm = 73.f;
        
        int[] a = new int[]{0, 1, 2, 4, 4, 4, 4, 6, 8, 12, 12, 16};
        
        float[] expected = new float[]{
            (3.f/73.f), (12.f/73.f), (18.f/73.f), (40.f/73.f)};
    
        float[] b = MiscMath0.calcQuartiles(a, true);
        
        assertTrue(Arrays.equals(b, expected));
        
        a = new int[]{4, 4, 4, 4, 0, 1, 2, 6, 8, 12, 12, 16};
        
        b = MiscMath0.calcQuartiles(a, false);
        
        assertTrue(Arrays.equals(b, expected));
    }
    
    public void testMean() {
        
        double[] data;
        int nDimensions;
        double[] expected, c;
        double tol = 0.001;
        int nc;
        
        //(1, 1), (3, 3)                           --> (2,2)        
        nDimensions = 2;
        data = new double[]{1, 1,  3, 3};
        expected = new double[]{2, 2};
        
        c = MiscMath0.mean(data, nDimensions);
        assertEquals(expected.length, c.length);
        nc = 0;
        for (int i = 0; i < c.length; ++i) {
            if (Math.abs(expected[i] - c[i]) < tol) {
                nc++;
            }
        }
        assertEquals(expected.length, nc);
        
        //(0, 0), (0, 0), (0, 12)    
        nDimensions = 2;
        data = new double[]{0, 0, 0, 0, 0, 12};
        expected = new double[]{0, 4};
        
        c = MiscMath0.mean(data, nDimensions);
        assertEquals(expected.length, c.length);
        nc = 0;
        for (int i = 0; i < c.length; ++i) {
            if (Math.abs(expected[i] - c[i]) < tol) {
                nc++;
            }
        }
        assertEquals(expected.length, nc);
        
        //(-20, 48), (-20, -48), (20, 0), (59, 0)
        nDimensions = 2;
        data = new double[]{-20, 48, -20, -48, 20, 0, 59, 0};
        expected = new double[]{9.75, 0};
        
        c = MiscMath0.mean(data, nDimensions);
        assertEquals(expected.length, c.length);
        nc = 0;
        for (int i = 0; i < c.length; ++i) {
            if (Math.abs(expected[i] - c[i]) < tol) {
                nc++;
            }
        }
        assertEquals(expected.length, nc);
    }
    
    public void testStandardDeviation() {
        
        double[] data;
        int nDimensions;
        double[][] expected, c;
        expected = new double[2][2];
        double tol = 0.001;
        int nc;
        
        //(1, 1), (3, 3)                           --> (2,2)        
        nDimensions = 2;
        data = new double[]{1, 1,  3, 3};
        expected[0] = new double[]{2, 2};
        expected[1] = new double[]{1.414, 1.414};
        
        c = MiscMath0.standardDeviation(data, nDimensions);
        assertEquals(expected.length, c.length);
        assertEquals(expected[0].length, c[0].length);
        assertEquals(expected[1].length, c[1].length);
        nc = 0;
        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[i].length; ++j) {
                if (Math.abs(expected[i][j] - c[i][j]) < tol) {
                    nc++;
                }
            }
        }
        assertEquals(expected.length*expected[0].length, nc);
        
        //(0, 0), (0, 0), (0, 12)    
        nDimensions = 2;
        data = new double[]{0, 0, 0, 0, 0, 12};
        expected[0] = new double[]{0, 4};
        expected[1] = new double[]{0, 6.928};    
        
        c = MiscMath0.standardDeviation(data, nDimensions);
        assertEquals(expected.length, c.length);
        nc = 0;
        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[i].length; ++j) {
                if (Math.abs(expected[i][j] - c[i][j]) < tol) {
                    nc++;
                }
            }
        }
        assertEquals(expected.length*expected[0].length, nc);
        
        //(-20, 48), (-20, -48), (20, 0), (59, 0)
        nDimensions = 2;
        data = new double[]{-20, 48, -20, -48, 20, 0, 59, 0};
        expected[0] = new double[]{9.75, 0};
        expected[1] = new double[]{37.863, 39.192};
                
        c = MiscMath0.standardDeviation(data, nDimensions);
        assertEquals(expected.length, c.length);
        nc = 0;
        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[i].length; ++j) {
                if (Math.abs(expected[i][j] - c[i][j]) < tol) {
                    nc++;
                }
            }
        }
        assertEquals(expected.length*expected[0].length, nc);
    }
    
    public static void testCumulativeSum() {
        
        double[] a = new double[]{2.1, 3.2, 9, 12.1};
        
        double[] eS = new double[]{2.1, 5.3, 14.3, 26.4};
        
        double[] s = MiscMath0.cumulativeSum(a);
        
        assertEquals(eS.length, s.length);
        
        double eps = 1.e-17;
        
        double diff;
        for (int i = 0; i < a.length; ++i) {
            diff = eS[i] - s[i];
            assertTrue(diff < eps);
        }
        
    }
    
    public void testFrequencyMap() {
        
        int[] a = new int[]{ 3, 4, 2, 1, 3, 9, 23, 2, 3};
        
        TIntIntMap expected = new TIntIntHashMap();
        expected.put(3, 3);
        expected.put(4, 1);
        expected.put(2, 2);
        expected.put(1, 1);
        expected.put(9, 1);
        expected.put(23, 1);
        
        TIntIntMap f = MiscMath0.makeFrequencyMap(a);
        
        assertEquals(expected.size(), f.size());
        
        int key, c, i;
        TIntIntIterator iter = expected.iterator();
        for (i = 0; i < expected.size(); ++i) {
            iter.advance();
            key = iter.key();
            c = iter.value();
            assertTrue(f.containsKey(key));
            assertEquals(c, f.get(key));
        }
    }
    
    public void test() {
        for (int i = 0; i < 31; ++i) {
            int v = 1 << i;
            String bs = Integer.toBinaryString(v);
            int nb = MiscMath0.numberOfBits(v);
            
            //System.out.println("i=" + i + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
        for (int i = 0; i < 31; ++i) {
            int v = 1 << i;
            String bs = Integer.toBinaryString(v);
            v *= -1;
            int nb = MiscMath0.numberOfBits(v);
            
            //System.out.println("i=" + i + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
    }
    
    public void test1() {
        for (int i = 0; i < 31; ++i) {
            int v = 1 << i;
            String bs = Integer.toBinaryString(v);
            int nb = MiscMath0.numberOfBitsWOB(v);
            
            //System.out.println("i=" + i + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
        for (int i = 0; i < 31; ++i) {
            int v = 1 << i;
            String bs = Integer.toBinaryString(v);
            v *= -1;
            int nb = MiscMath0.numberOfBitsWOB(v);
            
            //System.out.println("bs=" + bs + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
    }
    
    public void test2() {
        for (int i = 1; i < 63; ++i) {
            long v = 1L << i;
            int nb = MiscMath0.numberOfBits(v);
            String bs = Long.toBinaryString(v);
            //System.out.println("i=" + i + " bs=" + bs + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
    }

    public void test3() {
        for (int i = 0; i < 31; ++i) {
            int nb = MiscMath0.numberOfBits(i);
            String bs = Integer.toBinaryString(i);
            //System.out.println("v=" + bs + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
    }
    public void test4() {
        for (int i = 0; i < 31; ++i) {
            int nb = MiscMath0.numberOfBits((long)i);
            String bs = Integer.toBinaryString(i);
            //System.out.println("v=" + bs + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
    }
    
    public void testComputeNDivNMinusK() throws Exception {

        // n! / (n - k)!
        long result = MiscMath0.computeNDivNMinusK(2, 1);
        assertTrue(result == 2);

        result = MiscMath0.computeNDivNMinusK(4, 3);
        assertTrue(result == 4*3*2);

        result = MiscMath0.computeNDivNMinusK(4, 2);
        assertTrue(result == 4*3);

        result = MiscMath0.computeNDivNMinusK(6, 3);
        assertTrue(result == 6*5*4);
    }

    public void testFactorial() throws Exception {

        long result = MiscMath0.factorial(2);
        assertTrue(result == 2);

        result = MiscMath0.factorial(4);
        assertTrue(result == 4*3*2);

        result = MiscMath0.factorial(5);
        assertTrue(result == 5*4*3*2);

        result = MiscMath0.factorial(1);
        assertTrue(result == 1);

        result = MiscMath0.factorial(0);
        assertTrue(result == 0);
    }
    
}
