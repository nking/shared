package algorithms.misc;

import algorithms.signalProcessing.Util;
import algorithms.statistics.UnivariateNormalDistribution;
import algorithms.util.FormatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;

import java.security.NoSuchAlgorithmException;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.math.BigInteger;

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
    
    public void testGetAvgAndStDev_long() throws Exception {
        
        long[] x = new long[]{2, 2, 2, 4, 4, 4};
        
        double[] avgAndStDev = MiscMath0.getAvgAndStDev(x, x.length);
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

    public static void testNumberOfSetBits() {
        //    public static int numberOfSetBits(long x) {
        Random rand = new Random();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);

        int nTests = 100;
        int nb;
        long x;
        int ns;
        int r;
        int j;
        TIntSet set = new TIntHashSet();
        for (int i = 0; i < nTests; ++i) {
            x = 0;
            set.clear();
            // number of bits to set in x.  some may be the same bits
            nb = rand.nextInt(62);
            for (j = 0; j < nb; ++j) {
                r = rand.nextInt(62);
                set.add(r);
                x |= (1L << r);
            }
            ns = MiscMath0.numberOfSetBits(x);
            assertEquals(set.size(), ns);
        }
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

    public void testMSBRandom() {
        long seed = System.nanoTime();
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);
        int nTests = 100;
        for (int i = 0; i < nTests; ++i) {
            long r = rand.nextLong();
            int expAns = 63 - Long.numberOfLeadingZeros(Math.abs(r));
            int expAns2 = MiscMath0.numberOfBits(r) - 1;
            assert(expAns2 == expAns);
            int ans = MiscMath0.MSBWithoutBuiltIn(r);
            assertEquals(expAns, ans);
        }
    }
    public void testLSBRandom() {
        long seed = System.nanoTime();
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);
        int nTests = 100;
        for (int i = 0; i < nTests; ++i) {
            long r = rand.nextLong();
            String str = Long.toBinaryString(r);
            int expAns = 0;
            for (int j = str.length() - 1; j >= 0; --j) {
                if (str.charAt(j) == '1') {
                    break;
                }
                ++expAns;
            }
            assertEquals(expAns, MiscMath0.LSB(r));
            assertEquals(expAns, MiscMath0.LSBWithoutBuiltIn1(r));
        }
    }

    public void testMSB_LSB() {
        long t0 = (1L<<4)+1; // 17, 0b10001,  MSB=4, nBits=5
        assertEquals(4, MiscMath0.MSBWithoutBuiltIn(t0));
        assertEquals(5, MiscMath0.numberOfBits(t0));
        assertEquals(4, MiscMath0.MSBWithoutBuiltIn(-1*t0));
        assertEquals(5, MiscMath0.numberOfBits(-1*t0));

        assertEquals(62, MiscMath0.MSBWithoutBuiltIn(Long.MAX_VALUE));
        assertEquals(62, MiscMath0.MSBWithoutBuiltIn(Long.MIN_VALUE));

        assertEquals(62, MiscMath0.MSBWithoutBuiltIn(Long.MAX_VALUE-1));
        assertEquals(62, MiscMath0.MSBWithoutBuiltIn(Long.MIN_VALUE+1));

        assertEquals(61, MiscMath0.MSBWithoutBuiltIn(1L<<62));

        long t1 = 0b10010000;
        assertEquals(4, MiscMath0.LSB(t1));
        assertEquals(4, MiscMath0.LSBWithoutBuiltIn1(t1));

        assertEquals(4, MiscMath0.LSB(-1*t1));
        assertEquals(4, MiscMath0.LSBWithoutBuiltIn1(-1*t1));
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
    
    public void testcomputeNDivKTimesNMinusK() throws Exception {
        
        int n, k;
        long nComb;
        
        n = 6;
        k = 2;
        nComb = MiscMath0.computeNDivKTimesNMinusK(n, k);
        assertEquals(15, nComb);
        
        n = 2;
        k = 1;
        nComb = MiscMath0.computeNDivKTimesNMinusK(n, k);
        assertEquals(2, nComb);
        
        n = (int)Math.sqrt(Integer.MAX_VALUE);
        k = 2;
        nComb = MiscMath0.computeNDivKTimesNMinusK(n, k);
        assertEquals(((n*(n-1))/2), nComb);
        
        n = (int)Math.sqrt(Integer.MAX_VALUE);
        k = 2;
        nComb = MiscMath0.computeNDivKTimesNMinusK(n, k);
        assertEquals(((n*(n-1))/2), nComb);
        
        n = 20;
        k = 4;
        long expected = 4845;
        nComb = MiscMath0.computeNDivKTimesNMinusK(n, k);
        assertEquals(expected, nComb);
        
        n = 100;
        k = 90;
        expected = 17310309456440L;
        nComb = MiscMath0.computeNDivKTimesNMinusK(n, k);
        assertEquals(expected, nComb);
        
        n = 100;
        k = 80;
        BigInteger expected2 = new BigInteger("535983370403809682970");
        BigInteger nComb2 = MiscMath0.computeNDivKTimesNMinusKBigInteger(n,k);
        assertEquals(expected2, nComb2);
        
        n = 100;
        k = 20;
        expected2 = new BigInteger("535983370403809682970");
        nComb2 = MiscMath0.computeNDivKTimesNMinusKBigInteger(n,k);
        assertEquals(expected2, nComb2);
    }
    
    public void testFactorialBigInteger() throws Exception {

        BigInteger expected = new BigInteger("87178291200");
        
        BigInteger result = MiscMath0.factorialBigInteger(14);

        assertTrue(result.compareTo(expected) == 0);
    }
    
    public void testMAD() {
        // test from https://en.m.wikipedia.org/wiki/Median_absolute_deviation
        double[] x = new double[]{1, 1, 2, 2, 4, 6, 9};
        double[] mad = MiscMath0.calculateMedianOfAbsoluteDeviation(x);
        double tol = 1e-3;
        double expectedMAD = 1;
        double expectedMedian = 2;
        double expectedMin = 1;
        double expectedMax = 9;
        assertTrue(Math.abs(mad[0] - expectedMAD) < tol);
        assertTrue(Math.abs(mad[1] - expectedMedian) < tol);
        assertTrue(Math.abs(mad[2] - expectedMin) < tol);
        assertTrue(Math.abs(mad[3] - expectedMax) < tol);
    }
    
    public void testQuartiles() {
        // test is from https://en.m.wikipedia.org/wiki/Interquartile_range
        double[] x = new double[]{
            7, 7, 31, 31, 47, 75, 87, 115, 116, 119, 119, 155, 177
        };
        
        double[] expected = new double[]{
            31, 87, 119, 7, 177
        };
        // returns q1, q2, q3, min, max
        double[] q = MiscMath0.calculateQuartiles(x);
        
        assertEquals(expected.length, q.length);
        
        double tol = 1e-16;
        for (int i = 0; i < q.length; ++i) {
            assertTrue(Math.abs(expected[i] - q[i]) < tol);
        }
        
        // test tukey while at it
        // iqr= 119 - 31 = 88.
        // r0 = 31 - k*88 = 
        // r1 = 119 + k*88;
        
        int[] expectedInliers;
        int n = x.length;
        
        double k = 0.2; // very low value just for tests
        expectedInliers = new int[]{2, 3, 4, 5, 6, 7, 8, 9, 10};
        
        int[] inliers = MiscMath0.findInliersUsingTukeyFences(x, k);
        
        assertEquals(expectedInliers.length, inliers.length);
        
        for (int i = 0; i < q.length; ++i) {
            assertTrue(Math.abs(expectedInliers[i] - inliers[i]) < tol);
        }
    }

    public void testAreColinear() {

        //note that the z-axis is ignored.

        double s0 = 2.;
        double s2 = 3.;
        // these columns are perpendicular to one another, so not colinear
        double[][] x = new double[3][3];
        x[0] = new double[]{s0*1, 0, 0};
        x[1] = new double[]{0, 1, 0};
        x[2] = new double[]{0, 0, s2*1};

        boolean a = MiscMath0.areColinear(x, 1E-6);
        assertFalse(a);

        // make 2 columns parallel.
        x = new double[3][4];
        x[0] = new double[]{s0*1, 0,    0, 0};
        x[1] = new double[]{0,    1,    0, 0};
        x[2] = new double[]{0,    0, s2*1, 1};

        a = MiscMath0.areColinear(x, 1E-6);
        assertTrue(a);

        // make 3 columns parallel.  they should be colinear
        x = new double[3][4];
        x[0] = new double[]{s0*1, 0,    0, 0, 0};
        x[1] = new double[]{0,    1,    0, s0*1, s2*1};
        x[2] = new double[]{0,    0, s2*1, 0, 0};
        a = MiscMath0.areColinear(x, 1E-6);
        assertTrue(a);

        /*
        -5.144e-02, 1.089e+00, 1.685e+00, 1.235e+00, 9.944e-01, 1.311e+00, 7.789e-01
        1.026e+00, 1.026e+00, 8.093e-01, -6.703e-01, -7.966e-01, -2.372e-01, 3.402e-01
        1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00

        1.235e+00, 9.944e-01, 1.311e+00, 1.102e+00, 9.437e-01, 7.789e-01, 7.092e-01
        -6.703e-01, -7.966e-01, -2.372e-01, 8.454e-01, -4.718e-01, 3.402e-01, -1.230e+00
        1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00

        8.930e-01, -5.144e-02, 1.089e+00, 1.235e+00, 1.311e+00, 7.789e-01, 7.092e-01
        7.372e-01, 1.026e+00, 1.026e+00, -6.703e-01, -8.327e-01, 3.402e-01, -1.230e+00
        1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00
         */
    }

    public void testEqualWithinTolerance() {
        double[] a = new double[]{1,2,3,4,5};
        assertTrue(MiscMath0.equalsWithinTolerance(a, a, 1E-11));
        double[] b = Arrays.copyOf(a, a.length);
        b[3] += 1E-4;
        assertFalse(MiscMath0.equalsWithinTolerance(a, b, 1E-5));
    }

    public void testUniformRandomDoubleTest() {

        // test that using same seed produces same Unif(0,1) numbers and that the numbers fit a uniform distribution
        long seed = 1234567L;
        Random rand = new Random(seed);

        // Unif(0,1); E[x]=0.5, Var(x)=0.25
        int n = 1000;
        double[] a1 = MiscMath0.uniformRandomDouble(n, rand);
        double[] s1 = MiscMath0.getAvgAndStDev(a1);
        assertTrue(Math.abs(s1[0] - 0.5) < 0.01);
        assertTrue(Math.abs(s1[1] - 0.25) < 0.1);

        // check that seeded rand produces repeatable numbers for use in tests
        rand = new Random(seed);
        double[] a2 = MiscMath0.uniformRandomDouble(n, rand);

        assertTrue(MiscMath0.equalsWithinTolerance(a1, a2, 1E-11));
    }

    public void testLog() {
        // simple test that each number is log of x
        double[] x = MiscMath0.uniformRandomDouble(10);
        double[] logX = MiscMath0.log(x);
        for (int i = 0; i < x.length; ++i) {
            assertTrue(Math.abs(x[i] - Math.exp(logX[i])) < 1E-11);
        }
    }

    public void testCalcMean() {
        int n = 10000;
        double[] x = MiscMath0.uniformRandomDouble(n);
        assertTrue(Math.abs(0.5 - MiscMath0.mean(x)) < 1E-2);
    }

    public void testCalcMeanAndSSD_double() throws NoSuchAlgorithmException {
        int n = 100000;
        double[] x = UnivariateNormalDistribution.randomSampleOfUnitStandard(n);
        double[] ms = MiscMath0.calcMeanAndSSD(x);
        double stDev = Math.sqrt(ms[1]/(n-1));
        assertTrue(Math.abs(ms[0]) < 1E-2);//0
        assertTrue(Math.abs(stDev - 1) < 1E-2);//1
    }

    public void testCalcMeanAndSSD_int() throws Exception {
        int n = 100000;
        double[] x = UnivariateNormalDistribution.randomSampleOfUnitStandard(n);
        int[] xInt = new int[n];
        int f = 1;
        for (int i = 0; i < n; ++i) {
            xInt[i] = (int)Math.round(f * x[i]);
        }
        double[] ms = MiscMath0.calcMeanAndSSD(xInt);
        double stDev = Math.sqrt(ms[1]/(n-1));
        assertTrue(Math.abs(ms[0]) < 1E-1);//0
        assertTrue(Math.abs(stDev - 1) < 1E-1);//1
    }

    public void testCalcGeometricMean() {
        double[] x = new double[]{1,3,7,5,3,11};
        double exAns = 3.8900084061965643;
        //(product(x_i))^(1/ n)
        double gm = MiscMath0.calcGeometricMean(x);
        assertTrue(Math.abs(exAns - gm) < 1E-11);

    }

    public void testCalcHarmonicMean() {
        double[] x = new double[]{1,3,7,5,3,11};
        double exAns = 2.8565539983511954;
        //n/(sum(1/ a_i))
        double hm = MiscMath0.calcHarmonicMean(x);
        assertTrue(Math.abs(exAns - hm) < 1E-11);
    }

    public void testMultiDimensionToSingleIndex() {
        int[] indexes, dims, out;
        int idx;

        indexes = new int[]{1, 4};
        dims = new int[]{10, 11};
        out = new int[2];
        idx = MiscMath0.multiDimensionToSingleIndex(indexes, dims);
        MiscMath0.singleIndexToMultiDimension(idx, dims, out);
        assertTrue(Arrays.equals(indexes, out));

        long seed = System.nanoTime();
        //seed = 23539291542956L;
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);

        int nTests = 100;
        int maxDim = 5;
        // maxVal^(maxDim) <= Integer.MAX_VALUE
        int maxIdx = (int)Math.pow(Integer.MAX_VALUE, 1./maxDim);

        int nDim, maxWidth;
        for (int i = 0; i < nTests; ++i) {
            nDim = 1 + rand.nextInt(maxDim - 1);
            maxWidth = 2 + rand.nextInt(maxIdx - 2);

            indexes = new int[nDim];
            dims = new int[nDim];
            out = new int[nDim];

            for (int k = 0; k < nDim; ++k) {
                dims[k] = 1 + rand.nextInt(maxWidth - 1);
                indexes[k] = rand.nextInt(dims[k]);
            }

            idx = MiscMath0.multiDimensionToSingleIndex(indexes, dims);
            MiscMath0.singleIndexToMultiDimension(idx, dims, out);
            assertTrue(Arrays.equals(indexes, out));

        }
    }

    public void testCalcQuartilesMinMax() {
        double[] x, qM, expAns;

        x = new double[]{6, 7, 15, 36, 39, 40, 41, 42, 43, 47, 49};
        expAns = new double[]{6, 15., 40., 43., 49.};
        qM = MiscMath0.calcQuartiles(x);
        assertEquals(expAns.length, qM.length);
        for (int i = 0; i < expAns.length; ++i) {
            assertTrue(Math.abs(expAns[i] - qM[i]) < 1E-7);
        }

        x = new double[]{7, 15, 36, 39, 40, 41.};
        expAns = new double[]{7, 15., 37.5, 40., 41.};
        qM = MiscMath0.calcQuartiles(x);
        assertEquals(expAns.length, qM.length);
        for (int i = 0; i < expAns.length; ++i) {
            assertTrue(Math.abs(expAns[i] - qM[i]) < 1E-7);
        }
    }
}
