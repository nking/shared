package algorithms.misc;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.logging.Logger;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
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
    }
}
