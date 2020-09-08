package algorithms.misc;

import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;
import java.util.Random;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CDFRandomSelectTest extends TestCase {
    
    public CDFRandomSelectTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
        double[] cdf;
        double srch;
        int idx;
        
        double tol = 3.e-5;
        
        cdf = new double[]{ 1., 2., 3., 3.1, 4.};
        srch = 3.04;
        idx = CDFRandomSelect.binarySearchForNearest(cdf, srch, tol);
        assertEquals(2, idx);
        
        srch = 1.04;
        idx = CDFRandomSelect.binarySearchForNearest(cdf, srch, tol);
        assertEquals(0, idx);
        
        srch = 0.04;
        idx = CDFRandomSelect.binarySearchForNearest(cdf, srch, tol);
        assertEquals(0, idx);
        
        srch = 3.15;
        idx = CDFRandomSelect.binarySearchForNearest(cdf, srch, tol);
        assertEquals(3, idx);
        
        srch = 3.6;
        idx = CDFRandomSelect.binarySearchForNearest(cdf, srch, tol);
        assertEquals(4, idx);
        
        srch = 5.0;
        idx = CDFRandomSelect.binarySearchForNearest(cdf, srch, tol);
        assertEquals(4, idx);
        
        for (int i = 0; i < cdf.length; ++i) {
            srch = cdf[i] + tol/2.;
            idx = CDFRandomSelect.binarySearchForNearest(cdf, srch, tol);
            assertEquals(i, idx);
        }
    }
    
    public void test1() {
        
        double[] cdf;
        int[] selected;
                
        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        //seed = 180328550254112L;
        System.out.println("seed=" + seed);
        System.out.flush();
        rand.setSeed(seed);
        
        cdf = new double[100];
        for (int i = 0; i < cdf.length; ++i) {
            cdf[i] = i + 1;
        }
        // normalize so that the last bin is "1".
        double norm = cdf[cdf.length - 1];
        for (int i = 0; i < cdf.length; ++i) {
            cdf[i] /= norm;
        }
        
        int k = 1 * cdf.length;
        
        selected = CDFRandomSelect.chooseKFromBinarySearch(cdf, k, rand);
        
        assertEquals(k, selected.length);
        
        TIntIntMap freq = new TIntIntHashMap();
        int key;
        int val;
        int minKey = Integer.MAX_VALUE;
        int maxKey = Integer.MIN_VALUE;
        for (int i = 0; i < k; ++i) {
            key = selected[i];
            if (freq.containsKey(key)) {
                val = freq.get(key) + 1;
            } else {
                val = 1;
            }
            freq.put(key, val);
            if (key < minKey) {
                minKey = key;
            }
            if (key > maxKey) {
                maxKey = key;
            }
        }
                 
        double nBins = cdf.length;
        double eUniform = (double)k/nBins;//nDraws/(double)nBins;
        
        int dOF = cdf.length - 1 - 1;
        double chiSqStat = ChiSquaredCriticalValues.upperTailStatistic(
            ChiSquaredCriticalValues.PROB_UT.Z_90, dOF);
                
        double expectedSigmaSq = (nBins*nBins)/12;
        double expectedSigma = Math.sqrt(expectedSigmaSq);
        
        double[] meanAndStDev = MiscMath0.getAvgAndStDev(selected);
        assertTrue(Math.abs(meanAndStDev[1] - expectedSigma) < (expectedSigma/5.));
        
        double diff;
        double statSq = 0;
        double frequency;
        int totVals = 0;
        //for (int key2 : sortedKeys) {
        for (int i = 0; i < cdf.length; ++i) {
            if (freq.containsKey(i)) {
                frequency = freq.get(i);
            } else {
                frequency = 0;
            }
            totVals += frequency;
            diff = frequency - eUniform;
            diff *= diff;
            statSq += (diff/eUniform);
        }
        
        double pVal = ChiSquaredCriticalValues.approxPValueLin(statSq, dOF);
        double chiSqStat3 = ChiSquaredCriticalValues.approxChiSqStatLin(
            1./(double)nBins, dOF);
        double chiSqStat4 = ChiSquaredCriticalValues.approxChiSqStatLin(
            1./(double)k, dOF);
        
        // if statSq < chiSqStat, there is no evidence to suggest the
        //    distribution is not uniform.
        //assertTrue(statSq < chiSqStat);
        
    }
    
    public void test2() {
        
        double[] cdf;
        int[] selected;
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        seed = 182153812288260L;
        System.out.println("*seed=" + seed);
        System.out.flush();
        rand.setSeed(seed);
        
        cdf = new double[100];
        for (int i = 0; i < cdf.length; ++i) {
            cdf[i] = i + 1;
        }
        // normalize so that the last bin is "1".
        double norm = cdf[cdf.length - 1];
        for (int i = 0; i < cdf.length; ++i) {
            cdf[i] /= norm;
        }
        
        int k = 5 * cdf.length;
        
        selected = CDFRandomSelect.chooseKFromBinarySearch(cdf, k, rand);
        
        assertEquals(k, selected.length);
        
        TIntIntMap freq = new TIntIntHashMap();
        int key;
        int val;
        int minKey = Integer.MAX_VALUE;
        int maxKey = Integer.MIN_VALUE;
        for (int i = 0; i < k; ++i) {
            key = selected[i];
            if (freq.containsKey(key)) {
                val = freq.get(key) + 1;
            } else {
                val = 1;
            }
            freq.put(key, val);
            if (key < minKey) {
                minKey = key;
            }
            if (key > maxKey) {
                maxKey = key;
            }
        }
                 
        double nBins = cdf.length;
        double eUniform = (double)k/nBins;//nDraws/(double)nBins;
        
        int dOF = cdf.length - 1;
        double chiSqStat = ChiSquaredCriticalValues.upperTailStatistic(
            ChiSquaredCriticalValues.PROB_UT.Z_95, dOF);
        
        double expectedSigmaSq = (nBins*nBins)/12;
        double expectedSigma = Math.sqrt(expectedSigmaSq);
        
        double[] meanAndStDev = MiscMath0.getAvgAndStDev(selected);
        assertTrue(Math.abs(meanAndStDev[1] - expectedSigma) < (expectedSigma/5.));
        
        double diff;
        double statSq = 0;
        double frequency;
        int totVals = 0;
        //for (int key2 : sortedKeys) {
        for (int i = 0; i < cdf.length; ++i) {
            if (freq.containsKey(i)) {
                frequency = freq.get(i);
            } else {
                frequency = 0;
            }
            totVals += frequency;
            diff = frequency - eUniform;
            diff *= diff;
            statSq += (diff/eUniform);
        }
        
        double pVal = ChiSquaredCriticalValues.approxPValueLin(statSq, dOF);
        double chiSqStat3 = ChiSquaredCriticalValues.approxChiSqStatLin(
            1./nBins, dOF);
        double chiSqStat4 = ChiSquaredCriticalValues.approxChiSqStatLin(
            1./(double)k, dOF);
        
        // if statSq < chiSqStat, there is no evidence to suggest the
        //    distribution is not uniform.
        assertTrue(statSq < chiSqStat || statSq < chiSqStat4);
        
    }

}
