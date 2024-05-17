package algorithms.statistics;

import algorithms.misc.Misc0;
import algorithms.misc.MiscMath0;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
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

        double _n=4.0;
        cdf = new double[]{ 1./_n, 2./_n, 3./_n, 3./_n, 3./_n, 3.1/_n, 4./_n};
        srch = 3.0/_n;
        idx = CDFRandomSelect.binarySearchForQuantile(cdf, srch);
        assertEquals(2, idx);

        cdf = new double[]{ 1./4., 2./4., 3./4., 3.1/4., 4./4.};
        srch = 3.04/4.;
        idx = CDFRandomSelect.binarySearchForQuantile(cdf, srch);
        assertEquals(3, idx);

        srch = 3.0/4.;
        idx = CDFRandomSelect.binarySearchForQuantile(cdf, srch);
        assertEquals(2, idx);
        
        srch = 1.0/4.;
        idx = CDFRandomSelect.binarySearchForQuantile(cdf, srch);
        assertEquals(0, idx);

        srch = 1.04/4.;
        idx = CDFRandomSelect.binarySearchForQuantile(cdf, srch);
        assertEquals(1, idx);
        
        srch = 0.04/4.;
        idx = CDFRandomSelect.binarySearchForQuantile(cdf, srch);
        assertEquals(0, idx);
        
        srch = 3.15/4.;
        idx = CDFRandomSelect.binarySearchForQuantile(cdf, srch);
        assertEquals(4, idx);
        
        srch = 3.6/4.;
        idx = CDFRandomSelect.binarySearchForQuantile(cdf, srch);
        assertEquals(4, idx);
        
        srch = 5.0/4.;
        idx = CDFRandomSelect.binarySearchForQuantile(cdf, srch);
        assertEquals(4, idx);
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
        
        /*
        goodness of fit tests:
            the test only rejects the assumed distribution when there is 
                definite evidence that the distribution is incorrect.
            in hypothesis testing, 
                type I error = rejecting a true null hypothesis
                type II error = accepting a false null hypothesis.
            The power of a goodness-of-fit is the probability that the test
                will rcjecL the null hypothesis.
        */
        
        double pVal = ChiSquaredCriticalValues.approxPValueLin(statSq, dOF);
        double chiSqStat3 = ChiSquaredCriticalValues.approxChiSqStatLin(
            1./nBins, dOF);
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
        //seed = 182153812288260L;
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
        
        System.out.printf("statSq=%.4e\nchiSqStat=%.4e\nchiSqStat3 = %.4e\nchiSqStat4 = %.4e\n ",
            statSq, chiSqStat,chiSqStat3, chiSqStat4);
        System.out.flush();
        
        // if statSq < chiSqStat, there is no evidence to suggest the
        //    distribution is not uniform.
        assertTrue(statSq < chiSqStat || statSq < chiSqStat4);
        
    }

}
