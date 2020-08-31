package algorithms.misc;

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
    
    public void est0() {
        
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
        double srch;
        int idx;
        int[] selected;
        
        double tol = 3.e-5;
        int k = 3;
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        System.out.println("seed=" + seed);
        System.out.flush();
        rand.setSeed(seed);
        
        cdf = new double[99];
        for (int i = 0; i < 99; ++i) {
            cdf[i] = i + 1;
        }
        // normalize so that the last bin is "1".
        double norm = cdf[cdf.length - 1];
        for (int i = 0; i < cdf.length; ++i) {
            cdf[i] /= norm;
        }
        
        k = 5 * cdf.length;
        
        selected = CDFRandomSelect.chooseKFromBinarySearch(cdf, k, rand);
        
        assertEquals(k, selected.length);
        
        TIntIntMap freq = new TIntIntHashMap();
        int key, val;
        for (int i = 1; i < k; ++i) {
            key = selected[i];
            if (freq.containsKey(key)) {
                val = freq.get(key) + 1;
            } else {
                val = 1;
            }
            freq.put(key, val);
        }
        // fits a uniform distribution ?
        
        //int[] chooseKFromBinarySearch(double[] cdf, int k, Random rand, double tolerance) {
        //    public static int[] chooseKFromIntegerTransformAndTrie(double[] cdf, int k, Random rand) {       
    }

}
