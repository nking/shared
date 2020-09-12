package algorithms.misc;

import algorithms.misc.ROC.ROCResults;
import algorithms.util.PairFloatArray;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ROCTest extends TestCase {
    
    public ROCTest(String testName) {
        super(testName);
    }

    
    public void testCalcAUCAndPoints() {
        double[] fScores = new double[]
        {.1, .3, .33, .34, .35, .36, .37, .38, .39, .4,
        .9, .8, .7, .6, .55, .54, .53, .52, .51, .505};
        boolean[] labels = new boolean[]
        {false, true, false, true, false, false, false, true, false, true,
        true, true, false, true, true, true, false, false, true, false};
        
        assertEquals(fScores.length, labels.length);
        
        ROCResults rr = ROC.calcAUCAndPoints(fScores, labels);
        
        assertEquals(labels.length, rr.sortedFScores.length);
        assertEquals(labels.length, rr.indexes.length);
        assertEquals(labels.length, rr.pts.getN() - 1);
        
        // test for descending sort
        int idx;
        double diff;
        double tol = 1.e-15;
        double prev = Double.POSITIVE_INFINITY;
        for (int i = 0; i < rr.indexes.length; ++i) {
            assertTrue(rr.sortedFScores[i] <= prev);
            idx = rr.indexes[i];
            diff = Math.abs(rr.sortedFScores[i] - fScores[idx]);
            assertTrue(diff < tol);
            prev = rr.sortedFScores[i];
        }
        
        PairFloatArray epts = new PairFloatArray(21);
        epts.add(0, 0);
        epts.add(0, .1f);
        epts.add(0, .2f);
        
        epts.add(0.1f, 0.2f);
        epts.add(0.1f, 0.3f);
        epts.add(0.1f, 0.4f);
        epts.add(0.1f, 0.5f);
        
        epts.add(0.2f, 0.5f);
        
        epts.add(0.3f, 0.5f);
        epts.add(0.3f, 0.6f);
        
        epts.add(0.4f, 0.6f);
        epts.add(0.4f, 0.7f);
        
        epts.add(0.5f, 0.7f);
        epts.add(0.5f, 0.8f);
        
        epts.add(0.6f, 0.8f);
        
        epts.add(0.7f, 0.8f);
        
        epts.add(0.8f, 0.8f);
        epts.add(0.8f, 0.9f);
        
        epts.add(0.9f, 0.9f);
        epts.add(0.9f, 1.0f);
        
        epts.add(1.0f, 1.0f);
        
        assertEquals(epts.getN(), rr.pts.getN());
        
        for (int i = 0; i < epts.getN(); ++i) {
            diff = Math.abs(rr.pts.getX(i) - epts.getX(i));
            assertTrue(diff < tol);
            diff = Math.abs(rr.pts.getY(i) - epts.getY(i));
            assertTrue(diff < tol);
        }
    }
    
}
