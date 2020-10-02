package algorithms.statistics;

import algorithms.misc.Misc0;
import algorithms.misc.MiscMath0;
import java.util.Random;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CDFStandardNormalTest extends TestCase {
    
    public CDFStandardNormalTest(String testName) {
        super(testName);
    }
    
    public void testApproxShore() {
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        //long seed = 1595442111027L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        int n = 1000;
        double[] g = new double[n];
        
        for (int i = 0; i < n; ++i) {
            g[i] = CDFStandardNormal.approxInverseShort(rand.nextDouble());
        }
        
        double tol = 0.1;
        double expectedSigma = 1;
        double expectedMean = 0;
        
        double[] meanAndStDev = MiscMath0.getAvgAndStDev(g);
        assertTrue(Math.abs(meanAndStDev[0] - expectedMean) < tol);
        assertTrue(Math.abs(meanAndStDev[1] - expectedSigma) < tol);
        
    }
}
