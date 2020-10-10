package algorithms.misc;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DistancesTest extends TestCase {
    
    public DistancesTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
        double tol = 1e-5;
        
        double[] p1, p2;
        
        p1 = new double[]{Math.sqrt(2), Math.sqrt(2)};
        p2 = new double[]{Math.sqrt(2), 0};
        
        double dsq = Distances.calcEuclideanSquared(p1, p2);
        
        assertTrue(Math.abs(2. - dsq) < tol);
    }
}
