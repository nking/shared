package algorithms.matrix;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CURDecompositionTest extends TestCase {
    
    public CURDecompositionTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
        // from Chap 11 of book "Mining of Massive 
        // Datasets" by Jure Leskovec, Anand Rajaraman, Jeff Ullman
        // http://www.mmds.org/
        
        double[][] a = new double[7][5];
        a[0] = new double[]{1, 1, 1, 0, 0};
        a[1] = new double[]{3, 3, 3, 0, 0};
        a[2] = new double[]{4, 4, 4, 0, 0};
        a[3] = new double[]{5, 5, 5, 0, 0};
        a[4] = new double[]{0, 0, 0, 4, 4};
        a[5] = new double[]{0, 0, 0, 5, 5};
        a[6] = new double[]{0, 0, 0, 2, 2};
        
        double k = 2;
        
        
        
        
    }
}
