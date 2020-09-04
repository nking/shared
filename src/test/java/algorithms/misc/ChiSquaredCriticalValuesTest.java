package algorithms.misc;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ChiSquaredCriticalValuesTest extends TestCase {
    
    public ChiSquaredCriticalValuesTest(String testName) {
        super(testName);
    }
    
    public void testLin() {
        
        double p;
        int dOF;
        
        double[] x2 = new double[]{2.706, 3.841, 5.024, 6.635, 10.828};
        dOF = 1;
                
        for (double chiSquared : x2) {
            
            p = ChiSquaredCriticalValues.approxUpperTailPValueLin(chiSquared, dOF);
            
        }
                  
    }
    
    public void test0() {
                
        double p;
        int dOF;
        
        double[] x2 = new double[]{2.706, 3.841, 5.024, 6.635, 10.828};
        dOF = 1;
                
        for (double chiSquared : x2) {
            
            p = ChiSquaredCriticalValues.approxUpperTailPValueBeh(chiSquared, dOF);
            
        }
    }
}
