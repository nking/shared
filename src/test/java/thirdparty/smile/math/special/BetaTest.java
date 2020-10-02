package thirdparty.smile.math.special;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BetaTest extends TestCase {
    
    public BetaTest(String testName) {
        super(testName);
    }

    /**
     * Test of beta method, of class Beta.
     */
    public void testRegularizedIncompleteBeta() {
        
// from https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainc.html#scipy.special.betainc
        
        double tol = 1e-7;
        
        double a = 1.4;
        double b = 3.1;
        double x = 0.5;
        
        double p = Beta.regularizedIncompleteBetaFunction(a, b, x);
                
        assertTrue(Math.abs(p - 0.8148904036225296) < tol);
        
        a = 2.2;
        b = 3.1;
        x = 0.4;
        
        p = Beta.regularizedIncompleteBetaFunction(a, b, x);
                
        assertTrue(Math.abs(p - 0.49339638807619446) < tol);
        
        a = 2.2;
        b = 3.1;
        x = 0.4;
        
        p = Beta.regularizedIncompleteBetaFunction(b, a, 1. - x);
        double p2 = 1. - p;
                
        assertTrue(Math.abs(p2 - 0.49339638807619446) < tol);
    }

    /**
     * Test of inverseRegularizedIncompleteBetaFunction method, of class Beta.
     */
    public void testInverseRegularizedIncompleteBetaFunction() {
        // from https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betaincinv.html#scipy.special.betaincinv
        
        double tol = 1e-15;
        
        double a = 1.2;
        double b = 3.1;
        double x = 0.2;
        
        double p = Beta.regularizedIncompleteBetaFunction(a, b, x);
        
        double x2 = Beta.inverseRegularizedIncompleteBetaFunction(a, b, p);
        
        assertTrue(Math.abs(x2 - x) < tol);
        
        a = 7.5;
        b = 0.4;
        x = 0.5;
        p = Beta.regularizedIncompleteBetaFunction(a, b, x);
        x2 = Beta.inverseRegularizedIncompleteBetaFunction(a, b, p);     
        
        assertTrue(Math.abs(x2 - x) < tol);
    }
    
}
