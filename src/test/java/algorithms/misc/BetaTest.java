/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package algorithms.misc;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BetaTest extends TestCase {
    
    public BetaTest(String testName) {
        super(testName);
    }
 
    public void test0() {
        
        double eps = 1.e-15;
        double tol = 1e-9;
        double a, b, r;
        
        a = 0.5; b = 0.5;
        r = Beta.beta(a, b);
        assertTrue(Math.abs(r - Math.PI) < eps);
        
        a = 1.; b = 1;
        r = Beta.beta(a, b);
        assertTrue(Math.abs(r - 1.) < eps);
                
        a = 0.85e2; b = 0.85e2;
        r = Beta.beta(a, b);
        assertTrue(Math.abs(r - 0) < eps);
        
        a = 0.5; b = 1;
        r = Beta.beta(a, b);
        assertTrue(Math.abs(r - 2.0) < tol);
        
        a = 1; b = 170;
        r = Beta.beta(a, b);
        assertTrue(Math.abs(r - 2.0) < tol);
    }
}
