package algorithms.misc;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ComplexTest extends TestCase {
    
    public ComplexTest(String testName) {
        super(testName);
    }

    public void test0() {
        
        Complex a = new Complex(5, 6);
        Complex b = new Complex(-3, 4);
        
        double eps = 1e-5;
        
        assertTrue(Math.abs(5 - a.re()) < eps);
        assertTrue(Math.abs(6 - a.im()) < eps);
        
        Complex t = a.plus(b);
        assertTrue(Math.abs(2 - t.re()) < eps);
        assertTrue(Math.abs(10 - t.im()) < eps);
        
        t = a.minus(b);
        assertTrue(Math.abs(8 - t.re()) < eps);
        assertTrue(Math.abs(2 - t.im()) < eps);
        
        t = a.times(b);
        assertTrue(Math.abs(-39 - t.re()) < eps);
        assertTrue(Math.abs(2 - t.im()) < eps);
        
        t = b.times(a);
        assertTrue(Math.abs(-39 - t.re()) < eps);
        assertTrue(Math.abs(2 - t.im()) < eps);
        
        t = a.divided(b);
        assertTrue(Math.abs(0.36 - t.re()) < eps);
        assertTrue(Math.abs(-1.52 - t.im()) < eps);
        
        t = t.times(b);
        assertTrue(Math.abs(5 - t.re()) < eps);
        assertTrue(Math.abs(6 - t.im()) < eps);
        
        t = a.conjugate();
        assertTrue(Math.abs(5 - t.re()) < eps);
        assertTrue(Math.abs(-6 - t.im()) < eps);
        
        double eps2 = 1e-7;
        
        double r = a.abs();
        assertTrue(Math.abs(7.810249675906654 - r) < eps2);
        
        t = a.tan();
        assertTrue(Math.abs(-6.685231390246571E-6 - t.re()) < eps2);
        assertTrue(Math.abs(1.0000103108981198 - t.im()) < eps2);
        
        a = new Complex(8, 0);
        t = a.nthRoot(3);
        assertTrue(Math.abs(2 - t.abs()) < eps2);
        // phases 0, 360/3., 720./3, 1080/3. ...
        assertTrue(Math.abs(0 - t.phase()) < eps2);
        
        a = new Complex(1, 1);
        t = a.nthRoot(1./8.);
        assertTrue(Math.abs(16 - t.re()) < eps2);
        assertTrue(Math.abs(0 - t.im()) < eps2);
        
        a = new Complex(-1, 0);
        t = a.naturalLog();
        assertTrue(Math.abs(0 - t.re()) < eps2);
        assertTrue(Math.abs(Math.PI - t.im()) < eps2);
        double magn = t.abs();
        double phase = t.phase();
        
        a = new Complex(0, 1);
        b = new Complex(0, -2);
        t = a.power(b);
        assertTrue(Math.abs(23.14 - t.abs()) < 0.01);
        
        a = new Complex(0, 1);
        b = new Complex(0.5, 0);
        t = a.power(b);
        double s = 1./Math.sqrt(2);
        // (1+i)/sqrt(2)
        assertTrue(Math.abs(s - t.re()) < eps2);
        assertTrue(Math.abs(s - t.im()) < eps2);
    }
    
}
