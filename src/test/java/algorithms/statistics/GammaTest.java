package algorithms.statistics;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GammaTest extends TestCase {
    
    public GammaTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        // compare to python's math.gamma
        //Math.exp(logGamma(x)
        double[] expectedLogGamma = new double[] {
            0, 0, 0, 
            0.693147, 1.791759469, 3.1780538, 4.78749174,
            6.579251212, 8.52516136, 10.6046029, 12.801827
        };
        double[] expectedGamma = new double[] {
            0, 1, 1,
            2, 6, 24, 120, 720.0, 5040.0, 40320.0, 362880.0
        };
        double e, el, g0, gl0, g1, gl1; 
        double diff;
        double tol = 1e-3;
        for (int i = 1; i < 11; ++i) {
            e = expectedGamma[i];
            el = expectedLogGamma[i];
            
            g0 = Gamma.lanczosGamma9(i);
            gl0 = Gamma.lanczosLGamma9(i);
            
            g1 = Gamma.stirlingGamma(i);
            gl1 = Gamma.stirlingLGamma(i);
            
            //System.out.printf("%d | %11.4f: %11.4f  %11.4f | %11.4f: %11.4f  %11.4f%n",
            //    i, e, g0, g1, el, gl0, gl1);
            //System.out.flush();
            
            diff = Math.abs(e - g0);
            assertTrue(diff < tol);            
            diff = Math.abs(el - gl0);
            assertTrue(diff < tol);
            
            diff = Math.abs(e - g1);
            assertTrue(diff < 1.e-1);            
            diff = Math.abs(el - gl1);
            assertTrue(diff < tol);
            
        }
    }
    
}
