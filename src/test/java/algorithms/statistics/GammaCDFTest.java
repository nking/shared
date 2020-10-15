package algorithms.statistics;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GammaCDFTest extends TestCase {
    
    public GammaCDFTest(String testName) {
        super(testName);
    }
    
    public void testCDF() {
        
        double tol=1.e-4;
        
        // example values from https://www.mathworks.com/help/stats/gamcdf.html
        //   which agree with scipy.stats.gamma.cdf
        double[] a = new double[]{1, 2, 3, 4, 5, 6};
        double[] b = new double[]{5, 6, 7, 8, 9, 10};

        double[] expected = new double[]{0.6321, 0.5940, 0.5768, 0.5665, 0.5595, 0.5543};
        
        double p, diff;
        for (int i = 0; i < a.length; ++i) {
            p = GammaCDF.cdf(a[i]*b[i], a[i], b[i]);
            diff = Math.abs(p - expected[i]);
            //System.out.printf("diff=%.7e  x=%.1e; a=%.1e;  b=%.1e;  p=%.4e\n", diff, 
            //    a[i]*b[i], a[i], b[i], p);
            assertTrue(diff < tol);
        }
    }
    
    public void testInverseCDF() {
        
        double shape = 2;
        double scale = 2;
        
        double alpha = 0.999;
        double expectedX = 18.55;
        
        double x = GammaCDF.inverseCdf(shape, scale, alpha);
        
        double tol = 1.e-1;
        assertTrue(Math.abs(x - expectedX) < tol);
        
    }
}
