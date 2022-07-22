package algorithms.statistics;

import junit.framework.TestCase;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

/**
 *
 * @author nichole
 */
public class GumbelCDFTest extends TestCase {

    public GumbelCDFTest(String testName) {
        super(testName);
    }
    
    public void testCDF() {
        
        double tol=1.e-3;
        
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
    
    public void testInverseCDF0() {
        
        double shape = 2;
        double scale = 2;
        
        double alpha = 0.999;
        double expectedX = 17.6;
        
        double tol = 1.e-3;
        double x, p;
        
        x = GammaCDF.inverseCdf(shape, scale, alpha);
        
        assertTrue(x >= (expectedX - tol));
        p = GammaCDF.cdf(x, shape, scale);
        assertTrue(Math.abs(alpha - p) < tol);
        
        //------
        alpha = 0.1;
        expectedX = 1.063;
        
        x = GammaCDF.inverseCdf(shape, scale, alpha);
        
        assertTrue(Math.abs(x - expectedX) < 0.1);
        p = GammaCDF.cdf(x, shape, scale);
        assertTrue(Math.abs(alpha - p) < tol);
    }
    
    public void testInverseCDF1() throws NoSuchAlgorithmException {
        
        int nTests = 100;
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //seed = 1318454033349663L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        double shape, scale, alpha, x, p;
        final double tolP = 1;//1.e-1;
        final double tolX = 1.e-3;
        
        double alphaMin = 0.001;
        double alphaRange = 1. - 2.*alphaMin;
        
        for (int i = 0; i < nTests; ++i) {
            
            alpha = alphaMin + alphaRange * rand.nextDouble();            
            shape = tolX + rand.nextInt(1000) * rand.nextDouble();
            scale = tolX + rand.nextInt(1000) * rand.nextDouble();
            
            System.out.printf("alpha=%.4e, scale=%.4e, shape=%.4e\n",
                alpha, scale, shape);
            
            x = GammaCDF.inverseCdf(shape, scale, alpha);
            
            p = GammaCDF.cdf(x, shape, scale);
            
            System.out.printf("   x=%.4e, p=%.4e  |(alpha - p)|=%.4e\n", x, p, alpha - p);
            
            assertTrue(Math.abs(alpha - p) < tolP);
        }
        
    }
}
