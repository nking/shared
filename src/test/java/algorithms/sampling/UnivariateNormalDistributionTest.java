package algorithms.sampling;

import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath0;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class UnivariateNormalDistributionTest extends TestCase {
    
    public UnivariateNormalDistributionTest(String testName) {
        super(testName);
    }

    public void testRandomSampleOfUnitStandard_int() throws Exception {
    }

    public void testRandomSampleOfUnitStandard_SecureRandom_int() {
    }
    
    public void testRandomSampleOf() throws NoSuchAlgorithmException {
        
        double mean = 123.;
        double sigma = 23.45;
        int n = 100000;
        double h = 0.1;
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        // histogram for testing that null hypothesis is true or can't
        //    be proven...
        
        // K-S test
        
        double[] u = UnivariateNormalDistribution.randomSampleOf(mean, sigma,
            rand, n);
        
        // histogram from -3.1 to +3.1 in X of N(0,1) is 50.3 to 195.7
        //    
        //    with bin size = h
        HistogramHolder fh = Histogram.createSimpleHistogram(u, h, 100., 150);
        
        // create CDF from it
        int[] fhSum = MiscMath0.cumulativeSum(fh.getYHist());
        double[] cdf = new double[fhSum.length];
        for (int i = 0; i < fhSum.length; ++i) {
            cdf[i] = (double)fhSum[i]/fhSum[fhSum.length - 1];
        }
        
        // largest difference w/ expected gaussian(mean, sigma)
        double ksStat = Double.NEGATIVE_INFINITY;
        double diff, y;
        for (int i = 0; i < fhSum.length; ++i) {
            y = gaussianCDF(fh.getXHist()[i], mean, sigma);
            diff = Math.abs(cdf[i] - y);
            if (diff > ksStat) {
                ksStat = diff;
            }
        }
        
        // 95% level, ks statistic critical value = 1.36/sqrt(n)
        // if ksStat < crit, do not reject null hypothesis (which is that
        //    the generated samples come from a gaussian distribution of mean, sigma.
        double crit = 1.36/Math.sqrt(n);
        
        boolean t = ksStat < crit;
    }
    
    private double gaussianCDF(double x, double mean, double sigma) {
        // from wikipedia:
        // for gaussian, CDF for given x, mean, sigma
        //    = (0.5) * (1 + erf( (x-mean)/(sigma*sqrt(2)) ) )
        double a = Math.abs((x - mean)/(sigma*Math.sqrt(2.)));
        
        double cdf = 0.5 * (1. + erf(a));
      
        return cdf;
    }
    
    private double erf(double x) {
        
        // from wikipedia
        //     an Abramowitz and Stegun algorithm:
        double a1 = 0.278393;   double a2 = 0.230389;
        double a3 = 0.000972; double a4 = 0.078108;
        
        double xsq = x * x;
        double xcb = xsq * x;
        double xq = xcb * x;
        
        double denom = 1. + a1*x + a2*xsq + a3*xcb + a4*xq;
                
        return 1./denom;
    }
    
}
