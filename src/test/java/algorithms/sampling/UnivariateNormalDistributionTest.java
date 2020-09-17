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
        int n = 1000;
        double h = 0.5;
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        // histogram for testing that null hypothesis is true or can't
        //    be proven...
        
        // K-S test
        
        double[] u = UnivariateNormalDistribution.randomSampleOf(mean, sigma,
            rand, n);
        
        double[] avgAndStDev = MiscMath0.getAvgAndStDev(u);
        
        // 68% in [-1, 1], 95% in [-2, 2], 99.7% within [-3, 3]
        // [-2, 2] ==> [123 - 2*23.45, 123 + 2*23.45] = 
        //    
        //    with bin size = h
        HistogramHolder fh = Histogram.createSimpleHistogram(u, h, 76, 170);
        
        // create CDF from it
        int[] fhSum = MiscMath0.cumulativeSum(fh.getYHist());
        double[] cdf = new double[fhSum.length];
        for (int i = 0; i < fhSum.length; ++i) {
            cdf[i] = (double)fhSum[i]/(double)fhSum[fhSum.length - 1];
        }
        
        double[] expectedCDF = generateGaussianCDF(fh.getXHist(), mean, sigma);
        
        // largest difference w/ expected gaussian(mean, sigma)
        double ksStat = Double.NEGATIVE_INFINITY;
        double diff, y;
        for (int i = 0; i < fhSum.length; ++i) {
            y = expectedCDF[i];
            diff = Math.abs(cdf[i] - y);
            if (diff > ksStat) {
                ksStat = diff;
            }
        }
        
        //H0:  the data are normally distributed
        //Ha:  the data are not normally distributed
        //    reject H0 if ksstat > crit.
        
        // for N>35, can use  Smirnov (1948)
        //     which is summarized in https://blogs.sas.com/content/iml/2019/05/20/critical-values-kolmogorov-test.html
        // 
        // For 95% level, ks statistic critical value = 1.36/sqrt(n)
        //     99% level, ks statistic critical value = 1.63/sqrt(n)
        //     90% level, ks statistic critical value = 1.22/sqrt(n)
        //     85% level, ks statistic critical value = 1.14/sqrt(n)
        //     80% level, ks statistic critical value = 1.07/sqrt(n)
        
        // if ksStat < crit, do not reject null hypothesis (which is that
        //    the generated samples come from a gaussian distribution of mean, sigma.
        double crit = 1.36/Math.sqrt(n);
        
        assert(ksStat < crit);   
    }

    private double[] generateGaussianCDF(float[] xHist, double mean, double sigma) {

        double inv2 = -(1./(sigma*sigma*2.));
        double inv1 = (1./(sigma*Math.sqrt(2.*Math.PI)));
        
        float x;
        double[] y = new double[xHist.length];
        for (int i = 0; i < y.length; ++i) {
            x = xHist[i];
            y[i] = inv1*(Math.exp(inv2*(x - mean)*(x - mean)));
        }
        
        double[] ySum = MiscMath0.cumulativeSum(y);
        double[] cdf = new double[ySum.length];
        for (int i = 0; i < ySum.length; ++i) {
            cdf[i] = ySum[i]/ySum[ySum.length - 1];
        }
        return cdf;
    }
    
}
