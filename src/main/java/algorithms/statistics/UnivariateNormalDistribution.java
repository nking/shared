package algorithms.statistics;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Random;

/**
 *
 * @author nichole
 */
public class UnivariateNormalDistribution {
    
    /**
     * generate u = nx1 vector where each element u_j is independently sampled from N (0, 1)
            normal distribution with mean=0 and variance=1
                           1             ( -(x - mu)^2 )
             f = ------------------ * exp( ----------- )
                 sigma * sqrt(2*pi)      (    2o~^2    )
        
                           1             ( -(x)^2 )
             f = ------------------ * exp( ------ )
                     1 * sqrt(2*pi)      (   2    )
         using the inverse CDF of the standard normal.
     @param n
     @return 
     * @throws java.security.NoSuchAlgorithmException 
     */
    public static double[] randomSampleOfUnitStandard(int n) throws NoSuchAlgorithmException {
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        return randomSampleOfUnitStandard(rand, n);
    }
    
    /**
     * generate u = nx1 vector where each element u_j is independently sampled from N (0, 1)
            normal distribution with mean=0 and variance=1
                           1             ( -(x - mu)^2 )
             f = ------------------ * exp( ----------- )
                 sigma * sqrt(2*pi)      (    2o~^2    )
        
                           1             ( -(x)^2 )
             f = ------------------ * exp( ------ )
                     1 * sqrt(2*pi)      (   2    )
          using the inverse CDF of the standard normal.
     @param rand
     @param n
     @return n randomly sampled x's
     */
    public static double[] randomSampleOfUnitStandard(Random rand, int n) {
        double t;
        
        double[] u = new double[n];
        
        int i;
        // u range is approximately -3.1 to +3.1
        for (i = 0; i < n; ++i) {
            t = rand.nextDouble();
            u[i] = CDFStandardNormal.approxInverseShort(t);
        }
        
        return u;
    }
    
    /**
     * return a random sample of size n of a gaussian distribution that has 
     * the given mean and sigma.
     @param mean the location parameter of the gaussian.  it's the mean.
     @param sigma the shape parameter of the gaussian.  it's the standard deviation.
     @param rand
     @param n
     @return 
     */
    public static double[] randomSampleOf(double mean, double sigma,
                                          Random rand, int n) {
                
        double[] u = randomSampleOfUnitStandard(rand, n);
               
        int i;
        for (i = 0; i < n; ++i) {
            u[i] *= sigma;
            u[i] += mean;
        }
        
        return u;
    }
}
