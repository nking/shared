package algorithms.sampling;

import algorithms.misc.CDFStandardNormal;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

/**
 *
 * @author nichole
 */
public class Gaussian {
    
    /**
     * generate u = nx1 vector where each element u_j is independently sampled from N (0, 1)
            normal distribution with mean=0 and covariance=1
                           1             ( -(x - mu)^2 )
             f = ------------------ * exp( ----------- )
                 sigma * sqrt(2*pi)      (    2o~^2    )
        
                           1             ( -(x)^2 )
             f = ------------------ * exp( ------ )
                     1 * sqrt(2*pi)      (   2    )
     * @return 
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
          using the CDF of the standard normal.
     * @param rand
     * @param n
     * @return 
     */
    public static double[] randomSampleOfUnitStandard(SecureRandom rand, int n) {
        double norm = Math.pow(2.*Math.PI, -0.5);
        double t;
        
        double[] u = new double[n];
        
        int i;
        // u range is -3.1 to +3.1
        for (i = 0; i < n; ++i) {
            t = rand.nextDouble();
            u[i] = CDFStandardNormal.approxShort(t);
        }
        
        return u;
    }
}
