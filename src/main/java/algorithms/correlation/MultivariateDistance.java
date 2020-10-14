package algorithms.correlation;

import algorithms.matrix.MatrixUtil;
import algorithms.statistics.Gamma;
import algorithms.statistics.MultivariateUniformDistribution;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

/**
 *
 * @author nichole
 */
public class MultivariateDistance {
    
    /**
     * calculate the distance covariance using the average of 
     * fast univariate distance covariances of random projections.
     * 
     <pre>
      following the algorithm
      “A Statistically And Numerically Efficient Independence Test Based On 
      Random Projections And Distance Covariance”, 
      2017, Cheng Huang, And Xiaoming Huo, Annals of Statistics
     </pre>
     <pre>
     runtime complexity is O(k*n*log_2(n))
     where k is the number of random projections and n is the sample size.
     memory requirement is O(max{n, k}).
     </pre>
     * @param x multivariate variable where the columns are the variates and 
     * rows are the samples.
     * @param y multivariate variable where columns are the variates and 
     * rows are the samples.
     * @return 
     */
    public static double efficientDCov(double[][] x, double[][] y, double k) throws NoSuchAlgorithmException {
        // number of columns in X is p
        // number of columns in Y is q.
        int p = x[0].length;
        int q = y[0].length;
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        double[] u, v, xu, yv;
        for (int i = 0; i < k; ++i) {
            u = MultivariateUniformDistribution.generateOnUnitStandardNSphere(p, rand);
            v = MultivariateUniformDistribution.generateOnUnitStandardNSphere(q, rand);
            
            // x * u^T = nX1
            xu = MatrixUtil.multiply(x, u);
            yv = MatrixUtil.multiply(y, v);
        }
        
        throw new UnsupportedOperationException("not implemented yet");
    }
    
    public static double _calcC(double a) {
        double b = (a + 1.)/2.;
        double numer = Math.pow(Math.PI, b);
        double denom = Gamma.lanczosGamma9(b);
        return numer/denom;
    }
    
    public static double _calcCapitalC(double a) {
        double b1 = (a + 1.)/2.;
        double b2 = a/2.;
        double numer = Math.sqrt(Math.PI) * Gamma.lanczosGamma9(b1);
        double denom = Gamma.lanczosGamma9(b2);
        return numer/denom;
    }
}
