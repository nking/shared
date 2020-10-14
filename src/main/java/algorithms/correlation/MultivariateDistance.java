package algorithms.correlation;

import algorithms.correlation.UnivariateDistance.DCor;
import algorithms.correlation.UnivariateDistance.DCov;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.Shuffle;
import algorithms.statistics.Gamma;
import algorithms.statistics.MultivariateUniformDistribution;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;

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
     * @param k the number of random projections
     * @return 
     */
    public static double efficientDCov(double[][] x, double[][] y, int k) throws NoSuchAlgorithmException {
        // number of columns in X is p
        // number of columns in Y is q.
        int p = x[0].length;
        int q = y[0].length;
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        double CpCq = _calcCapitalC(p) * _calcCapitalC(q);
        
        double[] u, v, xu, yv;
        DCov dcov;
        double t;
        double meanT = 0;
        for (int i = 0; i < k; ++i) {
            u = MultivariateUniformDistribution.generateOnUnitStandardNSphere(p, rand);
            v = MultivariateUniformDistribution.generateOnUnitStandardNSphere(q, rand);
            
            // x * u^T = nX1
            xu = MatrixUtil.multiply(x, u);
            yv = MatrixUtil.multiply(y, v);
            
            dcov = UnivariateDistance.fastDcov(xu, yv);
            
            t = CpCq * dcov.covsq;
            meanT += t;
        }
        meanT /= (double)k;
        
        return meanT;
    }
    
    /**
     * test for independence of x and y using permutations of y and the efficient
     * dCov as a statistic.
     * 
     * @param x
     * @param y
     * @param k the number of random projections for each test statistic.
     * @param nIterations the number of iterations for statistic calculations
     * (note that each iteration constructs a new permutation of y, so this
     * step has runtime complexity O(y.length * y[0].length)
     * @param alpha significance level for testing null hypothesis
     * @return
     * @throws NoSuchAlgorithmException 
     */
    public static boolean areIndependent1(double[][] x, double[][] y, 
        int k, int nIterations, double alpha) throws NoSuchAlgorithmException {
           
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        int p = x.length;
        int q = y.length;
        
        int i, j;
        
        int[] yIdx = new int[q];
        for (i = 0; i < q; ++i) {
            yIdx[i] = i;
        }
        double[][] y2;
        
        double t = efficientDCov(x, y, k);
        
        double s = 0;
        
        double[] t2 = new double[nIterations];
        
        for (i = 0; i < nIterations; ++i) {
            
            // permute each column of a copy of y
            y2 = MatrixUtil.transpose(y);
            for (j = 0; j < y2.length; ++j) {
                Shuffle.fisherYates(y2[j], rand);
            }
            y2 = MatrixUtil.transpose(y2);
            
            t2[i] = efficientDCov(x, y2, k);
            
            if (t > t2[i]) {
                s++;
            }
        }
        
        s = (1. + s)/(1. + nIterations);
        
        return (s > alpha);
    }
        
    /**
     * test for independence of x and y using threshold of an approximate 
     * asymptotic distribution
     * 
     * @param x
     * @param y
     * @param k the number of random projections for each test statistic.
     * @param alpha significance level for testing null hypothesis
     * @return
     * @throws NoSuchAlgorithmException 
     */
    public static boolean areIndependent2(double[][] x, double[][] y, 
        int k, double alpha) throws NoSuchAlgorithmException {
        
        if (x.length != y.length) {
            throw new IllegalArgumentException("x.lenght must equal y.length");
        }
           
        // number of columns in X is p
        // number of columns in Y is q.
        int p = x[0].length;
        int q = y[0].length;
        
        int n = x.length;
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        double Cp = _calcCapitalC(p);
        double Cq = _calcCapitalC(q);   
        
        double[] u, v, xu, yv, uPrime, vPrime, xuPrime, yvPrime;
        DCor dcor;
        DCov dcov;
        
        int i, j;
                        
        // txy = Ω(k) = C_p * C_q * Ω(u_k*X, v_k*Y)
        // sxxyy = S_{n,1}(k) = (C_p * C_q)^2 * Ω(u_k*X, v_k*X) * * Ω(u_k*Y, v_k*Y)
        // s2 = S_{n,2}(k) = C_p * (aDotDot_{u_k))/(n*(n-1))
        // s3 = S_{n,3}(k) = C_q * (bDotDot_{v_k))/(n*(n-1))
        // txx = ==> Ω_X(k) = C_p^2 * Ω(u_k*X, u_k_prime*X) <====
        // tyy = ==> Ω_Y(k) = C_q^2 * Ω(v_k*Y, v_k_prime*Y) <====
        
        double txy = 0;
        double sxxyy = 0;
        double s2 = 0;
        double s3 = 0;
        double txx = 0;
        double tyy = 0;
        
        for (i = 0; i < k; ++i) {
            u = MultivariateUniformDistribution.generateOnUnitStandardNSphere(p, rand);
            v = MultivariateUniformDistribution.generateOnUnitStandardNSphere(q, rand);
            
            // x * u^T = nX1
            xu = MatrixUtil.multiply(x, u);
            yv = MatrixUtil.multiply(y, v);
            
            dcor = UnivariateDistance.fastDcor(xu, yv);
            
            uPrime = MultivariateUniformDistribution.generateOnUnitStandardNSphere(p, rand);
            vPrime = MultivariateUniformDistribution.generateOnUnitStandardNSphere(q, rand);
            xuPrime = MatrixUtil.multiply(x, uPrime);
            yvPrime = MatrixUtil.multiply(y, vPrime);
            
            // txy = Ω(k) = C_p * C_q * Ω(u_k*X, v_k*Y)
            // sxxyy = S_{n,1}(k) = (C_p * C_q)^2 * Ω(u_k*X, u_k*X) * * Ω(v_k*Y, v_k*Y)
            // s2 = S_{n,2}(k) = C_p * (aDotDot_{u_k))/(n*(n-1))
            // s3 = S_{n,3}(k) = C_q * (bDotDot_{v_k))/(n*(n-1))
            // txx = ==> Ω_X(k) = C_p^2 * Ω(u_k*X, u_k_prime*X) <====
            // tyy = ==> Ω_Y(k) = C_q^2 * Ω(v_k*Y, v_k_prime*Y) <====  
            
            txy += (dcor.covXYSq.covsq);
            sxxyy += (dcor.covXXSq.covsq * dcor.covYYSq.covsq);
            s2 += dcor.covXYSq.aDotDot;
            s3 += dcor.covXYSq.bDotDot;
            
            dcov = UnivariateDistance.fastDcov(xu, xuPrime);
            txx += dcov.covsq;
            dcov = UnivariateDistance.fastDcov(yv, yvPrime);
            tyy += dcov.covsq;
        }
        
        double invK = 1./k;
        
        txy *= Cp*Cq;
        sxxyy *= Math.pow(Cp*Cq, 2.);
        s2 *= (Cp/((double)(n*(n-1.))));
        s3 *= (Cq/((double)(n*(n-1.))));
        txx *= Cp*Cp;
        tyy *= Cq*Cq;
        
        txy *= invK;
        sxxyy *= invK;
        s2 *= invK;
        s3 *= invK;
        txx *= invK;
        tyy *= invK;
        
        double numer = s2 * s3;
        double denom = ((2.*(k-1.)/(double)k) * txx * tyy) + ((sxxyy)/(double)k);
        
        double betaT = numer / denom;        
        double alphaT = numer * betaT;
        
        double t = efficientDCov(x, y, k);
        
        //Reject H0 if n*t + s2*s3 > Gamma(alphaT, betaT; 1 - alpha);
        
        throw new UnsupportedOperationException("not yet implemented");
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
