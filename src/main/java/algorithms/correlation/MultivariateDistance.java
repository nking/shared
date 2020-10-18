package algorithms.correlation;

import algorithms.correlation.UnivariateDistance.DCor;
import algorithms.correlation.UnivariateDistance.DCov;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import algorithms.misc.Shuffle;
import algorithms.statistics.Gamma;
import algorithms.statistics.GammaCDF;
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
      https://arxiv.org/pdf/1701.06054.pdf
     </pre>
     <pre>
     runtime complexity is O(n * K * log_2(n))
     (more specifically: O(n * K * (log_2(n) + p + q)))
     where k is the number of random projections and n is the sample size.
     memory requirement is O(max{n, K}).
     </pre>
     NOTE: there is material from the paper by Huang and Huo, and more that is 
     related to it that is in the thesis of Huang:
     "Some computationally efficient methods in statistics and their 
     applications in parameter estimation and hypotheses testing"
     https://smartech.gatech.edu/bitstream/handle/1853/60120/HUANG-DISSERTATION-2017.pdf
     * @param x multivariate variable where the columns are the variates and 
     * rows are the samples.
     * @param y multivariate variable where columns are the variates and 
     * rows are the samples.
     * @param k the number of random projections
     * @return 
     */
    public static double efficientDCov(double[][] x, double[][] y, int k) throws NoSuchAlgorithmException {
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        return efficientDCov(x, y, k, rand);
    }
    
    /**
     * calculate the distance covariance using the average of 
     * fast univariate distance covariances of random projections.
     * 
     <pre>
      following the algorithm
      “A Statistically And Numerically Efficient Independence Test Based On 
      Random Projections And Distance Covariance”, 
      2017, Cheng Huang, And Xiaoming Huo, Annals of Statistics
      https://arxiv.org/pdf/1701.06054.pdf
     </pre>
     <pre>
     runtime complexity is O(n * K * log_2(n))
     (more specifically: O(n * K * (log_2(n) + p + q)))
     where k is the number of random projections and n is the sample size.
     memory requirement is O(max{n, K}).
     </pre>
     NOTE: there is material from the paper by Huang and Huo, and more that is 
     related to it that is in the thesis of Huang:
     "Some computationally efficient methods in statistics and their 
     applications in parameter estimation and hypotheses testing"
     https://smartech.gatech.edu/bitstream/handle/1853/60120/HUANG-DISSERTATION-2017.pdf
     * @param x multivariate variable where the columns are the variates and 
     * rows are the samples.
     * @param y multivariate variable where columns are the variates and 
     * rows are the samples.
     * @param k the number of random projections
     * @param rand instance of secure random number generator
     * @return 
     */
    public static double efficientDCov(double[][] x, double[][] y, int k,
        SecureRandom rand) {
        
        if (k < 1) {
            throw new IllegalArgumentException("k must be a positive number grater than 0");
        }
        // number of columns in X is p
        // number of columns in Y is q.
        int p = x[0].length;
        int q = y[0].length;
        
        double CpCq = _calcCapitalC(p) * _calcCapitalC(q);
        
        double[] u, v, xu, yv;
        DCov dcov;
        double meanT = 0;
        for (int i = 0; i < k; ++i) {
            u = MultivariateUniformDistribution.generateOnUnitStandardNSphere(p, rand);
            v = MultivariateUniformDistribution.generateOnUnitStandardNSphere(q, rand);
            
            // x * u^T = nX1
            xu = MatrixUtil.multiply(x, u);
            yv = MatrixUtil.multiply(y, v);
            
            dcov = UnivariateDistance.fastDcov(xu, yv);
            
            meanT += dcov.covsq;
        }
        meanT *= CpCq/(double)k;
        
        return meanT;
    }
    
    /**
     * test for independence of x and y using permutations of y (approximating the null distribution) 
     and the efficient  dCov as a statistic.
     <pre>
      following the algorithm
      “A Statistically And Numerically Efficient Independence Test Based On 
      Random Projections And Distance Covariance”, 
      2017, Cheng Huang, And Xiaoming Huo, Annals of Statistics
     </pre>
     <pre>
     runtime complexity is 
     </pre>
     * @param x.  x.length must be >= 20
     * @param y.  x.length must be >= 20
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
        
        return areIndependent1(x, y, k, nIterations, alpha, rand);
    }
    
    /**
     * test for independence of x and y using permutations of y (approximating the null distribution) 
     and the efficient  dCov as a statistic.
     <pre>
      following the algorithm
      “A Statistically And Numerically Efficient Independence Test Based On 
      Random Projections And Distance Covariance”, 
      2017, Cheng Huang, And Xiaoming Huo, Annals of Statistics
     </pre>
     <pre>
     runtime complexity is 
     </pre>
     * @param x.  x.length must be >= 20
     * @param y.  x.length must be >= 20
     * @param k the number of random projections for each test statistic.
     * @param nIterations the number of iterations for statistic calculations
     * (note that each iteration constructs a new permutation of y, so this
     * step has runtime complexity O(y.length * y[0].length)
     * @param alpha significance level for testing null hypothesis
     * @param rand
     * @return if true, x and y are consistent with independent, else if false
     * x and y are not consistent with independent.
     */
    public static boolean areIndependent1(double[][] x, double[][] y, 
        int k, int nIterations, double alpha, SecureRandom rand) {
        
        if (nIterations < 2) {
            throw new IllegalArgumentException("nIterations should be > 1");
        }
        
        int p = x.length;
        int q = y.length;
        
        int i, j;
        
        double[][] y2;

        double t = efficientDCov(x, y, k, rand);
        double t2;
        double s = 0;
        
        for (i = 0; i < nIterations; ++i) {
            
            // permute each column of a copy of y
            y2 = MatrixUtil.transpose(y);
            for (j = 0; j < y2.length; ++j) {
                Shuffle.fisherYates(y2[j], rand);
            }
            y2 = MatrixUtil.transpose(y2);
            
            t2 = efficientDCov(x, y2, k, rand);
            
            // t2 should have cov ~ 0
            // if dependent, t > 0
            System.out.printf("   t=%.4e, t2=%.4e\n", t, t2);
            if (t > t2) {
                s++;
            }
        }
        
        s = (1. + s)/(1. + nIterations);
        
        System.out.printf("t=%.4e, s=%.4e,  1.-alpha=%.4e\n", t, s, 1.-alpha);
        System.out.flush();
        //reject the independence hypothesis (H0) when s is smaller than critical level α (which is 1-α in this case).
        return (s < (1.-alpha));
    }
        
    /**
     * NOT READY FOR USE.
     * test for independence of x and y using threshold of an approximate 
     * asymptotic distribution
     <pre>
      following the algorithm
      “A Statistically And Numerically Efficient Independence Test Based On 
      Random Projections And Distance Covariance”, 
      2017, Cheng Huang, And Xiaoming Huo, Annals of Statistics
     </pre>
     <pre>
      The authors note that this assymptotic dependence test has less power for 
      low dimensional dependency in high dimensional data.
     </pre>
     * @param x.  x.length must be >= 20
     * @param y.  x.length must be >= 20
     * @param k the number of random projections for each test statistic.
     * @param alpha significance level for testing null hypothesis
     * @return
     * @throws NoSuchAlgorithmException 
     */
    public static boolean areIndependent2(double[][] x, double[][] y, 
        int k, double alpha) throws NoSuchAlgorithmException {
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        return areIndependent2(x, y, k, alpha, rand);
    }
    
    /**
     * NOTE READY FOR USE.
     * 
     * test for independence of x and y using threshold of an approximate 
     * asymptotic distribution
     <pre>
      following the algorithm
      “A Statistically And Numerically Efficient Independence Test Based On 
      Random Projections And Distance Covariance”, 
      2017, Cheng Huang, And Xiaoming Huo, Annals of Statistics
     </pre>
     <pre>
      The authors note that this assymptotic dependence test has less power for 
      low dimensional dependency in high dimensional data.
     </pre>
     * @param x.  x.length must be >= 20
     * @param y.  x.length must be >= 20
     * @param k the number of random projections for each test statistic.
     * @param alpha significance level for testing null hypothesis
     * @param rand
     * @return
     */
    public static boolean areIndependent2(double[][] x, double[][] y, 
        int k, double alpha, SecureRandom rand) {
        
        if (x.length != y.length) {
            throw new IllegalArgumentException("x.lenght must equal y.length");
        }
           
        // number of columns in X is p
        // number of columns in Y is q.
        int p = x[0].length;
        int q = y[0].length;
        
        int n = x.length;
        
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
        
        double dcorsq = 0;
        
        for (i = 0; i < k; ++i) {
            u = MultivariateUniformDistribution.generateOnUnitStandardNSphere(p, rand);
            v = MultivariateUniformDistribution.generateOnUnitStandardNSphere(q, rand);
            
            // x * u^T = nX1
            xu = MatrixUtil.multiply(x, u);
            yv = MatrixUtil.multiply(y, v);
            
            dcor = UnivariateDistance.fastDcor(xu, yv);
            dcorsq += dcor.corSq;
            
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
        
        System.out.printf("   k=%d txy=%.4e sxxyy=%.4e s2=%.4e s3=%.4e txx=%.4e tyy=%.4e\n", 
            k, txy, sxxyy, s2, s3, txx, tyy);
        
        double numer = s2 * s3;
        double denom = ((k-1.)/(double)k) * txx * tyy;
        denom += sxxyy/(double)k;
        denom *= 2.;
        
  // temporary fudge that may be introducing a Type III error.  alphaT seems too large, so exploring the normalization of denom first:
  denom *= k;
  
        System.out.printf("   =>numer=%4e : s2=%.4e s3=%.4e\n", numer, s2, s3);
        System.out.printf("   =>denom=%4e : 1st=%.4e 2nd=%.4e\n", denom, 
            ((k-1.)/(double)k) * txx * tyy, sxxyy/(double)k);
        
        //NOTE: alphaT seems too large
        double betaT = numer / denom;        
        double alphaT = numer * betaT;
        
        //double t = efficientDCov(x, y, k, rand);
        
        // see 3.3 and 3.5 in The Distance Correlation Chi-Square Test of Shen and Vogelstein
        //   and note that when corsq = 1, univariate samples are dependent
        
        
        //Reject independence (==H0) if n*t + s2*s3 > Gamma(alphaT, betaT; 1 - alpha);
                
        double g = GammaCDF.inverseCdf(alphaT, betaT, 1. - alpha);
        
        System.out.printf("?? dcor=%.4e  dcorsq=%.4e  dcorsq/k=%.4e  n*dcorsq/k = %.4e,\n   gamma.inverseCDF(%.3e, %.3e, %.3e) = (%.3e)\n",
            Math.sqrt(dcorsq), dcorsq, dcorsq*invK, n*dcorsq*invK, 
            alphaT, betaT, 1.-alpha, g);
        
        double stat = n*(txy + s2*s3);        
        
        System.out.printf("Cp=%.4e Cq=%.4e t=%.4e, n=%d k=%d s2=%.4e s3=%.4e\n   (stat=%.4e)  gamma.inverseCDF(%.3e, %.3e, %.3e) = (%.3e)\n",
            Cp, Cq, txy, n, k, s2, s3, stat, alphaT, betaT, 1.-alpha, g);
        System.out.flush();
        
        //if (stat > g) {
        if (Math.sqrt(dcorsq) > g) {
            return false;
        }
        return true;        
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
