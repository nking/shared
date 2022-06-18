package algorithms.statistics;

import algorithms.correlation.BruteForce;
import algorithms.matrix.MatrixUtil;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;

import algorithms.misc.MiscMath0;
import no.uib.cipr.matrix.*;

/**
 * sample from a multivariate normal distribution N(m⃗ ,K) following
 * Strang's SVD in machine learning.
 * 
 * @author nichole
 */
public class MultivariateNormalDistribution {

    /**
     * sample from a multivariate normal distribution N(m⃗, K) following
     * Strang's SVD in machine learning, that is using sqrt(K) for a factor to
     * the generated unit standard normal distribution.
     * @param m vector of means for multivariate distribution.
     * @param k double array of covariance matrix for multivariate distribution.
     * must be a symmetric positive definite matrix.
     * @param nSamples number of samples to generate of length m.length.  e.g. if m is length 2 and
     *                   sampleSize is 3, the method returns 3 rows of length 2.
     * @return a fair sampling from a normal distribution N(M, K) os length sampleSize.
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static double[][] sampleRandomlyFrom0(double[] m, double[][] k, int nSamples) throws NotConvergedException, NoSuchAlgorithmException {

        double[][] out = new double[nSamples][];
        for (int i = 0; i < nSamples; ++i) {
            out[i] = sampleRandomlyFrom0(m, k);
        }
        return out;
    }

    /**
     * sample from a multivariate normal distribution N(m⃗, K) following
     * Strang's SVD in machine learning, that is using sqrt(K) for a factor to
     * the generated unit standard normal distribution.
     * @param m vector of means for multivariate distribution.
     * @param k double array of covariance matrix for multivariate distribution.
     * must be a symmetric positive definite matrix.
     * @return a fair sampling from a normal distribution N(M, K).
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static double[] sampleRandomlyFrom0(double[] m, double[][] k) throws NotConvergedException, NoSuchAlgorithmException {

        int n = m.length;

        if (n != k.length) {
            throw new IllegalArgumentException("length of m must equal length of k");
        }
        if (!MatrixUtil.isPositiveDefinite(k)) {
            throw new IllegalArgumentException("k must be a positive definite matrix");
        }

        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);

        int i, j;
        double r;

        /*
        Wasserman "All of Statistics", eqn 2.10
        multivariate PDF:
           ((2*pi)^(-k/2)) * (det(k)^(-1/2)) * exp(-0.5*( (x-m)^T * (k^-1) * (x-m) ))
        */

        double[] u = UnivariateNormalDistribution.randomSampleOfUnitStandard(rand, n);
        //System.out.println("u="+Arrays.toString(u));

        return _sampleFrom0(u, m, k);
    }

    /**
     * sample from a multivariate normal distribution N(m⃗, K) following
     * Strang's SVD in machine learning, that is using sqrt(K) for a factor to
     * the generated unit standard normal distribution.
     * The implementation of the sqrt is from Strang's SVD algorithm and
     * from Wasserman's "All of Statistics", 2.43 Theorem:
     * <pre>
     *     X = mu + sqrt(sigma)*Z where mu is the mean vector,
     *     sigma is the covariance, and Z is N(O, I) which is the
     *     unit standard normal distribution.
     * </pre>
     * @param m vector of means for multivariate distribution.
     * @param k double array of covariance matrix for multivariate distribution.
     * must be a symmetric positive definite matrix.
     * @return a fair sampling from a normal distribution N(M, K).
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    static double[] _sampleFrom0(double[] u, double[] m, double[][] k) throws NotConvergedException, NoSuchAlgorithmException {
        
        int n = m.length;
        
        if (n != k.length) {
            throw new IllegalArgumentException("length of m must equal length of k");
        }
        if (!MatrixUtil.isPositiveDefinite(k)) {
            throw new IllegalArgumentException("k must be a positive definite matrix");
        }
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        int i, j;
        double r;
        
        /*
        Wasserman "All of Statistics", eqn 2.10
        multivariate PDF:
           ((2*pi)^(-k/2)) * (det(k)^(-1/2)) * exp(-0.5*( (x-m)^T * (k^-1) * (x-m) ))
        */
        
        // 1D: x = m + stDev * N(0, 1)
        
        /*
        To sample from the Normal distribution N(m⃗, K) we do the following:
        (K is square, symmetric, and positive definite)

        1)  generate u = nx1 vector where each element u_j is independently sampled from N (0, 1)
            normal distribution with mean=0 and variance=1

        2)  Compute the matrix square root of K
        3)  Then d = sqrt(K) * u generates a fair sample from N(0,K).
        4)  Then⃗ x= m⃗ + d, is a fair sample from N(m⃗,K)
        
        */
        //System.out.println("u="+Arrays.toString(u));

        //from Wasserman's "All of Statistics", 2.43 Theorem:
        //N(0, I) ~ Σ^(−1/2) * (X−μ)
        //Σ^(1/2) * N(0, I) ~ (X−μ)
        // d ~ X−μ
        double[][] ksq = MatrixUtil.squareRoot(k);
        double[] d = MatrixUtil.multiplyMatrixByColumnVector(ksq, u);
        assert(d.length == m.length);
        
        /*{
            // check
            System.out.printf("check k=\n");
            for ( i = 0; i < k.length; ++i) {
                for ( j = 0; j < k[i].length; ++j) {
                    System.out.printf("%11.3e  ", k[i][j]);
                }
                System.out.printf("\n");
            }
            System.out.printf("check ksq=\n");
            for ( i = 0; i < ksq.length; ++i) {
                for ( j = 0; j < ksq[i].length; ++j) {
                    System.out.printf("%11.3e  ", ksq[i][j]);
                }
                System.out.printf("\n");
            }
            double[][] cov = MatrixUtil.multiply(MatrixUtil.transpose(ksq), 
                ksq);
            System.out.printf("check cov(ksq^T*ksq)=\n");
            for ( i = 0; i < cov.length; ++i) {
                for ( j = 0; j < cov[i].length; ++j) {
                    System.out.printf("%11.3e  ", cov[i][j]);
                }
                System.out.printf("\n");
            }
            
            System.out.printf("u=\n");
            for ( i = 0; i < u.length; ++i) {
                System.out.printf("%11.3e  ", u[i]);
            }
            System.out.printf("\n");
            System.out.printf("d=\n");
            for ( i = 0; i < d.length; ++i) {
                System.out.printf("%11.3e  ", d[i]);
            }
            System.out.printf("\n");
            System.out.flush();            
        }
        */

        //Σ^(1/2) * N(0, I) ~ (X−μ)
        // d ~ X−μ
        // X ~ d + μ
        double[] x = new double[n];
        for (i = 0; i < n; ++i) {
            x[i] = d[i] + m[i];
        }
        
        return x;
    }
    
    /**
     * sample from a multivariate normal distribution N(m⃗, K) using 
     * the Cholesky decomposition for a factor to
     * the generated unit standard normal distribution.
     * @param m vector of means for multivariate distribution.
     * @param k double array of covariance matrix for multivariate distribution.
     * must be a symmetric positive definite matrix.
     * @return a fair sampling from a normal distribution N(M, K).
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static double[] sampleRandomlyFrom1(double[] m, double[][] k) throws NotConvergedException, NoSuchAlgorithmException {
        
        int n = m.length;
        
        if (n != k.length) {
            throw new IllegalArgumentException("length of m must equal length of k");
        }
        if (!MatrixUtil.isPositiveDefinite(k)) {
            throw new IllegalArgumentException("k must be a positive definite matrix");
        }
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        int i, j;
        double r;

        double[] u = UnivariateNormalDistribution.randomSampleOfUnitStandard(rand, n);

        return _sampleFrom1(u, m, k);
    }

    /**
     *
     * @param u random selection of x from randomSampleOfUnitStandard
     * @param m
     * @param k
     * @return
     * @throws NotConvergedException
     * @throws NoSuchAlgorithmException
     */
    static double[] _sampleFrom1(final double[] u, double[] m, double[][] k) throws NotConvergedException, NoSuchAlgorithmException {

        int n = m.length;

        if (n != k.length) {
            throw new IllegalArgumentException("length of m must equal length of k");
        }
        if (!MatrixUtil.isPositiveDefinite(k)) {
            throw new IllegalArgumentException("k must be a positive definite matrix");
        }

        int i, j;
        double r;

        /*
        multivariate PDF:
           ((2*pi)^(-k/2)) * (det(k)^(-1/2)) * exp(-0.5*( (x-m)^T * (k^-1) * (x-m) ))
        */

        // 1D: x = m + stDev * N(0, 1)

        /*<pre>
        To sample from the Normal distribution N(m⃗, K) we do the following:
        (K is square, symmetric, and positive definite)

        1)                 1             ( -(x - mu)^2 )
             f = ------------------ * exp( ----------- )
                 sigma * sqrt(2*pi)      (    2o~^2    )
        generate u = nx1 vector where each element u_j is independently sampled from N (0, 1)
            normal distribution with mean=0 and covariance=1
                           1             ( -(x)^2 )
             f = ------------------ * exp( ------ )
                     1 * sqrt(2*pi)      (   2    )

        2) add a perturbation to K for numerical stability for the Cholesky
           decomposition (next)
             K = K + 0.0001 * I

        3) calculate L from Cholesky decomposition of K
           (can assert that K = L * L^T
        4) Then x = m + L * u

        NOTE: can use in place of L, the inverse of the upper triangular matrix
        from the Cholesky decomposition.
        </pre>
        */

        //System.out.println("u="+Arrays.toString(u));

        double[][] k2 = MatrixUtil.nearestPositiveSemidefiniteToASymmetric(k, 1e-7);
        /*
        for (i = 0; i < n; ++i) {
            k[i][i] += 0.0001;
        }*/

        DenseCholesky c = new DenseCholesky(n, false);
        LowerSPDDenseMatrix lt = new LowerSPDDenseMatrix(new DenseMatrix(k2));
        c.factor(lt);
        double[] d = MatrixUtil.multiplyMatrixByColumnVector(lt, u);

        assert(d.length == m.length);

        /*{
            // check
            double[][] lu = MatrixUtil.convertToRowMajor(lt);
            double[][] cov = MatrixUtil.multiply(MatrixUtil.transpose(lu),
                    lu);
            System.out.printf("cov check of decomp  L^T * L =\n");
            for ( i = 0; i < cov.length; ++i) {
                for ( j = 0; j < cov[i].length; ++j) {
                    System.out.printf("%11.3e  ", cov[i][j]);
                }
                System.out.printf("\n");
            }
            System.out.flush();
        }*/

        double[] x = new double[n];
        for (i = 0; i < n; ++i) {
            x[i] = m[i] + d[i];
        }

        return x;
    }
}
