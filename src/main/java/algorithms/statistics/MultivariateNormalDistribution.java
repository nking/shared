package algorithms.statistics;

import algorithms.matrix.MatrixUtil;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import no.uib.cipr.matrix.DenseCholesky;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.LowerSPDDenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.NotConvergedException;

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
     * @return a fair sampling from a normal distribution N(M, K).
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static double[] sampleFrom0(double[] m, double[][] k) throws NotConvergedException, NoSuchAlgorithmException {
        
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
        double[] u = UnivariateNormalDistribution.randomSampleOfUnitStandard(rand, n);
        //System.out.println("u="+Arrays.toString(u));
        
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
     * caveat is that the Cholesky decomposition may fail for ill-conditioned
     * matrices k.
     * @param m vector of means for multivariate distribution.
     * @param k double array of covariance matrix for multivariate distribution.
     * must be a symmetric positive definite matrix.
     * @return a fair sampling from a normal distribution N(M, K).
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static double[] sampleFrom1(double[] m, double[][] k) throws NotConvergedException, NoSuchAlgorithmException {
        
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
        multivariate PDF:
           ((2*pi)^(-k/2)) * (det(k)^(-1/2)) * exp(-0.5*( (x-m)^T * (k^-1) * (x-m) ))
        */
        
        // 1D: x = m + stDev * N(0, 1)
        
        /*
        To sample from the Normal distribution N(m⃗, K) we do the following:
        (K is square, symmetric, and positive definite)

        1)  generate u = nx1 vector where each element u_j is independently sampled from N (0, 1)
            normal distribution with mean=0 and covariance=1
                           1             ( -(x - mu)^2 )
             f = ------------------ * exp( ----------- )
                 sigma * sqrt(2*pi)      (    2o~^2    )
        
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
        */        
        double[] u = UnivariateNormalDistribution.randomSampleOfUnitStandard(rand, n);
        
        //System.out.println("u="+Arrays.toString(u));
        
        for (i = 0; i < n; ++i) {
            k[i][i] += 0.0001;
        }
        
        DenseCholesky c = new DenseCholesky(n, false);
        LowerSPDDenseMatrix lt = new LowerSPDDenseMatrix(new DenseMatrix(k));
        c.factor(lt);
        double[] d = MatrixUtil.multiplyMatrixByColumnVector(lt, u);
        
        assert(d.length == m.length);
        
        {
            // check
            double[][] lu = Matrices.getArray(lt);
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
        }
        
        double[] x = new double[n];
        for (i = 0; i < n; ++i) {
            x[i] = m[i] + d[i];
        }
        
        return x;
    }
    
}
