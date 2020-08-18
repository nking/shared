package algorithms.pca;

import algorithms.correlation.BruteForce;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import java.util.Arrays;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 * find the principal components, that is, the minimal residual variance bases
 * of vectors of samples following G. Strang's SVD in machine learning.
 * 
 * @author nichole
 */
public class PrincipalComponents {
    
    /**
     * 
     * @param x is a 2-dimensional array of k vectors of length n in format
     *    double[n][k]
     * @param nComponents the number of principal components to return.
     * @return a few statistics of the SVD of the covariance of A, up to the
     * nComponents dimension.  Note that if the rank of the SVD(cov(A)) is
     * less than nComponents, then only that amount is returned;
     */
    public static PCAStats calcPrincipalComponents(double[][] x, int nComponents) throws NotConvergedException {
        
        int n = x.length;
        int nDimensions = x[0].length;
        
        if (nComponents > nDimensions) {
            throw new IllegalArgumentException("x has " + nDimensions + " but " 
                + " nComponents principal components were requested.");
        }
        
        /*
        minimal residual variance basis:
           basis direction b.
           samples are vectors x_j where j = 0 to k
           sample mean m_s = (1/k) summation_over_j(x_j)

           looking for a p-dimensional basis that minimizes:
              SSD_p = min_{B}( summation_over_j( min_{a_j}(||x_j - (m_s + B*a_j)||^2) ) )
                 where B = (b1, . . . , bp) is the n × p matrix formed from the selected basis.

              choose the coefficients ⃗aj which minimize the least squares error
                 E_j^2 = ||x_j - (m_s + B*a_j)||^2
              then choose B to minimize the SSD = summation_over_j(E_j^2)

           SSD_p is called the minimum residual variance for any basis of dimension p.

           [U, S, V] = SVD( Cov )
           the 1st principal direction is the 1st column of U.
        */
    
        //TODO: replace this with the faster methods when those are ready for use
        //[x[0].length][x[0].length]
        double[][] cov = BruteForce.covariance(x);
        assert(nDimensions == cov.length);
        assert(nDimensions == cov[0].length);
        
        SVD svd = SVD.factorize(new DenseMatrix(cov));
        
        double[] s = svd.getS();
        
        System.out.println("eigenvalues of cov = " + Arrays.toString(s));
        System.out.flush();
        
        int rank = 0;
        for (int i = 0; i < s.length; ++i) {
            if (Math.abs(s[i]) > 1e-17) {
                rank++;
            }
        }
        
        if (nComponents > rank) {
            nComponents = rank;
        }
        
        // the first nComponents COLUMNS of u are the principal axes
        DenseMatrix u = svd.getU();
        assert(nDimensions == u.numRows());
        assert(nDimensions == u.numColumns());
        
        System.out.println("U of SVD(cov) = " + u.toString());
        System.out.flush();
                        
        double[][] pc = new double[nDimensions][nComponents];
        for (int row = 0; row < u.numRows(); ++row) {
            for (int col = 0; col < nComponents; ++col) {
                pc[row][col] = u.get(col, row);
            }
        }
        
        PCAStats stats = new PCAStats();
        stats.nComponents = nComponents;
        stats.principalDirections = pc;
        
        stats.s = new double[nComponents];
        for (int i = 0; i < nComponents; ++i) {
            stats.s[i] = s[i];
        }
                
        double total = 0;
        for (int j = 0; j < nComponents; ++j) {
            total += stats.s[j];
        }
        double sum = 0;
        int p;
        for (int j = (nComponents+1); j <= nDimensions; ++j) {
            p = j - 1;
            sum += s[p];
        }
        stats.ssdP = sum;
        stats.fractionVariance = (total - sum)/total;
        
        return stats;
    }
    
    /**
     * reconstruct an image from x, given principal directions.
     * <pre>
     *     ⃗r(⃗a_0) = m⃗_s + U_p * ⃗a 
     *         where U_p is the principalDirections matrix,
     *         where m_s is the sample x mean,
     *         where a_0 is
     *            a_0 = arg min_{a} ||⃗x − (m⃗_s + U _p * a)||^2
     * </pre>
     * @param x a 2-dimensional array of k vectors of length n in format
     *    double[n][k]
     * @param principalDirections the principal components derived from
     * SVD of the covariance of the training data.
     * format is [nDimensions][nComponents]
     * @return 
     */
    public static double[][] reconstruct(double[][] x, double[][] principalDirections) {
        
        int n = x.length;
        int nDimensions = x[0].length;
        
        int i, j;
        
        double[] mean = new double[nDimensions];
        double sum;
        for (i = 0; i < nDimensions; ++i) {
            sum = 0;
            for (j = 0; j < n; ++j) {
                sum += (x[j][i]);
            }
            mean[i] = sum/(double)n;
        }
        
        // find the 'a' which minimizes ||⃗x − (m⃗_s + U _p * a)||^2
        
        /*
        ||⃗x − (m⃗_s + U _p * a)||^2 
            where m⃗_s, U _p, and 'a' are k dimensional vectors
            and x is a matrix of size [n][dimensions possibly larger than k]

        f = (x_{i:n, 0} − (m⃗_0 + U_p_0 * a_0))^2
             + (x_{i:n, 1} − (m⃗_1 + U_p_1 * a_1))^2
             + ... up to dimension of principal components

        df/da_0 = 2 * U_p_0 * (x_{i:n, 0} − (m⃗_0 + U_p_0 * a_0))
        df/da_1 = 2 * U_p_1 * (x_{i:n, 1} − (m⃗_1 + U_p_1 * a_1))

        set df/da_0 = 0 ==> x_{i:n, 0} = (m⃗_0 + U_p_0 * a_0)
                            a_0 = (x_{i:n, 0} - m⃗_0) / U_p_0
        */
            
        throw new UnsupportedOperationException("not yet implemented");
        
    }
    
    public static class PCAStats {
        /** the number of components requested from the calculation, that is,
         the first p principal components from the U matrix of the SVD 
         decomposition of the covariance of x.
         * NOTE that x given to the algorithm, is a sample of vectors 
         *    {⃗x_j}_{j=1 to k}, with each ⃗x_j ∈ R^n (each x_j is a vector of 
         *    n real numbers) that form the matrix X ∈ R^(nxk).
. 
         */
        int nComponents;
        /**
         * the first nComponents principal directions of x (the bases that
         * minimize the variance).
         */
        double[][] principalDirections;
        
        /**
         * the fraction of the total variance Q_p (where p is the dimension
         * number, that is, nComponents):
         * <pre>
         *    Q_p = (SSD_0 − SSD_p)/SSD_0
         * 
         *    where SSD_p = summation_{j=p+1 to n}(s_j)
         *       where s_j is from the diagonal matrix S of the SVD
         *          of the covariance of x.
         * </pre>
         */
        double fractionVariance;
        
        /**
         * Minimum residual variance
         */
        double ssdP;
        
        /**
         * the first nComponents singular values from the diagonal matrix of
         * the SVD of covariance of x;
         */
        double[] s;
    }
}
