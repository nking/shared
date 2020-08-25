package algorithms.pca;

import algorithms.correlation.BruteForce;
import algorithms.matrix.MatrixUtil;
import java.util.Arrays;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 * find the principal components, that is, the minimal residual variance bases
 * of vectors of samples following G. Strang's SVD in machine learning in
 * http://www.cs.toronto.edu/~jepson/csc420/
 * 
 * @author nichole
 */
public class PrincipalComponents {
    
    /**
     * calculate the principal components of the unit standardized data x
     * using SVD.
     * NOTE: should standardize the data before using this method,
     * <pre>
     *      double[] mean = new double[data[0].length];
            double[] stDev = new double[data[0].length];
     * e.g. double[][] x = Standardization.standardUnitNormalization(data, mean, stDev);
     * </pre>
     * 
     * from http://www.cs.toronto.edu/~jepson/csc420/
     * combined with the book by Strang "Introduction to Linear Algebra" 
     * and the book by Leskovec, Rajaraman, and Ullman "Mining of Massive Datasets".
     * also useful:
     * https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
     * and https://online.stat.psu.edu/stat505/book/export/html/670
     * 
     * NOTE: SVD has the drawbacks that a singular vector specifies a linear
       combination of all input columns or rows, and there is a
       a Lack of sparsity in the U, V, and s matrices.
       See CUR decomposition method in contrast.
     * @param x is a 2-dimensional array of k vectors of length n in format
     *    double[n][k].  n is the number of samples, and k is the number of
     *    variables, a.k.a. dimensions.
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
        
        double eps = 1.e-15;
           
        /*
        NOTE that the Strang approach does not assume the data has been
        standardized to mean 0 and unit variance.
        
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
    
        //TODO: replace this covariance with the faster methods when those are ready for use
        //[x[0].length][x[0].length]
        double[][] cov = BruteForce.covariance(x);
        assert(nDimensions == cov.length);
        assert(nDimensions == cov[0].length);
                        
        //NOTE: we know that cov is symmetric positive definite, so there are
        //   many operations special to it one could explore for diagonalization
        
        SVD svd = SVD.factorize(new DenseMatrix(cov));
        // U is mxm
        // S is mxn
        // V is nxn
        
        // singular values are square roots of the eigenvalues:
        double[] s = svd.getS();
        
        int rank = 0;
        for (int i = 0; i < s.length; ++i) {
            if (Math.abs(s[i]) > eps) {
                rank++;
            }
        }
        
        if (nComponents > rank) {
            nComponents = rank;
        }
        
        //eigenvalue_i = lambda_i = (s_i^2)/(n-1)
        double[] eV = Arrays.copyOf(s, rank);
        for (int i = 0; i < eV.length; ++i) {
            eV[i] *= eV[i];
            eV[i] /= ((double)n - 1.);
        }
        
        // COLUMNS of u are the principal axes, a.k.a. principal directions
        // size is nDimensions x nDimensions
        DenseMatrix u = svd.getU();
        assert(nDimensions == u.numRows());
        assert(nDimensions == u.numColumns());
       
        // extract the nComponents columns of U as the principal axes, a.k.a. 
        //  principal directions.  array format: [nDimensions][nComponents]
        double[][] pa = new double[u.numRows()][nComponents];                                
        for (int row = 0; row < u.numRows(); ++row) {
            pa[row] = new double[nComponents];
            for (int col = 0; col < nComponents; ++col) {
                pa[row][col] = u.get(row, col);
            }
        }
        
        // NOTE: the principal scores is a vector Y_1 of length n:
        //     first principal score Y_1 = 
        //         sum of (first principal direction for each variable dot x[*][variable] )
        //  The magnitudes of the principal direction coefficients give the 
        //  contributions of each variable to that component. 
        //  Because the data have been standardized, they do not depend on the 
        //  variances of the corresponding variables.      
        
        // principal components are the projections of the principal axes on U
        //   the "1 sigma" lengths are:
        double[][] projections = new double[pa.length][pa[0].length];
        for (int row = 0; row < pa.length; ++row) {
            projections[row] = Arrays.copyOf(pa[row], pa[row].length);
            for (int col = 0; col < pa[row].length; ++col) {
                // note, if wanted "3 sigma" instead, use factor 3*s[col] here:
                projections[row][col] *= s[col];
            }
        }
        
        // extract the first nComponents of rows of V^T:
        double[][] v = Matrices.getArray(svd.getVt());
        v = MatrixUtil.copySubMatrix(v, 0, nComponents-1, 0, v[0].length-1);
        
        PCAStats stats = new PCAStats();
        stats.nComponents = nComponents;
        stats.principalDirections = pa;
        stats.eigenvalues = eV;
        stats.vTP = v;
        
        stats.s = new double[nComponents];
        for (int i = 0; i < nComponents; ++i) {
            stats.s[i] = s[i];
        }
                
        double total = 0;
        for (int j = 0; j < rank; ++j) {
            total += s[j];
        }
        double[] fracs = new double[rank];
        for (int j = 0; j < rank; ++j) {
            fracs[j] = s[j]/total;
        }
        double[] c = Arrays.copyOf(fracs, fracs.length);
        for (int j = 1; j < c.length; ++j) {
            c[j] += c[j-1];
        }
        stats.cumulativeProportion = c;
        
        double sum = 0;
        int p;
        for (int j = (nComponents+1); j <= nDimensions; ++j) {
            p = j - 1;
            sum += s[p];
        }
        stats.ssdP = sum;
        stats.fractionVariance = (total - sum)/total;
        
        System.out.println("singular values of SVD(cov) = sqrts of eigenvalues of cov = ");
        for (int i = 0; i < rank; ++i) {
            System.out.printf("%11.3e  ", s[i]);
        }
        System.out.println();
        System.out.println("eigenvalues of cov = ");
        for (int i = 0; i < eV.length; ++i) {
            System.out.printf("%11.3e  ", eV[i]);
        }
        System.out.println();
        System.out.println("U of SVD(cov) = ");
        for (int i = 0; i < u.numRows(); ++i) {
            for (int j = 0; j < u.numColumns(); ++j) {
                System.out.printf("%12.5e  ", u.get(i, j));
            }
            System.out.printf("\n");
        }
        System.out.println("V_p of SVD(cov) = ");
        for (int i = 0; i < v.length; ++i) {
            for (int j = 0; j < v[i].length; ++j) {
                System.out.printf("%12.5e  ", v[i][j]);
            }
            System.out.printf("\n");
        }
        
        System.out.println("eigenvalue fractions of total = ");
        for (int i = 0; i < fracs.length; ++i) {
            System.out.printf("%11.3e  ", fracs[i]);
        }
        System.out.println();
        System.out.println("eigenvalue cumulativeProportion= ");
        for (int i = 0; i < stats.cumulativeProportion.length; ++i) {
            System.out.printf("%11.3e  ", stats.cumulativeProportion[i]);
        }
        System.out.println();
        
        System.out.println("principal directions= ");
        for (int i = 0; i < stats.principalDirections.length; ++i) {
            for (int j = 0; j < stats.principalDirections[i].length; ++j) {
                System.out.printf("%12.5e  ", stats.principalDirections[i][j]);
            }
            System.out.printf("\n");
        }
        System.out.println("principal projections (=u_p*s)=");
        for (int i = 0; i < projections.length; ++i) {
            for (int j = 0; j < projections[i].length; ++j) {
                System.out.printf("%11.3e  ", projections[i][j]);
            }
            System.out.printf("\n");
        }
        
        System.out.flush();
        
        return stats;
    }
    
    /**
     * reconstruct an image from x, given principal directions.
     * from http://www.cs.toronto.edu/~jepson/csc420/
     * combined with the book by Strang "Introduction to Linear Algebra" 
     * and the book by Leskovec, Rajaraman, and Ullman "Mining of Massive Datasets"
     * <pre>
     *     ⃗r(⃗a_0) = m⃗_s + U_p * ⃗a 
     *         where U_p is the principalDirections matrix,
     *         where m_s is the sample x mean,
     *         where a_0 is
     *            a_0 = arg min_{a} ||⃗x − (m⃗_s + U _p * a)||^2
     * </pre>
     * @param x a 2-dimensional array of k vectors of length n in format
     *    double[n][k] which is double[nSamples][nDimensions] for which to apply
     *    the principal axes transformation on.
     * @param stats the statistics of the principal axes derived from
     * SVD of the covariance of the training data.
     * @return 
     */
    public static double[][] reconstruct(double[][] x, PCAStats stats) {
        
        /* find the 'a' which minimizes ||⃗x − (m⃗_s + p * a)||^2
        ⃗
        then the reconstruction is r(⃗a ) = m⃗ + U ⃗a 
        */
       
        double[][] b = MatrixUtil.multiply(x, stats.principalDirections);
        b = MatrixUtil.multiply(b, stats.vTP);
        
        return b;
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
         * the first nComponents principal directions of x. 
         * this is obtained from the first p columns of the U matrix of SVD(cov(x)).
         * the U_p columns have been transposed into rows here:
         * principalDirections[0] holds 1st principal direction,
         * principalDirections[1] holds 2nd principal direction, etc.
         * <pre>
         *    x * principalDirections * v^T_p gives the principal axes which
         *       minimize the variance.
         * </pre>
         */
        double[][] principalDirections;
        
        /**
         * the first p rows of the V^T matrix of SVD(cov(x)).
         * <pre>
         *    x * principalDirections * v^T_p gives the principal axes which
         *       minimize the variance.
         * </pre>
         */
        double[][] vTP;
        
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
         * the cumulative addition of 
         *     (sum_of_eigenvalues - eigenvalue_i)/sum_of_eigeneigenvalues
         */
        double[] cumulativeProportion;
        
        /**
         * Minimum residual variance
         */
        double ssdP;
        
        /**
         * the first nComponents singular values from the diagonal matrix of
         * the SVD of covariance of x;
         */
        double[] s;
        double[] eigenvalues;
    }
}
