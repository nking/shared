package algorithms.dimensionReduction;

import algorithms.correlation.BruteForce;
import algorithms.matrix.MatrixUtil;
import java.util.Arrays;

import algorithms.statistics.Covariance;
import algorithms.util.FormatArray;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 * find the principal components, that is, the minimal residual variance bases
 * of vectors of samples following G. Strang's SVD in machine learning in
 * http://www.cs.toronto.edu/~jepson/csc420/
 * 
 * from wikipedia:
 * The principal components of a collection of points in a real coordinate space 
 * are a sequence of p unit vectors, where the i-th vector is the direction of 
 * a line that best fits the data while being orthogonal to the first i-1 vectors. 
 * Here, a best-fitting line is defined as one that minimizes the average 
 * squared distance from the points to the line. These directions constitute an 
 * orthonormal basis in which different individual dimensions of the data are 
 * linearly uncorrelated. Principal component analysis (PCA) is the process of 
 * computing the principal components and using them to perform a change of basis 
 * on the data, sometimes using only the first few principal components and 
 * ignoring the rest.
 * 
 * from wikipedia:
 * 
 * PCA is closely related to Fisher's Discriminant Anaylsisi a.k.a. 
 * Linear Discriminant analysis and factor analysis in that they both look for 
 * linear combinations of variables which best explain the data.[4] 
 * LDA explicitly attempts to model the difference between the classes of data. 
 * PCA, in contrast, does not take into account any difference in class, 
 * and factor analysis builds the feature combinations based on differences 
 * rather than similarities.
 * NOTE:  Fisher's original article[1] actually describes a slightly different 
 * discriminant, which does not make some of the assumptions of LDA such as 
 * normally distributed classes or equal class covariances.
 * 
 * LDA is closely related to analysis of variance (ANOVA) and regression analysis, 
 * which also attempt to express one dependent variable as a linear combination 
 * of other features or measurements.[1][2] However, ANOVA uses categorical 
 * independent variables and a continuous dependent variable, whereas 
 * discriminant analysis has continuous independent variables and a categorical 
 * dependent variable (i.e. the class label).[3] Logistic regression and 
 * probit regression are more similar to LDA than ANOVA is, as they also 
 * explain a categorical variable by the values of continuous independent 
 * variables. These other methods are preferable in applications where it is 
 * not reasonable to assume that the independent variables are normally 
 * distributed, which is a fundamental assumption of the LDA method.
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
       http://www.mmds.org/mmds/v2.1/ch11-dimred.pdf
     * @param x is a 2-dimensional array of k vectors of length n in format
     *    double[n][k].  n is the number of samples, and k is the number of
     *    variables, a.k.a. dimensions.
     *          x should be zero-centered (mean=0).
     *          if x is a correlation matrix, it should be standardized to unit normalization (which is mean=0, stdev=1)
     * @param nComponents the number of principal components to return.
     * @return a few statistics of the SVD of the covariance of A, up to the
     * nComponents dimension.  Note that if the rank of the SVD(cov(A)) is
     * less than nComponents, then only the number of components as rank is returned.
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

           SVD( Sample Covariance ) = SVD((1/(n-1)) * X^T*X)

           SVD(A^T*A).U and V are both == SVD(A).V
           SVD(A*A^T).U and V are both == SVD(A).U

             A = U * D * V^T  where all U and V^T are from SVD(A), and D is from SVD(A).s
             A^T A = V * D^2 * V^T
             SampleCov = (1/(n-1)) * A^T A
                       = (1/(n-1)) * V * D^2 * V^T

             we also have U * D = X * V from SVD

             Jepson:  B first vector is the first column of V in SVD(C_S)
                where C_S is the sample covariance and SVD(C_S).V == SVD(C_S).U ~ SVD(X).V

         SVD(C_S) = SVD((1/(n-1)) * X^T*X)
         SVD(X^T*X).U == SVD(X).V
            V gives the principal axes or principal directions of the dat
            X*V = U * D = principal components = projection of the data onto the principal axes
               and the coords of the newly transformed data are on row_0 for the first data point, etc.
        */

        SVD svd = SVD.factorize(new DenseMatrix(x));
        // U is mxm
        // S is mxn
        // V is nxn
        
        double[] s = svd.getS();
        int i;
        int j;
        int rank = 0;
        for (i = 0; i < s.length; ++i) {
            if (Math.abs(s[i]) > eps) {
                rank++;
            }
        }
        
        if (nComponents > rank) {
            nComponents = rank;
        }

        double sumEvTrunc = 0;
        double[] eig = new double[rank];
        double sumEVTotal = 0;
        for (i = 0; i < eig.length; ++i) {
            eig[i] = (s[i] * s[i])/((double)n - 1.);
            sumEVTotal += eig[i];
            if (i < rank) {
                sumEvTrunc = sumEVTotal;
            }
        }
        //residual fractional eigenvalue is 1-(sum_{i=0 to nComponents-1}(s[i]) / sum_{i=0 to n-1}(s[i]))
        double residFracEigen = 1 - (sumEvTrunc/sumEVTotal);

        // COLUMNS of v are the principal axes, a.k.a. principal directions
        double[][] vT = MatrixUtil.convertToRowMajor(svd.getVt());
        double[][] pA = MatrixUtil.copySubMatrix(vT, 0, nComponents - 1, 0, vT[0].length - 1);

        // // X*V = U * diag(s) = principal components
        DenseMatrix u = svd.getU();
        // principal components for "1 sigma".  if need 2 sigma, multiply pc by 2,...
        double[][] pC = MatrixUtil.zeros(u.numRows(), nComponents);
        for (i = 0; i < u.numRows(); ++i) {
            for (j = 0; j < nComponents; ++j) {
                // row u_i * diag(s[j])
                pC[i][j] = u.get(i, j) * svd.getS()[j];
            }
        }

        PCAStats stats = new PCAStats();
        stats.nComponents = nComponents;
        stats.principalAxes = pA;
        stats.principalComponents = pC;
        stats.eigenValues = eig;

        double total = 0;
        for (j = 0; j < eig.length; ++j) {
            total += eig[j];
        }
        double[] fracs = new double[rank];
        for (j = 0; j < rank; ++j) {
            fracs[j] = eig[j]/total;
        }
        double[] c = Arrays.copyOf(fracs, fracs.length);
        for (j = 1; j < c.length; ++j) {
            c[j] += c[j-1];
        }
        stats.cumulativeProportion = c;
        
        double ssdP = 0;
        for (j = nComponents; j < nDimensions; ++j) {
            ssdP += eig[j];
        }
        stats.ssdP = ssdP;
        stats.fractionVariance = (total - ssdP)/total;

        System.out.printf("ssd_p=%.5e\n", stats.ssdP);
        System.out.printf("fractionVariance=%.5e\n", stats.fractionVariance);

        System.out.println("first few lines U = ");
        int end = (3 < u.numRows()) ? 3 : u.numRows();
        for (i = 0; i < end; ++i) {
            for (j = 0; j < u.numColumns(); ++j) {
                System.out.printf("%12.5e  ", u.get(i, j));
            }
            System.out.printf("\n");
        }
        System.out.println("VT = ");
        for (i = 0; i < vT.length; ++i) {
            for (j = 0; j < vT[i].length; ++j) {
                System.out.printf("%12.5e  ", vT[i][j]);
            }
            System.out.printf("\n");
        }
        System.out.printf("s fractions of total = \n%s\n",
                FormatArray.toString(fracs, "%.4e"));
        System.out.printf("eigenvalue cumulativeProportion=\n%s\n",
                FormatArray.toString(stats.cumulativeProportion, "%.4e"));
        System.out.printf("principal directions=\n%s\n", FormatArray.toString(pA, "%.4e"));
        System.out.printf("principal projections (=u_p*s)=\n%s\n", FormatArray.toString(pC, "%.4e"));
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

        //TODO: revisit this

        double[][] b = stats.principalComponents;
        
        return b;
    }
    
    public static class PCAStats {

        /** the number of components requested from the calculation, unless larger than the matrix
         * rank in which case the value will be the rank.
         */
        int nComponents;

        /**
         * principal axes a.k.a. principal directions of the data.
         * These are the first nComponents columns of the SVD V matrix.
         * When considering X as sample covariance, Note that SVD(A^T*A).U and V are both == SVD(A).V.
         */
        double[][] principalAxes;
        
        /**
         * principal components are the data projected onto the principal axes.
         * principal components = X*V = U * D where X = SVD(X).U * diag(SVD(X).s) * SVD(X).vT
         */
        double[][] principalComponents;
        
        /**
         * the fraction of the total variance Q_p (where p is the dimension number, that is, nComponents):
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
         *     (sum_of_eigenvalues - eigenvalue_i)/sum_of_eigenvalues
         *     after truncation.  that is, the sum of eigenvalues is the sum of the truncated number of eigenvalues,
         *      not the entire number from the original decomposition.
         */
        double[] cumulativeProportion;
        
        /**
         * Minimum residual variance
         */
        double ssdP;
        
        /**
         * the first nComponents number of eigenValues
         */
        double[] eigenValues;
    }
}
