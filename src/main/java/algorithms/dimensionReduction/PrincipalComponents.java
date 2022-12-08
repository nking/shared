package algorithms.dimensionReduction;

import algorithms.correlation.BruteForce;
import algorithms.matrix.MatrixUtil;
import java.util.Arrays;

import algorithms.statistics.Covariance;
import algorithms.statistics.Standardization;
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
     * using Singular Value Decomposition.
     * NOTE: the data need to be zero centered first.
     *
     <pre>
     the method follows:
     http://www.cs.toronto.edu/~jepson/csc420/
     combined with the book by Strang "Introduction to Linear Algebra"
     also useful:
     https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
     and https://online.stat.psu.edu/stat505/book/export/html/670 for testing results,
     </pre>
     <pre>
     NOTE:
     variance-covariance(standardized data) == correlation(unstandardized data) == correlation(zero mean centered data).
     therefore, pca using the standardized data == pca using the correlation matrix.
     Also, the eigen of cov(standardized) == eigen of cor(unstandardized).
     </pre>
     @param x is a 2-dimensional array of k vectors of length n in format
     *    double[n][k].  n is the number of samples, and k is the number of
     *    variables, a.k.a. dimensions.
     *    x should be zero-centered (mean=0) OR standardized to unit normalization (which is mean=0, stdev=1).
     *    If the variance and scale of the variables are different, then unit standard normalization should be used.
     @param nComponents the number of principal components to return.
     @return the principal axes, the principal components, and
     a few statistics of the CUR decomposition of A, up to the
     nComponents dimension.  Note that if the rank of the A is
     less than nComponents, then only the number of components as rank is returned.
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static PCAStats calcPrincipalComponents(double[][] x, int nComponents) throws NotConvergedException {
        return calcPrincipalComponents(x, nComponents, false);
    }

    /**
     * calculate the principal components of the unit standardized data x
     * using Singular Value Decomposition or CUR Decomposition.
     * NOTE: the data need to be zero centered first.
     *
     <pre>
     the method follows:
     http://www.cs.toronto.edu/~jepson/csc420/
     combined with the book by Strang "Introduction to Linear Algebra"
     and the book by Leskovec, Rajaraman, and Ullman "Mining of Massive Datasets".
     also useful:
     https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
     and https://online.stat.psu.edu/stat505/book/export/html/670
     </pre>
     <pre>
     NOTE:
     variance-covariance(standardized data) == correlation(unstandardized data) == correlation(zero mean centered data).
     therefore, pca using the standardized data == pca using the correlation matrix.
     Also, the eigen of cov(standardized) == eigen of cor(unstandardized).
     </pre>
     @param x is a 2-dimensional array of k vectors of length n in format
     *    double[n][k].  n is the number of samples, and k is the number of
     *    variables, a.k.a. dimensions.
     *    x should be zero-centered (mean=0) OR standardized to unit normalization (which is mean=0, stdev=1).
     *    If the variance and scale of the variables are different, then unit standard normalization should be used.
     @param useEVProp if true, uses the cumulative fraction of eigenvalues in combination with the
     proportion prop to determine the number of principal components to calculate, else if false,
     uses the cumulative fraction of singular values in combination with the proportaion prop to
     determine the number of principal components.  Note that the closest cumulative fraction to prop is used
     rather than meeting and possibly exceeding prop.
     e.g.  prop of 80% to 90% w/ useEVProp = true.
     @param prop
     @return the principal axes, the principal components, and
     a few statistics of the CUR decomposition of A, up to the
     nComponents dimension.  Note that if the rank of the A is
     less than nComponents, then only the number of components as rank is returned.
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static PCAStats calcPrincipalComponents(double[][] x, boolean useEVProp,
                                                   double prop) throws NotConvergedException {

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
        int n = x.length;
        int nDimensions = x[0].length;

        double eps = 1.e-15;

        // U is mxm
        // S is mxn
        // V is nxn

        double[][] u;
        double[] s;
        double[][] vT;
        SVD svd = SVD.factorize(new DenseMatrix(x));
        s = svd.getS();
        vT = MatrixUtil.convertToRowMajor(svd.getVt());
        u = MatrixUtil.convertToRowMajor(svd.getU());

        int i;
        int j;
        int rank = 0;
        for (i = 0; i < s.length; ++i) {
            if (Math.abs(s[i]) > eps) {
                rank++;
            }
        }

        double[] eig = new double[rank];
        double sumEVTotal = 0;
        double sumSTotal = 0;
        for (i = 0; i < eig.length; ++i) {
            eig[i] = (s[i] * s[i])/((double)n - 1.);
            sumEVTotal += eig[i];
            sumSTotal += s[i];
        }

        double[] fracs = new double[rank];
        if (useEVProp) {
            for (j = 0; j < fracs.length; ++j) {
                fracs[j] = eig[j]/sumEVTotal;
            }
        } else {
            for (j = 0; j < fracs.length; ++j) {
                fracs[j] = s[j]/sumSTotal;
            }
        }

        // find smallest diff between cumulative fraction and prop

        double[] c = Arrays.copyOf(fracs, fracs.length);
        int idxM = 0;
        double diff = Math.abs(c[0] - prop);
        double min = diff;
        for (j = 1; j < c.length; ++j) {
            c[j] += c[j-1];
            diff = Math.abs(c[j] - prop);
            if (diff < min) {
                idxM = j;
                min = diff;
            }
        }
        int nComponents = idxM + 1;

        System.out.printf("nComponents=%d\n", nComponents);
        System.out.printf("eig = %s\n", FormatArray.toString(eig, "%.4e"));
        System.out.printf("s = %s\n", FormatArray.toString(s, "%.4e"));
        System.out.printf("fractions (eig or s) of total = %s\n", FormatArray.toString(fracs, "%.4e"));
        System.out.printf("cumulative proportion=\n%s\n", FormatArray.toString(c, "%.4e"));

        return calcPrincipalComponents(x, nComponents, u, s, vT, eps);
    }

    /**
     * calculate the principal components of the unit standardized data x
     * using Singular Value Decomposition or CUR Decomposition.
     * NOTE: the data need to be zero centered first.
     *
     <pre>
     the method follows:
     http://www.cs.toronto.edu/~jepson/csc420/
     combined with the book by Strang "Introduction to Linear Algebra"
     and the book by Leskovec, Rajaraman, and Ullman "Mining of Massive Datasets".
     also useful:
     https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
     and https://online.stat.psu.edu/stat505/book/export/html/670
     </pre>
     <pre>
     NOTE:
     variance-covariance(standardized data) == correlation(unstandardized data) == correlation(zero mean centered data).
     therefore, pca using the standardized data == pca using the correlation matrix.
     Also, the eigen of cov(standardized) == eigen of cor(unstandardized).
     </pre>
     @param x is a 2-dimensional array of k vectors of length n in format
     *    double[n][k].  n is the number of samples, and k is the number of
     *    variables, a.k.a. dimensions.
     *    x should be zero-centered (mean=0) OR standardized to unit normalization (which is mean=0, stdev=1).
     *    If the variance and scale of the variables are different, then unit standard normalization should be used.
     @param nComponents the number of principal components to return.
     @param useCUR if true, uses CUR decomposition instead of Singular Value Decomposition.  CUR decomposition
     is useful for very large datasets.
     @return the principal axes, the principal components, and
     a few statistics of the CUR decomposition of A, up to the
     nComponents dimension.  Note that if the rank of the A is
     less than nComponents, then only the number of components as rank is returned.
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static PCAStats calcPrincipalComponents(double[][] x, int nComponents, boolean useCUR) throws NotConvergedException {

        int n = x.length;
        int nDimensions = x[0].length;

        if (nComponents > nDimensions) {
            throw new IllegalArgumentException("x has " + nDimensions + " but "
                    + " nComponents principal components were requested.");
        }

        double eps = 1.e-15;

        // U is mxm
        // S is mxn
        // V is nxn

        double[][] u;
        double[] s;
        double[][] vT;
        if (!useCUR) {
            SVD svd = SVD.factorize(new DenseMatrix(x));
            s = svd.getS();
            vT = MatrixUtil.convertToRowMajor(svd.getVt());
            u = MatrixUtil.convertToRowMajor(svd.getU());
        } else {
            CURDecomposition.CUR cur = CURDecomposition.calculateDecomposition(x, nComponents);
            MatrixUtil.SVDProducts svd = cur.getApproximateSVD();
            s = svd.s;
            vT = svd.vT;
            u = svd.u;
        }

        return calcPrincipalComponents(x, nComponents, u, s, vT, eps);
    }

    private static PCAStats calcPrincipalComponents(double[][] x, int nComponents,
        double[][] u, double[] s, double[][] vT, double eps) {

        int n = x.length;
        int nDimensions = x[0].length;

        if (nComponents > nDimensions) {
            throw new IllegalArgumentException("x has " + nDimensions + " but "
                    + " nComponents principal components were requested.");
        }

        // U is mxm
        // S is mxn
        // V is nxn

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

        // re-doing these for output stats

        double sumEvTrunc = 0;
        double[] eig = new double[rank];
        double sumEVTotal = 0;
        double sumSTotal = 0;
        for (i = 0; i < eig.length; ++i) {
            eig[i] = (s[i] * s[i])/((double)n - 1.);
            sumEVTotal += eig[i];
            sumSTotal += s[i];
            if (i < nComponents) {
                sumEvTrunc = sumEVTotal;
            }
        }

        double[] eigFracs = new double[rank];
        double[] sFracs = new double[rank];
        for (j = 0; j < eig.length; ++j) {
            eigFracs[j] = eig[j]/sumEVTotal;
            sFracs[j] = s[j]/sumSTotal;
        }
        //residual fractional eigenvalue is 1-(sum_{i=0 to nComponents-1}(s[i]) / sum_{i=0 to n-1}(s[i]))
        double residFracEigen = 1 - (sumEvTrunc/sumEVTotal);

        double[] cEig = Arrays.copyOf(eigFracs, eigFracs.length);
        for (j = 1; j < cEig.length; ++j) {
            cEig[j] += cEig[j-1];
        }
        double[] cSingular = Arrays.copyOf(sFracs, sFracs.length);
        for (j = 1; j < cSingular.length; ++j) {
            cSingular[j] += cSingular[j-1];
        }

        // COLUMNS of v are the principal axes, a.k.a. principal directions

        double[][] pA = MatrixUtil.copySubMatrix(vT, 0, nComponents - 1, 0, vT[0].length - 1);

        // X*V = U * diag(s) = principal components
        double[][] uP = MatrixUtil.zeros(u.length, nComponents);
        // principal components for "1 sigma".  if need 2 sigma, multiply pc by 2,...
        double[][] pC = MatrixUtil.zeros(u.length, nComponents);
        for (i = 0; i < u.length; ++i) {
            for (j = 0; j < nComponents; ++j) {
                // row u_i * diag(s[j])
                pC[i][j] = u[i][j] * s[j];
                uP[i][j] = u[i][j];
            }
        }
        //checked, same as pC = U * diag(s)
        //double[][] xv = MatrixUtil.multiply(x, MatrixUtil.transpose(vT));
        //System.out.printf("xV=\n%s\n", FormatArray.toString(xv, "%.4e"));

        PCAStats stats = new PCAStats();
        stats.nComponents = nComponents;
        stats.principalAxes = pA;
        stats.principalComponents = pC;
        stats.eigenValues = eig;
        stats.uP = uP;
        stats._s = Arrays.copyOf(s, s.length);
        stats.cumulativeProportion = cEig;

        double ssdP = 0;
        for (j = nComponents; j < eig.length; ++j) {
            ssdP += eig[j];
        }
        stats.ssdP = ssdP;
        stats.fractionVariance = (sumEVTotal - ssdP)/sumEVTotal;

        System.out.printf("ssd_p=%.5e\n", stats.ssdP);
        System.out.printf("fractionVariance=%.5e\n", stats.fractionVariance);
        System.out.printf("principal axes=\n%s\n", FormatArray.toString(pA, "%.4e"));
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
     *
     @param m the means of each column of x (used in the zero-correction of x).  the length is k
     *          where k is the number of columns in x before
     *          dimension reduction.
     @param stats the statistics of the principal axes derived from
     * SVD of the covariance of the training data.
     @return the reconstruction of x.  size is [n x k]
     */
    public static double[][] reconstruct(double[] m, PCAStats stats) {

        int p = stats.nComponents;
        int k = stats._s.length;
        int n = stats.principalComponents.length;

        if (n != stats.uP.length) {
            throw new IllegalArgumentException("n must be stats.principalComponents.length which must be equal to " +
                    "stats.principalComponents.length");
        }
        if (k != m.length) {
            throw new IllegalArgumentException("m.length must equal k == stats.principalAxes.length");
        }

        //U_p is stats.uP

        //find the 'a' which minimizes ||⃗x − (m⃗_s + p * a)||^2
        //then the reconstruction is r(⃗a ) = m⃗ + U ⃗a

        // U_p is [n X p]
        // x_j is [n X 1]
        // a_j is [p X 1]
        // A is [p X k]

        // reconstruction of x0 is
        //   U_p * a_0
        //   [n X p] * [p X 1] = [n X 1];

        int i;
        int j;

        double[][] re = MatrixUtil.multiply(stats.principalComponents,
                stats.principalAxes);
        //reconstructions = stats.principalComponents dot stats.pA + column means
        //    which is u_P * s dot pA
        //[n X p] * [p X k] = [n X k]
        for (i = 0; i < re.length; ++i) {
            for (j = 0; j < re[i].length; ++j) {
                re[i][j] += m[j];
            }
        }
        return re;
    }
    
    /**
     *
     */
    public static class PCAStats {

        /** the number of components requested from the calculation, unless larger than the matrix
         * rank in which case the value will be the rank.
         */
        int nComponents;

        /**
         * principal axes a.k.a. principal directions of the data.
         * These are the first nComponents columns of the SVD(X).V matrix.
         * When considering X as sample covariance, Note that SVD(A^T*A).U and V are both == SVD(A).V.
         * the size is [k X k] where k is x[0].length.
         */
        double[][] principalAxes;

        /**
         * principal components are the data projected onto the principal axes.
         * principal components = X*V = U * D where X = SVD(X).U * diag(SVD(X).s) * SVD(X).vT
         * The size is [n X p] where n is x.length and p is the number of components.
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

        /**
         * the singular values of SVD(A)
         */
        double[] _s;

        /**
         * the first p eigenvectors of the left eigenvector matrix
         * (which is p columns of U from A = SVD(A).U * SVD(A).D * SVD(A).V^T
         */
        double[][] uP;

    }
}
