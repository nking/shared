package algorithms.matrix;

import java.util.Arrays;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.NotConvergedException;

/**
 a class of utility methods for linear algebra.
 <pre>
 note, for solving linear least squares problems, a good resource aside
 * from Strang, Golub & van Loan, and Cormen et al. is:
 “Algorithms for Linear Least Squares Problems”, 
 Bj̈orck 1991, published in Computer Algorithms for Solving Linear Algebraic Equations;
 The State of the Art., Vol. 77 of NATO-ASI Series F: Computer and Systems Sciences, pages 57–92.
 </pre>
   LUPSolve, __, and __ follow
   pseudocode from Cormen, Leiserson, Rivest, and Stein, "Introduction to
   Computer Algorithms".
*/
public class LinearEquations {

    /**
     solving for 'x' in equations:
       a_0_0*x_0   +  a_0_1*x_1   + ... a_0_n*x_n     = b_0
       a_1_0*x_0   +  a_1_1*x_1   + ... a_1_n*x_n     = b_1
          ...
       a_n-1_0*x_0 +  a_n-1_1*x_1 + ... a_{n-1}_n*x_{n-1} = b_{n-1}
     where a is a square matrix, that is, there are n equations and n
     unknowns.
     LUP decomposition is more numerically stable than x = A^-1*b.
         can find L, U, and P such that P*A=L*U
     where L is a lower triangular matrix, U is an upper triangular matrix,
     and P is a permutation matrix.
     (a permutation matrix is all 0's excepting a single 1 in each column, uniquely.
     multiplying a vector by a permutation rearranges the members of the vector.)
     Uses forward substitution then back substitution.
       L*U*x=P*b.
         let y=U*x.
       L*y=P*b  (uses forward substitution with these results).
       U*x=y    (uses back substitution with these results).
       A*x=b    (then solves x_i=(y_i - summation_j=i_to_{n-1}(u_i_j*x_j))/u_i_i).
     runtime complexity is O(n^2) for backward and forward substitutions.
     @param ell is an nxn lower triangular matrix using row major format.
     @param u is an nxn upper triangular matrix using row major format.
     @param p is an array of length n holding permutation vector columns. 
            e.g. p_0=0 states that column 0 contains a 1 for row=0.
            e.g. p_1=2 states that column 2 contains a 1 for row=1.
            full_permutation_matrix = zeros everywhere except the row, col pairs in p
     @param b is an array of length n.
     @return array of x
    */
    public static double[] LUPSolve(double[][] ell, double[][] u, int[] p, double[] b) {
        int n = ell.length;
        assertSquareMatrix(ell, "ell");
        assertSquareMatrix(u, "u");
        assertArrayLength(n, p, "p");
        assertArrayLength(n, b, "b");

        double[] y = new double[n];
        double[] x = new double[n];
        double tempsum = 0;
        for (int i = 0; i < n; ++i) {
            tempsum = 0;
            for (int j = 0; j < i; j++) {
                tempsum += ell[i][j]*y[j];
            }
            y[i] = b[p[i]] - tempsum;
        }
        tempsum = 0;
        for (int i = n-1; i >= 0; i--) {
            tempsum = 0;
            if (i != (n-1)) {
                for (int j = i+1; j < n; j++) {
                    tempsum += u[i][j]*x[j];
                }
            }
            x[i] = (y[i] - tempsum)/u[i][i];
        }
        return x;
    }
    
    /**
     * solve for x in A*x=b by LU decomposition
     * @param a
     * @param b
     * @return 
     */
    public static double[] solveXFromLUDecomposition(double[][] a, double[] b) {
        
        int n = a.length;
        
        DenseVector x = new DenseVector(n);
        DenseVector _b = new DenseVector(b);
        
        DenseMatrix aa = new DenseMatrix(a);
        x = (DenseVector)aa.solve(_b, x);
                        
        return x.getData();
    }

    /**
     * an efficient LUP decomposition for a being a square non-singular matrix and
     * P is the Identity matrix.   uses Gaussian elimination and the Schur
     * complement while making recursive subdivision subdivisions.
     * The runtime is O(n^3).
     * @param a two dimensional array in row major format.  
     * a is a non-singular matrix(i.e. has exactly one solution).  the rank of
     * a is n (it's dimensions are m x n).
     * @return LU a wrapper holding the 2 two-dimensional row major output arrays.
     * L and U.  they are both size nXn where n=a.length
     */
    public static LU LUDecomposition(double[][] a) {
        int n = a.length;
        assertSquareMatrix(a, "a");
        
        a = copy(a);
        
        //nXn for both:
        double[][] ell = new double[n][];
        double[][] u = new double[n][];
        for (int i = 0; i < n; ++i) {
            ell[i] = new double[n];
            u[i] = new double[n];
            ell[i][i] = 1;
            u[i][i] = 1;
        }
        
        for (int k = 0; k < n; k++) {
            u[k][k] = a[k][k];
            
            for (int i = k+1; i < n; i++) {
                ell[i][k] = a[i][k]/u[k][k];
                u[k][i] = a[k][i];
            }
            
            for (int i = k+1; i < n; i++) {
                for (int j = k+1; j < n; j++) {
                    a[i][j] -= (ell[i][k] * u[k][j]); 
                }
            }
        }
                
        return new LU(ell, u);
    }
    
    /**
     * for a symmetric matrix A (with real numbers in matrix of size nXn), 
     * perform and L-D-M decomposition which is a 
     * variation of L-U decomposition.  If all leading principal sub-matrices of 
     * A are non-singular then there exist lower unit triangular matrices
     * L and M and a diagonal matrix D = diag(d1, d2, ...d_n) where d_i = u_i_i.
     * A = L*D*M^T.
     * Note that D is non-singular and that M^T = D^-1*U is unit upper triangular.
     * (A=L*U = L*D*(D^-1*U = L*D*M^T).
     * <pre>
     * Once A = L*D*M^T is decomposed, one can solve A*x = b:
     * L*y = b   (in (n^2)/2 flops
     * D*z= y    (in n flops)
     * M^T*x = z (in (n^2)/2 flops)
     * 
     * References:
     * Golub and va Loan, "MAtrix Computations", Section 5.1
     * </pre>
     * 2n^3/3 flops.
     * @param a two dimensional array in row major format.  
     * a is a symmetric matrix with dimensions n x n.
     * @return LDM a wrapper holding the 2 two-dimensional row major output arrays.
     * L and M and the diagonal matrix D as a an array of the diagonal.
     */
    public static LDM LDMDecomposition(double[][] a) {
        
        int n = a.length;
        
        assertSquareMatrix(a, "a");
        
        // need to have a in row echelon reduced state for this:
        //assertSymmetrix(a);
        
        double[][] ell = MatrixUtil.zeros(n, n);
        double[][] m = MatrixUtil.zeros(n, n);
        double[] d = new double[n];
        
        int k, i, p;
        for (k = 1; k <= n; ++k) {
            d[k-1] = a[k-1][k-1];
            for (p = 1; p <= (k-1); ++p) {
                d[k-1] -= (ell[k-1][p-1]*d[p-1]*m[k-1][p-1]);
            }
            // error in the 1st edition of Matrix Computations; i=k to n, but they use i=k+1 to n
            for (i = k/*k+1*/; i <= n; ++i) {
                
                ell[i-1][k-1] = a[i-1][k-1];
                for (p = 1; p <= (k-1); ++p) {
                    ell[i-1][k-1] -= (ell[i-1][p-1]*d[p-1]*m[k-1][p-1]);
                }
                ell[i-1][k-1] /= d[k-1];
                
                m[i-1][k-1] = a[k-1][i-1];
                for (p = 1; p <= (k-1); ++p) {
                    m[i-1][k-1] -= (ell[k-1][p-1]*d[p-1]*m[i-1][p-1]);
                }
                m[i-1][k-1] /= d[k-1];
            }
        }
        
        LDM ldm = new LDM(ell, d, m);
        
        return ldm;
    }
    
    /**
     * compute the cholesky decomposition for symmetric positive definite matrix
     * a.
     * reference:
     * golub & van loan "matrix computations, theorem 5.2-3.
     * This method uses LDL decomposition to compute G in 
     * a = G*G^T where G is a lower triangular matrix with positive
    * diagonal entries.
    * TODO: consider implementing algorithm 5.2-1 also.
    * And consider re-compiling MTJ library to include the DenseCholesky.java
    * which seems to be unsupported in the jar in this project.  alternatively
    * could invoke the LAPACK method directly using the current build. 
    * see https://github.com/fommil/matrix-toolkits-java/blob/master/src/main/java/no/uib/cipr/matrix/DenseCholesky.java)
     * @param a symmetric positive definite matrix 
     * @param eps value for an error tolerance around zero used in the LDL decomposition.
     * @return lower triangular matrix G  which G is a lower triangular matrix with positive
    * diagonal entries.  a = G*G^T.
    */
    public static double[][] choleskyDecompositionViaLDL(double[][] a, double eps) {
        
        int n = a.length;
        
        assertSquareMatrix(a, "a");
        
        //assertPositiveDefinite(a, "a");
        
        LDL ldl = LDLDecomposition(a, eps);
        
        double[] d = Arrays.copyOf(ldl.getD(), ldl.getD().length);
        for (int i = 0; i < d.length; ++i) {
            d[i] = Math.sqrt(d[i]);
        }
        double[][] g = MatrixUtil.multiplyByDiagonal(ldl.getL(), d);
        
        return g;
    }
       
    /**
     * for a nonsingular symmetric matrix A (with real numbers in matrix of size nXn), 
     * perform and L-D-L decomposition which is a 
     * variation of L-U decomposition.  Computes a unit lower triangular matrix
     * L and a diagonal matrix D = diag(d1, d2, ...d_n) such that
     * A = L*D*L^T.
     * <pre>
     * 
     * References:
     * Golub and va Loan, "Matrix Computations", Algorithm 5.1.2
     * </pre>
     * n^3/6 flops.
     * @param a two dimensional array in row major format.  
     * a is a symmetric matrix with dimensions n x n.
     * @return LDM a wrapper holding the 2 two-dimensional row major output arrays.
     * L and M and the diagonal matrix D as a an array of the diagonal.
     */
    public static LDL LDLDecomposition(double[][] a) {
        double eps = 1.e-7;
        return LDLDecomposition(a, eps);
    }
        
    /**
     * for a nonsingular symmetric matrix A (with real numbers in matrix of size nXn), 
     * perform and L-D-L decomposition which is a 
     * variation of L-U decomposition.  Computes a unit lower triangular matrix
     * L and a diagonal matrix D = diag(d1, d2, ...d_n) such that
     * A = L*D*L^T.
     * <pre>
     * 
     * References:
     * Golub and va Loan, "Matrix Computations", Algorithm 5.1.2
     * </pre>
     * n^3/6 flops.
     * @param a two dimensional array in row major format.  
     * a is a symmetric matrix with dimensions n x n.
     * @param eps value for a tolerance of an error around 0.
     * @return LDM a wrapper holding the 2 two-dimensional row major output arrays.
     * L and M and the diagonal matrix D as a an array of the diagonal.
     */
    public static LDL LDLDecomposition(double[][] a, double eps) {
        
        int n = a.length;
        
        assertSquareMatrix(a, "a");
                
        // need to have a in row echelon reduced state for this:
        //assertSymmetrix(a);
        
        //NOTE: if need to conserve memory, can remove the copy statement,
        //    and add to javadoc comments that a is modified in place.
        a = MatrixUtil.copy(a);
        double[] r = new double[n];
        double[] d = new double[n];
        
        int k, i, p;
        for (k = 1; k <= n; ++k) {
            for (p = 1; p <= (k-1); ++p) {
                r[p-1] = d[p-1]*a[k-1][p-1];
            }
            d[k-1] = a[k-1][k-1];
            for (p = 1; p <= (k-1); ++p) {
                d[k-1] -= a[k-1][p-1]*r[p-1];
            }
            if (Math.abs(d[k-1]) < eps) {
                System.err.printf("Error: number %.3e is smaller than %.3e\n", 
                    Math.abs(d[k-1]), eps);
                return null;
            }
            for (i = k; i <= n; ++i) {
                for (p = 1; p <= (k-1); ++p) {
                    a[i-1][k-1] -= a[i-1][p-1]*r[p-1];
                }
                a[i-1][k-1] /= d[k-1];
            }
        }
        
        for (i = 0; i < n; ++i) {
            for (k = i+1; k < n; ++k) {
                a[i][k] = 0;
            }
        }
        
        LDL ld = new LDL(a, d);
        
        return ld;
    }
    
    
    
    /**
     * an LUP decomposition for a being a square non-singular matrix that tries 
     * to reduce errors due to division by small numbers.  
     * creates a permutation matrix to pivot rows so that the row reduction
     * divisions are by the largest numbers.
     * uses Gaussian elimination and the Schur
     * complement while making recursive subdivision subdivisions.
     * The runtime is O(n^3).
     * @param a two dimensional array in row major format.  
     * a is a non-singular matrix(i.e. has exactly one solution).  the rank of
     * a is n (it's dimensions are m x n).
     * @return LUP a wrapper holding the 2 two-dimensional row major output arrays.
     * L and U and the condensed permutation array p, where P*A=L*U.
     */
    public static LUP LUPDecomposition(double[][] a) {
        int n = a.length;
        assertSquareMatrix(a, "a");
        
        a = copy(a);
        
        int[] pi = new int[n];
        for (int i = 0; i < n; ++i) {            
            pi[i] = i;
        }
        
        double p;
        double temp;
        int k2 = -1;
        double swap;
        int swapI;
        for (int k = 0; k < n; k++) {
            p = 0;
            for (int i = k; i < n; i++) {
                temp = Math.abs(a[i][k]);
                if (temp > p) {
                    p = temp;
                    k2 = i;
                }
            }
            if (p == 0.) {
                throw new IllegalStateException("Error: a is a singular matrix");
            }
            assert(k2 >= 0);
            swapI = pi[k2];
            pi[k2] = pi[k];
            pi[k] = swapI;
            for (int i = 0; i < n; i++) {
                swap = a[k2][i];
                a[k2][i] = a[k][i];
                a[k][i] = swap;
            }
            for (int i = k+1; i < n; i++) {
                a[i][k] /= a[k][k];
                for (int j = k+1; j < n; j++) {
                    a[i][j] -= (a[i][k] * a[k][j]);
                }
            }
        }
       
        // setting the upper and lower triangles of a into ell and u to pass
        //   pack as arguments.
        //   the original pseudocode uses only a and modifies the original in place
        //     to conserve space.
        double[][] ell = new double[n][];
        double[][] u = new double[n][];
        for (int i = 0; i < n; ++i) {
            ell[i] = new double[n];
            u[i] = new double[n];
            u[i][i] = 1;
            ell[i][i] = 1;
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i > j) {
                    ell[i][j] = a[i][j];
                } else {
                    u[i][j] = a[i][j];
                }
            }
        }
        
        return new LUP(ell, u, pi);
    }
    
    /**
     * given data points xy and the assumption that measurement errors are small,
     * fit a polynomial of order polyOrder to the data points, minimizing the
     * error, i.e. solve for coefficients c in y_i = summation(c_i*x^i) + error.
     * calculated by c = pseudo-inverse of A * y where A is the components of
     * x as polynomial factors.
     * NOTE that a regularized linear least squares algorithm called Elastic-Net
     * is implemented as thirdparty.scipy.optimization.ElastticNet.
     * This method follows pseudocode in chapter 28 of Cormen et al. Introduction
     * To Algorithms.
     * @param xy two dimensional array of format row0=[x0,y0], row1=[x1,y1], etc.
     * @param polyOrder the order of a polynomial to fit.  should be .lte. the
     * number of rows.
     * @param solveForFullRank
     * when 'True' AX=b has no solution (e.g. xy.length is larger than xy[0].length)
     * and the algorithm uses (inverse(A^T*A) * A^T) for the pseudo-inverse
     * (see Chap 4.3 from the book "Introduction to Linear
     * Algebra" by W Gilbert Strang and Chap 28 from the book "Introductionvto 
     * Algorithms" by Cormen, Leiserson, Rivest, and Stein.)
     * when solveForFullRank is set to 'False' this method uses the SVD to
     * create a pseudoinverse (see Chap 7 from the book "Introduction to Linear
     * Algebra" by W Gilbert Strang.)
     *  
     * @return coefficients c where y_i = summation(c_i*x^i) + error
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static double[] leastSquaresPolynomial(double[][] xy, int polyOrder, boolean
        solveForFullRank) throws NotConvergedException {
        return _leastSquaresPolynomial(xy, polyOrder, solveForFullRank);
    }
    
    /**
     * given data points xy and the assumption that measurement errors are small,
     * fit a polynomial of order polyOrder to the data points, minimizing the
     * error, i.e. solve for coefficients c in y_i = summation(c_i*x^i) + error.
     * calculated by c = pseudo-inverse of A * y where A is the components of
     * x as polynomial factors.
     * NOTE that a regularized linear least squares algorithm called Elastic-Net
     * is implemented as thirdparty.scipy.optimization.ElastticNet.
     * This method follows pseudocode in chapter 28 of Cormen et al. Introduction
     * To Algorithms.
     * @param xy two dimensional array of format row0=[x0,y0], row1=[x1,y1], etc.
     * @param polyOrder the order of a polynomial to fit.  should be .lte. the
     * number of rows.
     * @return coefficients c where y_i = summation(c_i*x^i) + error
     */
    public static double[] leastSquaresPolynomial(double[][] xy, int polyOrder) throws NotConvergedException {
        return _leastSquaresPolynomial(xy, polyOrder, false);
    }
    
    /**
     * given data points xy and the assumption that measurement errors are small,
     * fit a polynomial of order polyOrder to the data points, minimizing the
     * error, i.e. solve for coefficients c in y_i = summation(c_i*x^i) + error.
     * calculated by c = pseudo-inverse of A * y where A is the components of
     * x as polynomial factors.
     * NOTE that a regularized linear least squares algorithm called Elastic-Net
     * is implemented as thirdparty.scipy.optimization.ElastticNet.
     * This method follows pseudocode in chapter 28 of Cormen et al. Introduction
     * To Algorithms.
     * @param xy two dimensional array of format row0=[x0,y0], row1=[x1,y1], etc.
     * @param polyOrder the order of a polynomial to fit.  should be .lte. the
     * number of rows.
     * @param solveFullRank
     * when 'True' AX=b has no solution (e.g. xy.length is larger than xy[0].length)
     * and the algorithm uses (inverse(A^T*A) * A^T) for the pseudo-inverse
     * (see Chap 4.3 from the book "Introduction to Linear
     * Algebra" by W Gilbert Strang and Chap 28 from the book "Introductionvto 
     * Algorithms" by Cormen, Leiserson, Rivest, and Stein.)
     * when solveForFullRank is set to 'False' this method uses the SVD to
     * create a pseudoinverse (see Chap 7 from the book "Introduction to Linear
     * Algebra" by W Gilbert Strang.)
     * @return coefficients c where y_i = summation(c_i*x^i) + error
     */
    public static double[] _leastSquaresPolynomial(double[][] xy, int polyOrder,
        boolean solveFullRank) throws NotConvergedException {
        
        if (polyOrder < 0) {
            throw new IllegalArgumentException("polyOrder must be 0 or larger");
        }
        int nRows = xy.length;
        
        if (polyOrder > nRows) {
            throw new IllegalArgumentException("polyOrder should be smaller "
                    + " than or equal to the number of data rows");
        }
        
        // create matrix A
        double[][] a = new double[nRows][];
        double x;
        for (int i = 0; i < nRows; ++i) {
            a[i] = new double[polyOrder + 1];
            a[i][0] = 1;
        }
        for (int i = 0; i < nRows; ++i) {
            x = xy[i][0];
            for (int k = 0; k < polyOrder; ++k) {
                a[i][k + 1] = x * a[i][k];
            }
        }
        
        double[][] aPInv;
        if (solveFullRank) {
            //AX=b has no solution
            //using A_pseudoinverse = inverse(A^T*A) * A^T
            aPInv = MatrixUtil.pseudoinverseFullRank(a);
        } else {
            // uses SVD of a in pseudo-inverse
            // Note that if A^-1 exists, then the pseudo-inverse of A is equal to the
            //inverse of A.
            // method is valid for  rows .lte. rank,
            // but for rows .gt. rank, the null-space needs to be solved instead.
            aPInv = MatrixUtil.pseudoinverseRankDeficient(a);
        }
        
        double[] y = new double[nRows];
        for (int i = 0; i < nRows; ++i) {
            y[i] = xy[i][1];
        }
        
        double[] c = MatrixUtil.multiplyMatrixByColumnVector(aPInv, y);
        
        return c;
    }
    
    public static class LU {
        double[][] ell;
        double[][] u;
        public LU(double[][] ell, double[][] u) {
            this.ell = ell;
            this.u = u;
        }
    }
    
    public static class LUP {
        double[][] ell;
        double[][] u;
        int[] p;
        public LUP(double[][] ell, double[][] u, int[] p) {
            this.ell = ell;
            this.u = u;
            this.p = p;
        }
    }
    
    /**
     * lower triangular portion of A = L * D * M^T
     */
    public static class LDM {
        double[][] ell;
        double[] d;
        double[][] m;
        public LDM(double[][] ell, double[] d, double[][] m) {
            this.ell = ell;
            this.d = d;
            this.m = m;
        }
    }
    
    public static class LDL {
        private double[][] ell;
        private double[] d;
        public LDL(double[][] ell, double[] d) {
            this.ell = ell;
            this.d = d;
        }

        /**
         * @return the L
         */
        public double[][] getL() {
            return ell;
        }

        /**
         * @return the d
         */
        public double[] getD() {
            return d;
        }
    }
    
    private static void assertSquareMatrix(double[][] a, String name) {
        int n = a.length;
        if (a[0].length != n) {
            throw new IllegalArgumentException(name + " must be a square matrix");
        }
    }
    
    private static void assertPositiveDefinite(double[][] a, String name) {
        int n = a.length;
        if (!MatrixUtil.isPositiveDefinite(a)) {
            throw new IllegalArgumentException(name + " must be a positive definite matrix");
        }
    }

    private static void assertArrayLength(int n, double[] a, String name) {
        if (a.length != n) {
            throw new IllegalArgumentException(name + " must be length " + n);
        }
    }

    private static void assertArrayLength(int n, int[] a, String name) {
        if (a.length != n) {
            throw new IllegalArgumentException(name + " must be length " + n);
        }
    }
    
    private static double[][] copy(double[][] a) {
        int nrows = a.length;
        double[][] c = new double[nrows][];
        for (int i = 0; i < nrows; ++i) {
            c[i] = Arrays.copyOf(a[i], a[i].length);
        }
        return c;
    }
} 
