package algorithms.matrix;

import algorithms.matrix.LinearEquations.LU;
import algorithms.matrix.LinearEquations.LUP;
import algorithms.misc.Misc0;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import algorithms.util.PairInt;
import gnu.trove.list.array.TDoubleArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;

import gnu.trove.set.TDoubleSet;
import gnu.trove.set.hash.TDoubleHashSet;
import no.uib.cipr.matrix.*;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectDoubleIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import no.uib.cipr.matrix.sparse.ArpackSym;
import no.uib.cipr.matrix.sparse.FlexCompColMatrix;
import no.uib.cipr.matrix.sparse.LinkedSparseMatrix;

/**
 
 <pre>
   misc notes:
 
   The eigenvalues can be determined in a few ways depending upon the matrix:
       If A is a positive definite matrix, can use the power method to find 
       the largest eigenvalue.  
          Caveat is that it performs best when the spectral gap is large 
          (diff between largest and 2nd largest eigenvalues).
       If A is symmetric, can diagonalize A.
       If A is not symmetric, can either diagonalize A^T*A and use the square root 
          of those eigenvalues or can use Singular Value Decomposition on A 
          (the singular values are the eigenvalues of A).
   Note that the eigenvectors are the same for the diagonalization of A, the 
   diagonalization of A^T*A, the SVD(A), and the same operations performed 
   on the CUR-Decompositions of A (=C) or on C^T*C.
   NOTE: for large-scale eigenvalue problems,consider cur decomposition, power
   * method, qr decomposition (which for general matrices might include the 
   * schur decomposition followed by back substitution), or
   * for hermetician matrices can use divide and conquer eigenvalue algorithms and they
   * are parallelizable.

 Some notes on the smallest eigenvalue of A:
 from Wirtz and Guhr 2014,
 “Distribution of the Smallest Eigenvalue in Complex and Real Correlated Wishart Ensembles”
 https://arxiv.org/pdf/1310.2467.pdf
 In linear discriminant analysis it [smallest eigenvalue] gives the leading contribution for the threshold estimate [21].
 It is most sensitive to noise in the data [18].
 In linear principal component analysis, the smallest eigenvalue determines the plane of closest fit [18].
 It is also crucial for the identification of single statistical outliers [17].
 In numerical studies involving large random matrices, the condition number is used, which depends
 on the smallest eigenvalue [22, 23]. In wireless communication the Multi–Input–Multi– Output (MIMO) channel matrix
 of an antenna system is modeled by a random matrix [24]. The smallest eigenvalue of C yields an estimate for the
 error of a received signal [25, 26, 27]. In finance, the optimal portfolio is associated with the eigenvector
 to the smallest eigenvalue of the covariance matrix, which is directly related to the correlation matrix [28].
 This incomplete list of examples shows the influence of the smallest eigenvalue in applications.
 Further information on the role of the smallest eigenvalue is given in Appendix A.
 [17] Barnett V, Lewis T. Outliers in Statistical Data. first edition ed. John Wiley & Sons; 1980.
 [18] Gnanadesikan R. Methods for Statistical Data Analysis of Multivariate Oberservations.
 Second edition ed. John Wiley & Sons; 1997.

 Note also, regarding finding the closest distance of a vector to a subspace in matrix A:
 One can use orthgoncal projection to solve for the distance, which is also called the residual
 or error vector.
     Given matrix A and vector b, we solve for x (i.e. x_est).
             A * x = b
        x_est is the subspace in A closest to b.   It's the best estimate for x.
        x_est = pseudoinverse(A) * b; [A=mxn, b=mx1, x_est=nx1]
        e.g. for full rank
             x_est = (A^T*A)^-1 * A^T * b

       The projection of b onto subspace x_est:
           p = A * x_est.  A=mxn, x_est=nx1, P=mx1

        The projection matrix in P = p*b is then
            P = A * pseudoinverse(A)
            e.g. for full rank
                P = A * (A^T*A)^-1 * A^T; A=mxn, P=mxm

        The vector that is perpendicular to the subspace of x_est is
            b - A * x_est
            this is called the residual, i.e. error vector

 Note that to estimate b - A * x_est, one could instead use eigenvectors (or singular value vectors):
    If an exact solution to b is possible (i.e. A is invertible, etc.), the eigenvectors are solved for matrix A,
    else the eigenvectors are solved for matrix A^T*A (which, if A is zero mean centered, is the covariance of A.
    see PCA for related, but different context).
    We want a vector of length n where A is mxn, so that means we use the right matrix of decomposed eigen/singular
    vectors (e.g. the V matrix in the SVD U, S, V matrices).
    The eigenvector associated with the smallest eigenvalue solves for the null space
    (eigenvalue is 0), and so it is perpendicular to the row space.
    This method is often used in optimization as total least squares regression.
    The x_est that minimizes ||b - A * x_est||^2 is the least squares solution for x_est when m > n.
 The projection matrix of A is A * pseudoinverse(A).


 Note that the eigenvector for the 2nd smallest eigenvalue of A, called the Fielder has distinct uses also,
 especially in community finding.
 Search for Fiedler vector in this code base.
 </pre>
 * 
 * TODO: implement Orthogonal Iteration from Morita and Kanade Sect 3.2.2
   "A straightforward generalization of the power method can be used to compute 
   several dominant eigenvectors of a symmetric matrix."
   A useful alternative to the singular value decomposition in situations 
   where B is a large matrix and only a few of its largest eigenvalues are needed.
   * for p eigenvalues out of n, if lambda__(p_1)/lambda_p is very small,
   * the computation of the eigenvectors should proceed quickly.

 see https://netlib.org/lapack/explore-html/index.html for details of
 MTJ toolkit methods and LAPACK.getInstance() methods that might not be in MTJ

 Note that any rectangular matrix can be made into a square matrix by adding zero rows or columns,
 without changing the nonzero singular values.
 Bjorck 1991, "Algorithms for Linear Least Squares Problems"

 * @author nichole
 */
public class MatrixUtil {
    
    /**
     * multiply the row vector v by matrix m.
     @param v vector
     @param m matrix
     @return result is size [1][m[0].length]
     */
    public static double[] multiplyRowVectorByMatrix(double[] v, double[][] m) {
        
        if (m.length != v.length) {
            throw new IllegalArgumentException("length of v and length of m must be the same");
        }
        
        int mRows = m.length;
        int mCols = m[0].length;
        
        double[] out = new double[mCols];
        for (int i = 0; i < mCols; ++i) {
            for (int j = 0; j < mRows; ++j) {
                out[i] += (v[j] * m[j][i]);
            }
        }
        return out;
    }
    
    /**
     * multiply matrix m by vector n
     @param m two dimensional array in row major format
     @param n one dimensional array
     @return vector of length m.length
     */
    public static double[] multiplyMatrixByColumnVector(double[][] m, double[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mcols = m[0].length;

        int mrows = m.length;

        int ncols = n.length;
        
        if (mcols != ncols) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the length of n");
        }
        
        double[] c = new double[mrows];

        multiplyMatrixByColumnVector(m, n, c);

        return c;
    }
    
    /**
     * multiply matrix m by vector n
     @param m two dimensional array in row major format
     @param n one dimensional array
     @return vector of length m.length
     */
    public static double[] multiplyMatrixByColumnVector(double[][] m, int[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mcols = m[0].length;

        int mrows = m.length;

        int ncols = n.length;
        
        if (mcols != ncols) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the length of n");
        }
        
        double[] c = new double[mrows];

        multiplyMatrixByColumnVector(m, n, c);

        return c;
    }
    
    /**
     * multiply matrix m by vector n and return results in given vector out
     @param m two dimensional array in row major format
     @param n one dimensional array
     @param out vector of length m.length to return results in
     */
    public static void multiplyMatrixByColumnVector(double[][] m, double[] n,
        double[] out) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        // identity check
        if (n.hashCode() == out.hashCode()) {
            System.err.println("warning: n cannot be the same as out. this warning can be a mistake when the " +
                    "hashcodes for truly different objects collide.");
        }
        
        int mcols = m[0].length;

        int mrows = m.length;

        int ncols = n.length;
        
        if (mcols != ncols) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the length of n");
        }
        if (mrows != out.length) {
            throw new IllegalArgumentException(
                "out.length must equal m.length");
        }
        
        Arrays.fill(out, 0);
        
        int cCol = 0;

        for (int row = 0; row < mrows; row++) {                        
            for (int col = 0; col < mcols; col++) {
                out[cCol] += (m[row][col] * n[col]);
            }
            cCol++;        
        }
    }
    
    /**
     * multiply matrix m by vector n and return results in given vector out
     @param m two dimensional array in row major format
     @param n one dimensional array
     @param out vector of length m.length to return results in
     */
    public static void multiplyMatrixByColumnVector(double[][] m, int[] n,
        double[] out) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        // identity check
        if (m.hashCode() == out.hashCode()) {
            System.err.println("warning: m cannot be the same as out. this warning can be a mistake when the " +
                    "hashcodes of truly different objects collide.");
        }
        
        int mcols = m[0].length;

        int mrows = m.length;

        int ncols = n.length;
        
        if (mcols != ncols) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the length of n");
        }
        if (mrows != out.length) {
            throw new IllegalArgumentException(
                "out.length must equal m.length");
        }
        
        Arrays.fill(out, 0);
        
        int cCol = 0;
        
        for (int row = 0; row < mrows; row++) {                        
            for (int col = 0; col < mcols; col++) {
                out[cCol] += (m[row][col] * n[col]);
            }
            cCol++;        
        }
    }
    
    /**
     * multiply matrix m by vector n
     @param m two dimensional array in row major format
     @param n one dimensional array
     @return result
     */
    public static float[] multiplyMatrixByColumnVector(float[][] m, float[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mcols = m[0].length;

        int mrows = m.length;

        int ncols = n.length;
        
        if (mcols != ncols) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the number of rows in n");
        }
        
        float[] c = new float[mrows];

        int cCol = 0;
        
        for (int row = 0; row < mrows; row++) {
            for (int col = 0; col < mcols; col++) {
                c[cCol] += (m[row][col] * n[col]);
            }            
            cCol++;        
        }

        return c;
    }
    
    /**
     *
     @param a
     @param b
     @return
     */
    public static double[] multiplyMatrixByColumnVector(Matrix a, double[] b) {
        
        if (a == null || a.numRows() == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (b == null || b.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = a.numRows();

        int mcols = a.numColumns();

        int nrows = b.length;
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of cols in a must equal the length of b");
        }
        
        double[] c = new double[mrows];

        int cCol = 0;
        
        /*
        a0 1 2     0
                   1
                   2
        */
        
        for (int row = 0; row < mrows; row++) {
            for (int col = 0; col < mcols; col++) {
                c[cCol] += (a.get(row, col) * b[col]);
            }
            cCol++;        
        }

        return c;
    }
    
    /**
     * multiply vector m by factor
     @param m one dimensional array that is input and output for result
     @param factor factor to multiply m by
     */
    public static void multiply(int[] m, int factor) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int len = m.length;
                
        for (int i = 0; i < len; i++) {                        
            m[i] *= factor;
        }
    }
    
    /**
     * multiply matrix m by factor
     @param a two dimensional array in that is input and output for result
     @param m factor to multiply m by
     */
    public static void multiply(float[][] a, float m) {

        if (a == null || a.length == 0) {
            throw new IllegalArgumentException("a cannot be null or empty");
        }
        
        int mcols = a.length;

        int mrows = a[0].length;
        
        for (int col = 0; col < mcols; col++) {
            for (int row = 0; row < mrows; row++) {
                a[col][row] *= m;
            }            
        }
    }
    
    /**
     *
     @param m
     @param n
     */
    public static void multiply(int[] m, int[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        if (m.length != n.length) {
            throw new IllegalArgumentException("m must be the same size as n");
        }
        
        int len = m.length;
                
        for (int i = 0; i < len; i++) {                        
            m[i] *= n[i];
        }
    }
    
    /*
    public static double[][] fastMultiply(double[][] m, double[][] n) {
        
        consider using Native BLAS from the included Netlib package
        BLAS.getInstance().dgemm.
        
        Can see example in this repository's related project called
        https://github.com/nking/curvature-scale-space-corners-and-transformations.git
        in the test class tests//algorithms/NetlibTest.java
        
        or in MTJ source code
        https://github.com/fommil/matrix-toolkits-java/blob/master/src/main/java/no/uib/cipr/matrix/DenseMatrix.java#L299
        
    }*/
    
    /**
     * multiply matrix m by matrix n
     @param m two-dimensional array in row major format
     @param n two-dimensional array in row major format
     @return multiplication of m by n.  resulting matrix is size mrows X ncols.
     */
    public static double[][] multiply(double[][] m, double[][] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = m.length;

        int mcols = m[0].length;

        int nrows = n.length;
        
        int ncols = n[0].length;
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m (=" + mcols + ") "
                + " must equal the number of rows in n (=" + nrows + ")");
        }
        
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0 + b*p3 + c*p6    a*p1 + b*p4 + c*p7    a*p2 + b*p5 + c*p8
        d*p0 + d*p3 + e*p6    d*p1 + d*p4 + e*p7    d*p2 + e*p5 + f*p8
        */
        
        // mrows X ncols
        double[][] c = MatrixUtil.zeros(mrows, ncols);
        
        multiply(m, n, c);

        return c;
    }
    
    /**
     * multiply matrix m by matrix n
     @param m two dimensional array in row major format
     @param n two dimensional array in row major format
     @return multiplication of m by n.  resulting matrix is size mrows X ncols.
     */
    public static double[][] multiply(double[][] m, int[][] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = m.length;

        int mcols = m[0].length;

        int nrows = n.length;
        
        int ncols = n[0].length;
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m (=" + mcols + ") "
                + " must equal the number of rows in n (=" + nrows + ")");
        }
        
        // mrows X ncols
        double[][] c = MatrixUtil.zeros(mrows, ncols);
        
        multiply(m, n, c);

        return c;
    }
    
     /**
     * multiply matrix m by matrix n
     @param m two dimensional array in row major format
     @param n two dimensional array in row major format
     @return multiplication of m by n.  resulting matrix is size mrows X ncols.
     */
    public static double[][] multiply(int[][] m, double[][] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = m.length;

        int mcols = m[0].length;

        int nrows = n.length;
        
        int ncols = n[0].length;
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m (=" + mcols + ") "
                + " must equal the number of rows in n (=" + nrows + ")");
        }
        
        // mrows X ncols
        double[][] c = MatrixUtil.zeros(mrows, ncols);
        
        multiply(m, n, c);

        return c;
    }


    /**
     * multiply matrix m by matrix n
     @param m two dimensional array in row major format
     @param n two dimensional array in row major format
     @return multiplication of m by n.  resulting matrix is size mrows X ncols.
     */
    public static String[][] multiply(String[][] m, String[][] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }

        int mrows = m.length;

        int mcols = m[0].length;

        int nrows = n.length;

        int ncols = n[0].length;

        if (mcols != nrows) {
            throw new IllegalArgumentException(
                    "the number of columns in m (=" + mcols + ") "
                            + " must equal the number of rows in n (=" + nrows + ")");
        }

        String[][] out = new String[mrows][ncols];

        for (int mrow = 0; mrow < mrows; mrow++) {
            for (int ncol = 0; ncol < ncols; ncol++) {
                StringBuilder sum = new StringBuilder();
                for (int mcol = 0; mcol < mcols; mcol++) {
                    if (sum.length() != 0) {
                        sum.append(" + ");
                    }
                    sum.append(String.format("(%s * %s)", m[mrow][mcol], n[mcol][ncol]));
                }
                out[mrow][ncol] = sum.toString();
            }
        }

        return out;
    }

    public static String[][] transpose(String[][] m) {
        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }

        int mRows = m.length;
        int mCols = m[0].length;

        String[][] t = new String[mCols][mRows];
        int i, j;
        for (i = 0; i < mRows; i++) {
            for (j = 0; j < mCols; j++) {
                t[j][i] = m[i][j];
            }
        }

        return t;
    }
    
    /**
     * multiply matrix m by matrix n
     @param m two-dimensional array in row major format
     @param n two-dimensional array in row major format
     @param out the results of multiplication of m by n.  the matrix should be size mrows X ncols.
     */
    public static void multiply(double[][] m, double[][] n, double[][] out) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        // identity check:
        if (m.hashCode() == out.hashCode() || n.hashCode() == out.hashCode()) {
            System.err.println("warning: neither m nor n can be the same as out. this warning can be a mistake when the " +
                    "hashcodes for truly different objects collide.");
        }
        
        int mrows = m.length;

        int mcols = m[0].length;

        int nrows = n.length;
        
        int ncols = n[0].length;
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m (=" + mcols + ") "
                + " must equal the number of rows in n (=" + nrows + ")");
        }
        
        if (out.length != mrows || out[0].length != ncols) {
            throw new IllegalArgumentException("out must be [m.length X n[0].length]");
        }
        
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0 + b*p3 + c*p6    a*p1 + b*p4 + c*p7    a*p2 + b*p5 + c*p8
        d*p0 + d*p3 + e*p6    d*p1 + d*p4 + e*p7    d*p2 + e*p5 + f*p8
        */

        for (int mrow = 0; mrow < mrows; mrow++) {
            for (int ncol = 0; ncol < ncols; ncol++) {
                double sum = 0;                
                for (int mcol = 0; mcol < mcols; mcol++) {
                    sum += (m[mrow][mcol] * n[mcol][ncol]);                    
                }
                out[mrow][ncol] = sum;
            }            
        }
    }
    
    /**
     * multiply matrix m by matrix n
     @param m tow dimensional array in row major format
     @param n two dimensional array in row major format
     @param out the results of multiplication of m by n.  the matrix should be size mrows X ncols.
     */
    public static void multiply(double[][] m, int[][] n, double[][] out) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        // identity check:
        if (m.hashCode() == out.hashCode()) {
            System.err.println("warning: m cannot be the same as out. this warning can be a mistake when the " +
                    "hashcodes for truly different objects collide.");
        }
        
        int mrows = m.length;

        int mcols = m[0].length;

        int nrows = n.length;
        
        int ncols = n[0].length;
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m (=" + mcols + ") "
                + " must equal the number of rows in n (=" + nrows + ")");
        }
        
        if (out.length != mrows || out[0].length != ncols) {
            throw new IllegalArgumentException("out must be [m.length X n[0].length]");
        }
        
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0 + b*p3 + c*p6    a*p1 + b*p4 + c*p7    a*p2 + b*p5 + c*p8
        d*p0 + d*p3 + e*p6    d*p1 + d*p4 + e*p7    d*p2 + e*p5 + f*p8
        */
        
        for (int mrow = 0; mrow < mrows; mrow++) {
            for (int ncol = 0; ncol < ncols; ncol++) {
                double sum = 0;                
                for (int mcol = 0; mcol < mcols; mcol++) {
                    sum += (m[mrow][mcol] * n[mcol][ncol]);                    
                }
                out[mrow][ncol] = sum;
            }            
        }
    }
    
    /**
     * multiply matrix m by matrix n
     @param m tow dimensional array in row major format
     @param n two dimensional array in row major format
     @param out the results of multiplication of m by n.  the matrix should be size mrows X ncols.
     */
    public static void multiply(int[][] m, double[][] n, double[][] out) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        // identity check:
        if (n.hashCode() == out.hashCode()) {
            System.err.println("warning: n cannot be the same as out. this warning can be a mistake when the " +
                    "hashcodes for truly different objects collide.");
        }
        
        int mrows = m.length;

        int mcols = m[0].length;

        int nrows = n.length;
        
        int ncols = n[0].length;
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m (=" + mcols + ") "
                + " must equal the number of rows in n (=" + nrows + ")");
        }
        
        if (out.length != mrows || out[0].length != ncols) {
            throw new IllegalArgumentException("out must be [m.length X n[0].length]");
        }
        
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0 + b*p3 + c*p6    a*p1 + b*p4 + c*p7    a*p2 + b*p5 + c*p8
        d*p0 + d*p3 + e*p6    d*p1 + d*p4 + e*p7    d*p2 + e*p5 + f*p8
        */
        
        for (int mrow = 0; mrow < mrows; mrow++) {
            for (int ncol = 0; ncol < ncols; ncol++) {
                double sum = 0;                
                for (int mcol = 0; mcol < mcols; mcol++) {
                    sum += (m[mrow][mcol] * n[mcol][ncol]);                    
                }
                out[mrow][ncol] = sum;
            }            
        }
    }
    
    /**
     * multiply matrix m by matrix n
     @param m tow dimensional array in row major format
     @param n two dimensional array in row major format
     @return out the results of multiplication of m by n.  the matrix should be size mrows X ncols.
     */
    public static int[][] multiply(int[][] m, int[][] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = m.length;

        int mcols = m[0].length;

        int nrows = n.length;
        
        int ncols = n[0].length;
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m (=" + mcols + ") "
                + " must equal the number of rows in n (=" + nrows + ")");
        }
        
        int[][] out = new int[mrows][];
        int mrow, ncol;
        for (mrow = 0; mrow < mrows; mrow++) {
            out[mrow] = new int[ncols];
        }
        
        for (mrow = 0; mrow < mrows; mrow++) {
            for (ncol = 0; ncol < ncols; ncol++) {
                int sum = 0;                
                for (int mcol = 0; mcol < mcols; mcol++) {
                    sum += (m[mrow][mcol] * n[mcol][ncol]);                    
                }
                out[mrow][ncol] = sum;
            }            
        }
        return out;
    }
    
    /**
     *
     @param a
     @return
     */
    public static boolean isOrthogonal(int[][] a) {
        
        int[][] aT = MatrixUtil.transpose(a);
        
        int[][] prod;
        int i, n;
        for (int k = 0; k < 2; ++k) {
            if (k == 0) {
                prod = MatrixUtil.multiply(a, aT);
            } else {
                prod = MatrixUtil.multiply(aT, a);
            }
            n = prod.length;
            assert(n == prod[0].length);
            for (i = 0; i < n; ++i) {
                if (prod[i][i] != 1) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    /**
     *
     @param a
     @return
     */
    public static boolean isAPermutationMatrix(int[][] a) {
        
        int i, j;
        int sum = 0;
        for (i = 0; i < a.length; ++i) {
            sum = 0;
            for (j = 0; j < a[i].length; ++j) {
                sum += a[i][j];
            }
            if (sum != 1) {
                return false;
            }
        }
        for (j = 0; j < a[0].length; ++j) {
            sum = 0;
            for (i = 0; i < a.length; ++i) {
                sum += a[i][j];
            }
            if (sum != 1) {
                return false;
            }
        }

        return true;
    }

    /**
     *
     @param a
     @return
     */
    public static double[][] createATransposedTimesA(double[][] a) {
        if (a == null || a.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }

        int m = a.length;
        int n = a[0].length;

        double[][] out = MatrixUtil.zeros(n, n);
        createATransposedTimesA(a, out);
        return out;
    }

    /**
     *
     @param a
     @param out
     */
    public static void createATransposedTimesA(double[][] a, double[][] out) {

        if (a == null || a.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int m = a.length;
        int n = a[0].length;
        if (out.length != n || out[0].length != n) {
            throw new IllegalArgumentException("out must be size a[0].length X a[0].length");
        }

        /*
        a00  a01     a00  a01 
        a10  a11     a10  a11 
        a20  a21     a20  a21 
        
                         i:[0,n)
        a00  a10  a20    a00  a01
        a01  a11  a21    a10  a11
                         a20  a21 
        
        a00*a00 + a10*a10 + a20*a20   a00*a01 + a10*a11 + a20*a21   
        a01*a00 + a11*a10 + a21*a20   a01*a01 + a11*a11 + a21*a21 
        
        col0_ dot col0
        col1_ dot col0
        */
        
        int outCol, i, j;
        double sum;
        for (i = 0; i < n; i++) {
            for (outCol = 0; outCol < n; outCol++) {
                sum = 0;
                for (j = 0; j < m; j++) {
                    sum += (a[j][outCol] * a[j][i]);
                }
                out[outCol][i] = sum;
            }
        }
    }
    
    /**
     *
     @param a
     @param f
     */
    public static void multiply(double[] a, double f) {
        for (int i = 0; i < a.length; ++i) {
            a[i] *= f;
        }
    }
    
    /**
     *
     @param m
     @param factor
     */
    public static void multiply(double[][] m, double factor) {

        int nrows = m.length;
        int ncols = m[0].length;

        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                m[i][j] = factor*m[i][j];
            }
        }
    }
    
    /**
     * multiply matrices and return matrix of size mrows X ncols
     @param m matrix
     @param n matrix
     @return result
     */
    public static DenseMatrix multiply(Matrix m, Matrix n) {

        if (m == null || m.numRows() == 0 || m.numColumns() == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.numRows() == 0 || n.numColumns() == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = m.numRows();

        int mcols = m.numColumns();

        int nrows = n.numRows();
        
        int ncols = n.numColumns();
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the number of rows in n");
        }
    
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0+... a*p a*p
        d*p0+... d*p d*p
        */
        
        no.uib.cipr.matrix.DenseMatrix c = new DenseMatrix(mrows, ncols);        
        for (int row = 0; row < mrows; row++) {
            for (int ncol = 0; ncol < ncols; ncol++) {
                double sum = 0;                
                for (int mcol = 0; mcol < mcols; mcol++) {
                    sum += (m.get(row, mcol) * n.get(mcol, ncol));                    
                }
                c.set(row, ncol, sum);
            }            
        }

        return c;
    }
    
    /**
     *
     @param a
     @param f
     */
    public static void multiply(TDoubleArrayList a, double f) {
        for (int i = 0; i < a.size(); ++i) {
            a.set(i, f * a.get(i));
        }
    }
    
    /**
     *
     @param a
     @param b
     */
    public static void multiply(Matrix a, double b) {
        
        if (a == null || a.numRows() == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        Iterator<MatrixEntry> iter = a.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            entry.set(entry.get() * b);
        }
        
    }
    /**
     * perform dot product of m and a diagonalized matrix of diag,
     * and return matrix of size mrows X mcols
     @param m matrix
     @param diag vector
     @return result
     */
    public static double[][] multiplyByDiagonal(double[][] m, double[] diag) {
        double[][] out = MatrixUtil.zeros(m.length, m[0].length);
        multiplyByDiagonal(m, diag, out);
        return out;
    }

    /**
     * multiply the diagonal matrix of diag by the matrix m and return the results in the given output matrix.
     @param diag an array holding the diagonal elements of the diagonal matrix.  size is length a.length.
     @param a an mxn matrix.
     @return
     */
    public static double[][] multiplyDiagonalByMatrix(double[] diag, double[][] a) {
        double[][] out = MatrixUtil.zeros(a.length, a[0].length);
        multiplyDiagonalByMatrix(diag, a, out);
        return out;
    }

    /**
     * multiply the diagonal matrix of diag by the matrix m and return the results in the given output matrix.
     @param diag an array holding the diagonal elements of the diagonal matrix.  size is length a.length.
     @param a an mxn matrix.
     @param out output matrix of size a.length X a[0].length
     */
    public static void multiplyDiagonalByMatrix(double[] diag, double[][] a, double[][] out) {
        int m = a.length;
        int n = a[0].length;
        if (diag.length != m) {
            throw new IllegalArgumentException("diag must be length a.length");
        }
        if (out.length != m) {
            throw new IllegalArgumentException("out.length must equal a.length");
        }
        if (out[0].length != n) {
            throw new IllegalArgumentException("out[0].length must equal a[0].length");
        }
        int i;
        int j;
        for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j) {
                out[i][j] = diag[i] * a[i][j];
            }
        }
    }
    
    /**
     * perform dot product of m and a diagonalized matrix of diag,
     * and return matrix of size mrows X mcols
     @param m matrix
     @param diag vector
     @param out
     */
    public static void multiplyByDiagonal(double[][] m, double[] diag, double[][] out) {

        if (m == null || m.length == 0 || m[0].length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (diag == null || diag.length == 0) {
            throw new IllegalArgumentException("diag cannot be null or empty");
        }
        
        int mrows = m.length;

        int mcols = m[0].length;

        int nrows = diag.length;
                
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the number of rows in n");
        }

        if (out.length != mrows || out[0].length != mcols) {
            throw new IllegalArgumentException(
                    "out size must be [" + mrows+ " x " + mcols + "]");
        }
        
        /*
        a b c      p0 0  0
        d e f      0  p1 0
                   0  0  p2        
        */

        for (int row = 0; row < mrows; row++) {
            out[row] = new double[mcols];
            for (int mcol = 0; mcol < mcols; mcol++) {
                out[row][mcol] = m[row][mcol] * diag[mcol];
            }            
        }
    }
    
    /**
     * perform dot product of m and a diagonalized matrix of diag,
     * and return matrix of size mrows X mcols
     @param m matrix
     @param diag diagonal array to use in diagonal matrix
     */
    public static void multiplyByDiagonal(DenseMatrix m, double[] diag) {

        if (m == null || m.numRows() == 0 || m.numColumns() == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (diag == null || diag.length == 0) {
            throw new IllegalArgumentException("diag cannot be null or empty");
        }
        
        int mrows = m.numRows();

        int mcols = m.numColumns();

        int nrows = diag.length;
                
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the number of rows in n");
        }
        
        /*
        a b c      p0 0  0
        d e f      0  p1 0
                   0  0  p2        
        */
                
        for (int row = 0; row < mrows; row++) {
            for (int mcol = 0; mcol < mcols; mcol++) {
                m.set(row, mcol,  m.get(row, mcol) * diag[mcol]);
            }            
        }
    }
    
    /**
     * calculates the inner product of a and b, which is a as a single row matrix
     * and b as a single column matrix, so is a^T * b.  it's also known as the
     * scalar product or dot product.
     @param a list
     @param b list
     @return scale result of a^T * b
     */
    public static double innerProduct(TDoubleArrayList a, 
        TDoubleArrayList b) {
        int sz0 = a.size();
        int sz1 = b.size();
        if (sz0 != sz1) {
            throw new IllegalArgumentException(
                "a and b must be same size");
        }
        double s = 0;
        for (int i = 0; i < a.size(); ++i) {
            s += a.get(i) * b.get(i);
        }
        return s;
    }
    
    /**
     * calculates the inner product of a and b, which is a as a single row matrix
     * and b as a single column matrix, so is a^T * b.
     @param a array
     @param b array
     @return scalar result of a^T * b
     */
    public static double innerProduct(double[] a, 
        double[] b) {
        int sz0 = a.length;
        int sz1 = b.length;
        if (sz0 != sz1) {
            throw new IllegalArgumentException(
                "a and b must be same size");
        }
        double s = 0;
        for (int i = 0; i < a.length; ++i) {
            s += a[i] * b[i];
        }
        return s;
    }
    
    /**
     *
     @param a
     @param b
     @return
     */
    public static double innerProduct(int[] a, double[] b) {

        if (a.length != b.length) {
            throw new IllegalArgumentException("a.length must == b.length");
        }
        
        double sum = 0;
        for (int i = 0; i < a.length; ++i) {
            sum += (a[i] * b[i]);
        }

        return sum;
    }
     
    /**
     *
     @param a
     @param b
     @return
     */
    public static TDoubleArrayList subtract(TDoubleArrayList a, TDoubleArrayList b) {
        int sz0 = a.size();
        int sz1 = b.size();
        if (sz0 != sz1) {
            throw new IllegalArgumentException(
                "a and b must be same size");
        }
        TDoubleArrayList c = new TDoubleArrayList(sz0);
        for (int i = 0; i < a.size(); ++i) {
            c.add(a.get(i) - b.get(i));
        }
        return c;
    }
    
    /**
     *
     @param m
     @param n
     @return
     */
    public static double[] subtract(double[] m, double[] n) {

        int len = m.length;

        double[] c = new double[len];

        subtract(m, n, c);

        return c;
    }
    
    /**
     * calculate m[i] - s for each i=[0, m.length).
     @param m array
     @param s scalar
     @return result
     */
    public static double[] subtract(double[] m, double s) {

        int len = m.length;

        double[] c = new double[len];

        for (int i = 0; i < len; ++i) {
            c[i] = m[i] - s;
        }

        return c;
    }
    /**
     * calculate s - m[i] for each i=[0, m.length).
     @param m array
     @param s scalar
     @return result
     */
    public static double[] subtract(double s, double[] m) {

        int len = m.length;

        double[] c = new double[len];

        for (int i = 0; i < len; ++i) {
            c[i] = s - m[i];
        }

        return c;
    }
    
    /**
     * calculate s + m[i] for each i=[0, m.length).
     @param m array
     @param s scalar to add to m
     @return resulting array
     */
    public static double[] add(double s, double[] m) {

        int len = m.length;

        double[] c = new double[len];

        for (int i = 0; i < len; ++i) {
            c[i] = s + m[i];
        }

        return c;
    }
    
    /**
     *
     @param m
     @param n
     @param output
     */
    public static void subtract(double[] m, double[] n,
        double[] output) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.length != n.length) {
            throw new IllegalArgumentException("m and n must be same length");
        }

        int len = m.length;

        for (int i = 0; i < len; i++) {
            output[i] = m[i] - n[i];
        }
    }
    
    /**
     *
     @param m
     @param n
     @return
     */
    public static double[] add(double[] m, double[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.length != n.length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        int len = m.length;
     
        double[] c = new double[len];
        
        for (int i = 0; i < len; i++) {
            c[i] = m[i] + n[i];
        }

        return c;
    }
    
    /**
     *
     @param m
     @param n
     @return
     */
    public static float[] add(float[] m, float[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.length != n.length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        int len = m.length;
     
        float[] c = new float[len];
        
        for (int i = 0; i < len; i++) {
            c[i] = m[i] + n[i];
        }

        return c;
    }
    
    /**
     *
     @param m
     @param n
     @return
     */
    public static float[][] subtract(float[][] m, float[][] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.length != n.length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        if (m[0].length != n[0].length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        float[][] c = new float[m.length][];

        for (int i = 0; i < m.length; ++i) {
            c[i] = new float[m[0].length];
            for (int j = 0; j < m[0].length; ++j) {
                c[i][j] -= m[i][j] - n[i][j];
            }
        }

        return c;
    }
    
    /**
     *
     @param m
     @param n
     @return
     */
    public static float[][] add(float[][] m, float[][] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.length != n.length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        if (m[0].length != n[0].length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        float[][] c = new float[m.length][];

        for (int i = 0; i < m.length; ++i) {
            c[i] = new float[m[0].length];
            for (int j = 0; j < m[0].length; ++j) {
                c[i][j] = m[i][j] + n[i][j];
            }
        }

        return c;
    }

    /**
     *
     @param m
     @param n
     @return
     */
    public static float[] subtract(float[] m, float[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.length != n.length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        int len = m.length;
     
        float[] c = new float[len];
        
        for (int i = 0; i < len; i++) {
            c[i] = m[i] - n[i];
        }

        return c;
    }
    
    /**
     *
     @param m
     @param n
     @return
     */
    public static DenseMatrix subtract(DenseMatrix m, DenseMatrix n) {

        if (m == null) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.numRows() != n.numRows() || m.numColumns() != n.numColumns()) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        DenseMatrix output = new DenseMatrix(m.numRows(), m.numColumns());
        
        for (int i = 0; i < m.numRows(); ++i) {
            for (int j = 0; j < m.numColumns(); ++j) {
                double v0 = m.get(i, j);
                double v1 = n.get(i, j);
                output.set(i, j, v0 - v1);
            }
        }
        
        return output;
    }
    
    /**
     *
     @param m
     @param n
     */
    public static void add(int[] m, int n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int len = m.length;
             
        for (int i = 0; i < len; i++) {
            m[i] += n;
        }
    }
    
    /**
     *
     @param m
     @return
     */
    public static float[][] transpose(float[][] m) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }

        int mRows = m.length;
        int mCols = m[0].length;

        float[][] t = new float[mCols][];
        for (int i = 0; i < mCols; i++) {
            t[i] = new float[mRows];
        }

        for (int i = 0; i < mRows; i++) {
            for (int j = 0; j < mCols; j++) {
                t[j][i] = m[i][j];
            }
        }

        return t;
    }
    
    /**
     *
     @param m
     @return
     */
    public static DenseMatrix transpose(DenseMatrix m) {
        
        int mRows = m.numRows();
        int mCols = m.numColumns();
        
        double[][] t = new double[mRows][];
        for (int i = 0; i < mRows; i++) {
            t[i] = new double[mCols];
            for (int j = 0; j < mCols; j++) {
                t[i][j] = m.get(i, j);
            }
        } 
        
        double[][] transposed = transpose(t);
        
        DenseMatrix mT = new DenseMatrix(transposed);
        
        return mT;
    }
    
    /**
     *
     @param m
     @return
     */
    public static double[][] transpose(double[][] m) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int mRows = m.length;
        int mCols = m[0].length;
        
        double[][] t = new double[mCols][];
        for (int i = 0; i < mCols; i++) {
            t[i] = new double[mRows];
        }
        
        for (int i = 0; i < mRows; i++) {
            for (int j = 0; j < mCols; j++) {
                t[j][i] = m[i][j];
            }
        }
        
        return t;
    }
    
    /**
     *
     @param m
     @return
     */
    public static int[][] transpose(int[][] m) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int mRows = m.length;
        int mCols = m[0].length;
        
        int[][] t = new int[mCols][];
        for (int i = 0; i < mCols; i++) {
            t[i] = new int[mRows];
        }
        int i, j;
        for (i = 0; i < mRows; i++) {
            for (j = 0; j < mCols; j++) {
                t[j][i] = m[i][j];
            }
        }
        
        return t;
    }
    
    /**
     * transpose matrix m into out matrix which must be size m[0].length X m.length.
     @param m matrix to transpose
     @param out output matrix to hold transposed m
     */
    public static void transpose(double[][] m, double[][] out) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int mRows = m.length;
        int mCols = m[0].length;
        
        if (out.length != mCols || out[0].length != mRows){
            throw new IllegalArgumentException("out size must be m[0].length X m.length");
        }
        
        for (int i = 0; i < mRows; i++) {
            for (int j = 0; j < mCols; j++) {
                out[j][i] = m[i][j];
            }
        }        
    }
    
    /**
     *
     @param a
     @return
     */
    public static double[][] convertToRowMajor(Matrix a) {
        int nc = a.numColumns();
        int nr = a.numRows();
        double[][] out = new double[nr][];
        for (int i = 0; i < nr; ++i) {
            out[i] = new double[nc];
            for (int j = 0; j < nc; ++j) {
                out[i][j] = a.get(i, j);
            }
        }
        return out;
    }

    /*
    columns of mxn matrix A are linearly independent only when rank r == n.
    there are n pivots and no free variables.  \in this case m is .geq. r).
    only x=0 is in the null space when n == r.
    
    the columns of A are independent if x=0 is the only solution to a*x=0.
    
    <pre>
    
      A is m x n matrix with rank r
    
                            r           n
        .                   .           .
        |                               |  row space is size        r   x  n
        |                               |  col space is size        m   x  r
    r ..|                               |  null space is size     (n-r) x  n
        |                               |  left null space is size  m   x  (m-r)
    m ..|                               |
    
    1) space r x n transformed    : C(A^T) is row space n x r           ==> all A^T*y
    2) space (n-r) x r            : N(A) is null space n x (n-r)        ==> A*x = 0
    3) space m x r                : C(A) is column space m x r          ==> all A*x
    4) space m x (m-r) transformed: N(A^T) is left null space m x (m-r) ==> A^T*y = 0
    
    </pre>
   .
    */

    /*
    (add MIT reference here for 
    Inverse and PseudoInverses

    A is an [mxn] matrix

    2-sided inverse:
       A * A^-1 = I = A^-1 * A
       rank r = m = n
       matrix has full rank

    Left inverse:
       A is full column rank matrix, that is, 
       rank r = n
       the null space of A contains just the zero vector
       A * x = b has exactly 1 solution or is not solvable
       the columns of A are independent
       A^T*A is [nxn] and invertible
          (A^T*A)^-1 * (A^T*A) = I
       left inverse = (A^T*A)^-1 * A^T
           note that there may be other left inverses

    Right inverse:
       A is full row rank matrix, that is
       rank r = m
       the null space of A^T contains just the zero vector
       the rows of A are independent
       A * x = b has at least 1 solution
       the null space has n - m free variables and if n>m has 
          infinitely many solutions
       right inverse = A^T * (A * A^T)^-1

    Pseudo inverse:
       = V * Sigma^+  * U^T
       where SVD(A) = U, Sigma, V^T  (that is, A = U * Sigma * V^T)
       and Sigma^* = Sigma with each element on the diagonal inverted,
           that is Sigma^* = diagonal matrix with [1/Sigma[0][0], 1./Sigma[1][1], ...]
    */
    
    /**
     * calculate the pseudo-inverse of matrix a for cases when the rank of a
     * is less than the width of matrix a.
     * Note that the term rank can be deceptive for cases when the 
     * original matrix a was rank deficient and then is perturbed to become
     * a full rank matrix leading to possibility of larger errors when treated as full rank.
     * the term rank deficient can be replaced by the ter numerically rank deficient.
     * (see Bjork 1991 Section 2, "Algorithms for linear least squares problems"
     * and Chap 6 of Golub and Van Loan).
     * This method uses the SVD of a,
     * specifically, V*R*U^T where R is 1/diagonal of S for cases where
     * rank .leq. m or rank .leq. n where mXn are the dimensions of matrix a.
     * Note that if A^-1 exists, then the pseudo-inverse of A is equal to the
     * inverse of A.
     * 
     * Following Gilbert Strang's "Introduction to Linear Algebra".
     * 
     * TODO: read "ALTERNATIVE METHODS OF CALCULATION OF THE PSEUDO INVERSE
       OF A NON FULL-RANK MATRIX" by M. A. Murray-Lasso, 2008
       http://www.scielo.org.mx/pdf/jart/v6n3/v6n3a4.pdf
       
      NOTE: if the rank is found to be equal to a[0].length and 
      a.length .geq. a[0].length, the full rank pseudo-inverse
      * is calculated instead of the rank-deficient;
     @param a and m X n matrix
     @return pseudo-inverse of matrix a. dimensions are those of a^T.
     * @throws NotConvergedException thrown by MTJ when SVD could not converge
     */
    public static double[][] pseudoinverseRankDeficient(double[][] a) throws NotConvergedException {
        return pseudoinverseRankDeficient(a, true);
    }
    
    /**
     * calculate the pseudo-inverse of matrix a for cases when the rank of a
     * is less than the width of matrix a.
     * Note that the term rank can be deceptive for cases when the 
     * original matrix a was rank deficient and then is perturbed to become
     * a full rank matrix leading to possibility of larger errors when treated as full rank.
     * the term rank deficient can be replaced by the term numerically rank deficient.
     * (see Bjork 1991 Section 2, "Algorithms for linear least squares problems"
     * and Chap 6 of Golub and Van Loan).
     * This method uses the SVD of a,
     * specifically, V*R*U^T where R is 1/diagonal of S for cases where
     * rank .leq. m or rank .leq. n where mXn are the dimensions of matrix a.
     * Note that if A^-1 exists, then the pseudo-inverse of A is equal to the
     * inverse of A.
     * 
     * Following Gilbert Strang's "Introduction to Linear Algebra".
     * 
     * TODO: read "ALTERNATIVE METHODS OF CALCULATION OF THE PSEUDO INVERSE
       OF A NON FULL-RANK MATRIX" by M. A. Murray-Lasso, 2008
       http://www.scielo.org.mx/pdf/jart/v6n3/v6n3a4.pdf
       
      NOTE: if the rank is found to be equal to a[0].length and 
      a.length .geq. a[0].length, the full rank pseudo-inverse
      * is calculated instead of the rank-deficient;
     @param a and m X n matrix
     @param checkForFullRank if true, the method looks for a[0],length == rank
     * and if it is full-rank, the method returns the results form pseudoinverseFullRank(a)
     * instead.
     @return pseudo-inverse of matrix a. dimensions are those of a^T.
     * @throws NotConvergedException thrown by MTJ when SVD could not converge
     */
    static double[][] pseudoinverseRankDeficient(double[][] a,
        boolean checkForFullRank) throws NotConvergedException {
        int m = a.length;
        int n = a[0].length;
        
        // limit for a number to be significant above 0
        double eps = 1e-16;
        
        // from Gilbert Strang's "Introduction to Linear Algebra", Chap 7:
        // uses SVD:
        //
        //   A_inverse = V * pseudoinverse(S) * U^T
        //        where pseudoinverse(S) is simply an empty matrix with the diagonal
        //        being the reciprocal of each singular value
        //      NOTE: compare number of ops w/ this factoring:
        //          V * (pseudoinverse(S) * U^T)
        //   A_inverse is n X m
        //
        //   V is n X n
        //   pseudoinverse(S) is n X m
        //   U^T is m X m
        
        DenseMatrix aMatrix = new DenseMatrix(a);
        SVD svd = SVD.factorize(aMatrix);
        
        //TODO: rewrite to use fewer data structures and multiply in place.
        
        // s is an array of size min(m,n)
        double[] s = svd.getS();
        int rank = 0;
        for (double sv : s) {
            if (sv > eps) {
                rank++;
            }
        }
        if (checkForFullRank && rank == n && m >= n) {
            // use full rank solution:
            return pseudoinverseFullColumnRank(a);
        }
        
        DenseMatrix vTM = svd.getVt();
        double[][] vT = MatrixUtil.convertToRowMajor(vTM);
        double[][] v = MatrixUtil.transpose(vT);
        double[][] uT = MatrixUtil.convertToRowMajor(svd.getU());
        
        /*
        U is mxm orthonormal columns
        S is mxn with non-negative singular values.  rank is number of non-zero entries
        V is  nxn
        */
        assert(v.length == n);
        assert(v[0].length == n);
        assert(uT.length == m);
        assert(uT[0].length == m);
        
        double[][] sInverse = zeros(n, m);
        double sI;
        for (int i = 0; i < s.length; ++i) {
            sI = s[i];
            if (sI > eps) {
                sInverse[i][i] = 1./sI;
            }
        }
        
        /*
        U is mxn orthonormal columns
        S is nxn with non-negative singular values.  rank is number of non-zero entries
        V is  nxn
        pseudoinverse(S) is nxm
        */
        
        //A_inverse = V * pseudoinverse(S) * U^T
        double[][] inv = MatrixUtil.multiply(v, sInverse);
        inv = MatrixUtil.multiply(inv, uT);
        
        return inv;
    }
    
    /**
     * calculate the pseudo-inverse of matrix a (dimensions mxn) which is a full
     * column rank matrix or overdetermined, 
     * n .geq. rank.
     * A_pseudoinverse for A being full column rank = inverse(A^T*A) * A^T.
     * This particular pseudoinverse constitutes a left inverse.
     * pseudoinv(A)*A = I.
     * 
     * If inverting A^T*A fails, the method returns results of pseudoinverseRankDeficient().
     * 
     * NOTE that (A^T*A) (or (A * A^T)) has to be invertible, that is, 
     * the reduced echelon form of A has linearly independent columns (rank==n).
     * following pseudocode from Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms.
     @param a two dimensional array in row major format with dimensions
     * m x n.  a is a full-rank matrix.
     * a is a non-singular matrix(i.e. has exactly one solution). 
     * 
     @return matrix of size [a[0].length][a.length]
     @throws NotConvergedException thrown by MTJ when SVD could not converge
     */
    public static double[][] pseudoinverseFullColumnRank(double[][] a) throws NotConvergedException {
        int m = a.length;
        int n = a[0].length;
        
        // limit for a number to be significant above 0 (precision of computer)
        double eps = 1e-16;
        
        //from Cormen, Leiserson, Rivest, and Stein: A_pseudoinverse = inverse(A^T*A) * A^T
        double[][] _aT = MatrixUtil.transpose(a);
        double[][] _aTA = MatrixUtil.multiply(_aT, a);
        DenseMatrix aTA = new DenseMatrix(_aTA);
        
        //NOTE that (A^T*A) has to be invertible, that is, the reduced echelon form
        // of A has linearly independent columns (no free variables, only pivots.
        // which also means rank==n).
        //
        // could invert (A^T*A) using the cofactor matrix/determinant
        //   or a convenience method from MTJ Matrix.solve
        
        DenseMatrix I = Matrices.identity(_aTA[0].length);
        DenseMatrix placeholder = I.copy();
             
        // the pseudoinverse needs (A^T*A)^-1:
        //    in general: if an inverse exists for a matrix V then  V * V^-1 = I
        //    if an inverse exists for (A^T*A) then (A^T*A) * (A^T*A)^-1 = I
                
        // matlab and other matrix notation uses left division symbol in this way:
        //     x = A\B solves the system of linear equations A*x = B for x.
        //        can substitute A = A^T*A,   X = (A^T*A)^-1,   B = I
        //  X = A\B in MTJ is X = A.solve(B, X), that is, inputs are A and B.
        try {
             DenseMatrix aTAI = (DenseMatrix)aTA.solve(I, placeholder);
             double[][] _aTAI = MatrixUtil.convertToRowMajor(aTAI);
            // _pseudoinverse = inverse(A^T*A) * A^T
            double[][] inv = MatrixUtil.multiply(_aTAI, _aT);
            return inv;
        } catch (no.uib.cipr.matrix.MatrixSingularException ex) {
            return pseudoinverseRankDeficient(a, false);
        }
    }
    
    /**
     * calculate the pseudo-inverse of matrix a (dimensions mxn) which is a full
     * row rank matrix, m .geq. rank.
     * A_pseudoinverse for A being full row rank = a^T*inverse(A*A^T).
     * This particular pseudoinverse constitutes a right inverse.
     * A*pseudoinv(A) = I.
     * 
     @param a two dimensional array in row major format with dimensions
     * m x n.  a is a full-rank matrix.
     * a is a non-singular matrix(i.e. has exactly one solution). 
     * 
     @return matrix of size [a[0].length][a.length]
     * @throws NotConvergedException thrown by MTJ when SVD could not converge
     */
    public static double[][] pseudoinverseFullRowRank(double[][] a) throws NotConvergedException {
        int m = a.length;
        int n = a[0].length;
        
        // limit for a number to be significant above 0 (precision of computer)
        double eps = 1e-16;
        
        //A^T(AA^T)^-1
        double[][] _aT = MatrixUtil.transpose(a);
        double[][] _aAT = MatrixUtil.multiply(a, _aT);
        DenseMatrix aAT = new DenseMatrix(_aAT);
        
        DenseMatrix I = Matrices.identity(_aAT[0].length);
        DenseMatrix placeholder = I.copy();
             
        // the pseudoinverse needs (A*A^T)^-1:
                
        // matlab and other matrix notation uses left division symbol in this way:
        //     x = A\B solves the system of linear equations A*x = B for x.
        //        can substitute A = A^T*A,   X = (A^T*A)^-1,   B = I
        //  X = A\B in MTJ is X = A.solve(B, X), that is, inputs are A and B.
        try {
             DenseMatrix aATI = (DenseMatrix)aAT.solve(I, placeholder);
             double[][] _aATI = MatrixUtil.convertToRowMajor(aATI);
            // _pseudoinverse = inverse(A^T*A) * A^T
            double[][] inv = MatrixUtil.multiply(_aT, _aATI);
            return inv;
        } catch (no.uib.cipr.matrix.MatrixSingularException ex) {
            return pseudoinverseRankDeficient(a, false);
        }
        
    }
    
    /**
      from Strang "Introduction to Linear Algebra":
      <pre>
       an inverse matrix may or may not exist.  
       (1) has to be a square matrix.
           A^-1 x A = I, where I is the identity matrix.
       (2) an inverse matrix has n pivots remaining after elimination,
                where pivot is the leftmost non-zero variable. i.e. the rank r
                is equal to the dimension of the square matrix which is n.
       (3) after elimination, next test is that the determinant is not zero

       if A is invertible, then A * x = b can be solved as x = A^-1 * b
       and (A * B)^-1 = B^-1 * A^-1
      </pre>
     NOTE: because this uses decomposition, each application using it should decide whether
     * to perform the exterior operations at same time to avoid recomputing
     * any matrices.
     @param a matrix
     @return true if is invertible
     */
    public static boolean isInvertible(double[][] a) {
        
        int m = a.length;
        int n = a[0].length;
        if (m != n) {
            throw new IllegalArgumentException("a is not a squae matrix");
        }
        
        //rank:
        // -- could use Cholesky factorization followed by check for positive definiteness
        // -- could use LUP decomposition and count the L diagonal 1's
        // -- could use SVD and count the unique singular values in D
        LUP lup = LinearEquations.LUPDecomposition(a);
        
        int rank = 0;
        for (int i = 0; i < lup.ell.length; ++i) {
            if (Math.abs(lup.ell[i][i] - 1.0) < 1.e-7) {
                rank++;
            }
        }
        if (rank != n) {
            return false;
        }
        
        // determinant (-1)^n * det(U) where det(U) is product of diagonal
        double det = 1;
        for (int i = 0; i < lup.ell.length; ++i) {
            det *= lup.u[i][i];
        }
        if ((n % 2) != 0) {
            det *= -1;
        }
        
        if (Math.abs(det) > 1e-15) {
            // determinant is not zero, so is invertible
            return true;
        }
        
        return false;
    }
    
    /**
      determine the rank of martix A
      <pre>
      If matrix A is a square matrix:
          uses LUP decomposition and counts the L diagonal 1's
      Else
          uses SVD and counts the non-zero diagonal singular values.
          
      NOTE: a positive tolerance level eps is used to find the number of singular values
      above eps instead of 0.
      </pre>
     @param a matrix
     @param eps a positive number for the tolerance above zero of the pivots or
     * singular values.
     @return the rank of A
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static int rank(double[][] a, double eps) throws NotConvergedException {
       
        int m = a.length;
        int n = a[0].length;

        if (m == n) {
            int rU = 0;
            LU lu = LinearEquations.LUDecomposition(a);
            for (int i = 0; i < lu.u.length; ++i) {
                if (Math.abs(lu.u[i][i] - 1.0) < eps) {
                    rU++;
                }
            }
            int rL = 0;
            for (int i = 0; i < lu.ell.length; ++i) {
                if (Math.abs(lu.ell[i][i] - 1.0) < eps) {
                    rL++;
                }
            }
            return Math.min(rU, rL);
        }
        
        /* or
        QRP qrp = QRP.factorize(new DenseMatrix(a));
        final double EPS = 1e-12;
        int rank2 =  qrp.getRank();
        */
        int rank = rank(a);
        return rank;
    }

    /**
     * given a map called adj having keys and values for each key, 
     * create a map where the keys are adj.values and the
     * values are the keys of adj.values.
     @param adj a map
     @return the reverse map, this is the relationships of key to value are reersed to value to key
     */
    public static TIntObjectMap<TIntSet> createReverseMap(TIntObjectMap<TIntSet> adj) {
        TIntObjectMap<TIntSet> r = new TIntObjectHashMap<TIntSet>();
        
        TIntObjectIterator<TIntSet> iter = adj.iterator();
        TIntIterator iterV;
        TIntSet vSet, rVSet;
        int u, v;
        while (iter.hasNext()) {
            iter.advance();
            
            u = iter.key();
            vSet = iter.value();
            
            if (vSet == null || vSet.isEmpty()) {
                continue;
            }
            
            iterV = vSet.iterator();
            while (iterV.hasNext()) {
                v = iterV.next();
                
                rVSet = r.get(v);
                
                if (rVSet == null) {
                    rVSet = new TIntHashSet();
                    r.put(v, rVSet);
                }
                rVSet.add(u);
            }
        }
        
        return r;
    }

    /**
     * class to hold the results of the Singular Value Decomposition
     */
    public static class SVDProducts {
        
        /**
         * the rank of the input matrix
         */
        public int rank;
        
        /**
         * given A as an mxn matrix, u is orthogonal column vectors in a matrix of size mxm.
         * These are eigenvectors of A as columns ordered by the eigenvalues s.
         */
        public double[][] u;
        
        /**
         * given A as an mxn matrix, v^T is orthogonal row vectors in a matrix of sizr nxn.
         * These are eigenvectors of A as rows ordered by the eigenvalues s.
         */
        public double[][] vT;
        
        /**
         * given a matrix A of size mxn, s holds the singular values, that is,
         * the eigenvalues of the SVD decomposition ordered from largest to smallest.
         * these are the diagonal of matrix sigma which may or may not be populated
         */
        public double[] s;
        
        /**
         * the diagonal matrix holding the eigenvalues.  this variable might not
         * be populated if s is.
         */
        public double[][] sigma = null;
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("u=\n");
            if (u != null) {
                sb.append(FormatArray.toString(u, "%.4e"));
            }
            sb.append("\ns=\n");
            if (s != null) {
                sb.append(FormatArray.toString(s, "%.4e"));
            }
            sb.append("\nsigma=\n");
            if (sigma != null) {
                sb.append(FormatArray.toString(sigma, "%.4e"));
            }
            sb.append("\nvT=\n");
            if (vT != null) {
                sb.append(FormatArray.toString(vT, "%.4e"));
            }
            return sb.toString();
        }

    }

    /**
     * performs SVD on matrix a and if fails to converge, performs SVD on
     * a*a^T and a^T*a separately to get the factorization components for a.
     * <pre>
          SVD(A).U == SVD(A^T).V == SVD(AA^T).U == SVD(AA^T).V

          SVD(A).V == SVD(A^T).U == SVD(A^TA).V == SVD(A^TA).U

          SVD(A) eigenvalues are the same as sqrt( SVD(AA^T) eigenvalues )
              and sqrt( SVD(A^TA) eigenvalues )
       </pre>
     @param a matrix
     @return result of SVD of a or aT*a with re-scaled diagonal
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static SVDProducts performSVD(double[][] a) throws NotConvergedException {
        return performSVD(new DenseMatrix(a));
    }

    /**
     * performs SVD on matrix a and if fails to converge, performs SVD on
     * a*a^T and a^T*a separately to get the factorization components for a.
     * <pre>
          SVD(A).U == SVD(A^T).V == SVD(AA^T).U == SVD(AA^T).V

          SVD(A).V == SVD(A^T).U == SVD(A^TA).V == SVD(A^TA).U

          SVD(A) eigenvalues are the same as sqrt( SVD(AA^T) eigenvalues )
              and sqrt( SVD(A^TA) eigenvalues )
       </pre>
     @param a matrix
     @return result of SVD of a or aT*a with re-scaled diagonal
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static SVDProducts performSVD(DenseMatrix a) throws NotConvergedException {
        SVD svd;
        DenseMatrix u = null;
        DenseMatrix vT = null;
        double[] sDiag = null;
        try {
            svd = SVD.factorize(a);
            vT = svd.getVt();
            u = svd.getU();
            sDiag = svd.getS();
        } catch (NotConvergedException e) {
            double[][] _a = MatrixUtil.convertToRowMajor(a);
            double[][] aTa = MatrixUtil.multiply(MatrixUtil.transpose(_a), _a);
            double[][] aaT = MatrixUtil.multiply(_a, MatrixUtil.transpose(_a));
            //SVD(A).U == SVD(A^T).V == SVD(AA^T).U == SVD(AA^T).V
            //SVD(A).V == SVD(A^T).U == SVD(A^TA).V == SVD(A^TA).U
            //SVD(A) eigenvalues are the same as sqrt( SVD(AA^T) eigenvalues )
            //    and sqrt( SVD(A^TA) eigenvalues )
            svd = SVD.factorize(new DenseMatrix(aTa));
            vT = svd.getVt();
            sDiag = svd.getS();
            for (int i = 0; i < sDiag.length; ++i) {
                if (sDiag[i] > 0) {
                    sDiag[i] = Math.sqrt(sDiag[i]);
                }
            }

            svd = SVD.factorize(new DenseMatrix(aaT));
            u = svd.getU();
        }
        SVDProducts out = new SVDProducts();
        out.u = (u != null) ? MatrixUtil.convertToRowMajor(u) : null;
        out.vT = (vT != null) ? MatrixUtil.convertToRowMajor(vT) : null;
        out.s = sDiag;
        out.rank = rank(svd);
        return out;
    }

    /**
     * use the MTJ arpack to partially solve symmetric eigensystems for eigenvalues.
     * You can get the eigenvectors using
     *
     <pre>
     see https://github.com/fommil/matrix-toolkits-java/blob/master/src/test/java/no/uib/cipr/matrix/sparse/ArpackSymTest.java
     Internally, this method does this:
     For example, to calculate the smallest 2 eigenvectors of matrix A,
     (1) calculates A^T * A to get the symmetric matrix:
         LinkedSparseMatrix A = ...;
         LinkedSparseMatrix AtA = A.transAmult(A, new LinkedSparseMatrix(A.numColumns(), A.numColumns()));
     (2) Then uses the Arpack solver with the enum for smallest eigenvalues.
         ArpackSym solver = new ArpackSym(AtA);
         Map Double, DenseVectorSub  results = solver.solve(2, ArpackSym.Ritz.SA);
     </pre>
     @param A
     @param nEigenValues
     @param ritz
        <pre>
        ArpackSym.Ritz types:
        BE : compute NEV eigenvalues, half from each end of the spectrum
        LA : compute the NEV largest (algebraic) eigenvalues.
        LM : compute the NEV largest (in magnitude) eigenvalues.
        SA : compute the NEV smallest (algebraic) eigenvalues.
        SM : compute the NEV smallest (in magnitude) eigenvalues.
        </pre>
     @return nEigenValues eigenvalues solved for the given sparse matrix A and ritz enum.
     */
     public static Map<Double, DenseVectorSub> sparseEigen(Matrix A, int nEigenValues, ArpackSym.Ritz ritz) {
         Matrix AtA;
         if (A instanceof LinkedSparseMatrix) {
             AtA = A.transAmult(A, new LinkedSparseMatrix(A.numColumns(), A.numColumns()));
         } else if (A instanceof FlexCompColMatrix) {
             AtA = A.transAmult(A, new FlexCompColMatrix(A.numColumns(), A.numColumns()));
         } else {
             throw new IllegalArgumentException("Incomplete Method Error.  This method currently can only handle" +
                     " LinkedSparseMatrix and FlexCompColMatrix");
         }

         ArpackSym solver = new ArpackSym(AtA);

         return solver.solve(nEigenValues, ritz);
     }

    /**
     * create matrix A^T*A then perform SVD on it.  NOTE that the singular values
     * returned in S will have the square of values of SVD(A).s.
     @param a the matrix a (internally, a^T*a will be calculated and used)
     @return result of SVD of aT*a
     * @throws NotConvergedException thrown by MTJ if SVD did not converge
     */
    public static SVDProducts performSVDATransposeA(double[][] a) throws NotConvergedException {
        double[][] aTa = MatrixUtil.createATransposedTimesA(a);
        return performSVDATransposeA(new DenseMatrix(aTa));
    }

    /**
     *
     @param aTa
     @return
     * @throws NotConvergedException
     */
    public static SVDProducts performSVDATransposeA(DenseMatrix aTa) throws NotConvergedException {

        SVD svd;
        DenseMatrix u = null;
        DenseMatrix vT = null;
        double[] sDiag = null;

        svd = SVD.factorize(aTa);
        vT = svd.getVt();
        u = svd.getU();
        sDiag = svd.getS();

        SVDProducts out = new SVDProducts();
        out.u = MatrixUtil.convertToRowMajor(u);
        out.vT = MatrixUtil.convertToRowMajor(vT);
        out.s = sDiag;
        return out;
    }

    /**
     * perform QR decomposition (a.k.a. Francis algorithm, a.k.a. Francis QR step)
     * on matrix a using the MTJ library.
     * It's a sophisticated version of the "power method".
     * A = Q*R of an orthonormal matrix Q and an upper triangular matrix R.
     * The columns of Q are the eigenvectors of A.
     * A*Q = Q * diag(eigenvalues of A).
     *
     * <pre>
     * To calculate the product of eigenvalue or singular values using
     * QR, see section "Connection to a determinant or a product of eigenvalues"
     * in wikipedia:
     * https://en.wikipedia.org/wiki/QR_decomposition
     * </pre>
     *
     @param a a square or rectangular matrix with independent columns.
     @return resulting QR decomposition of a
     */
    public static QR performQRDecomposition(double[][] a) {
        QR qr = QR.factorize(new DenseMatrix(a));
        return qr;
    }

    /**
     * given data points xy, want to create a matrix usable to transform
     * the data points by scaling and translation so that:
        a) points are translated so that their centroid is at the origin.
        b) points are scaled so that the average distance from the
           origin is sqrt(2).
       Can use the transformation matrix with dot operator: Misc.multiply(xy, tMatrix).
     @param xy points in format [n X 2]
     @return a matrix for use for canonical transformation of the points.
     * the format of the result is
     * <pre>
     *  t[0] = new double[]{scale,       0,     -centroidX*scale};
        t[1] = new double[]{0,           scale, -centroidY*scale};
       </pre>
     */
    public static double[][] calculateNormalizationMatrix2X3(double[][] xy) {

        int nRows = xy.length;

        if (nRows == 0) {
            System.err.println("xy.length must be > 0");
            return new double[2][0];
        }

        double cen0 = 0;
        double cen1 = 0;
        for (int i = 0; i < nRows; ++i) {
            cen0 += xy[i][0];
            cen1 += xy[i][1];
        }
        cen0 /= (double)nRows;
        cen1 /= (double)nRows;

        double mean = 0;
        for (int i = 0; i < nRows; ++i) {
            double diffX = xy[i][0] - cen0;
            double diffY = xy[i][1] - cen1;
            double dist = Math.sqrt((diffX * diffX) + (diffY * diffY));
            mean += dist;
        }
        mean /= (double)nRows;

        /*
        mean * factor = sqrt(2)
        */
        double scale = Math.sqrt(2)/mean;

        double[][] t = new double[2][];
        t[0] = new double[]{scale,       0,     -cen0*scale};
        t[1] = new double[]{0,           scale, -cen1*scale};

        return t;
    }

    /**
     *
     @param a
     @return
     */
    public static float[][] copy(float[][] a) {

        float[][] c = new float[a.length][a[0].length];

        for (int i = 0; i < a.length; ++i) {
            c[i] = new float[a[0].length];
            System.arraycopy(a[i], 0, c[i], 0, a[0].length);
        }

        return c;
    }
    public static int[][] copy(int[][] a) {

        int[][] c = new int[a.length][a[0].length];

        for (int i = 0; i < a.length; ++i) {
            System.arraycopy(a[i], 0, c[i], 0, a[0].length);
        }

        return c;
    }

    public static boolean equals(int[][] a, int[][] b) {
        if (a.length != b.length) return false;
        if (a[0].length != b[0].length) return false;
        for (int i = 0; i < a.length; ++i) {
            if (!Arrays.equals(a[i], b[i])) {
                return false;
            }
        }
        return true;
    }

    /**
     *
     @param a
     @return
     */
    public static double[][] copy(double[][] a) {

        double[][] m = new double[a.length][];

        for (int i = 0; i < a.length; ++i) {
            int n0 = a[i].length;
            m[i] = new double[n0];
            System.arraycopy(a[i], 0, m[i], 0, n0);
        }

        return m;
    }


    /**
     *
     @param a
     @return
     */
    public static float[] copy(double[] a) {

        float[] m = new float[a.length];

        for (int i = 0; i < a.length; ++i) {
            m[i] = (float)a[i];
        }

        return m;
    }

    public static float[] copyLongToFloat(long[] a) {

        float[] m = new float[a.length];

        for (int i = 0; i < a.length; ++i) {
            m[i] = (float)a[i];
        }

        return m;
    }

    public static double[] copyLongToDouble(long[] a) {

        double[] m = new double[a.length];

        for (int i = 0; i < a.length; ++i) {
            m[i] = (float)a[i];
        }

        return m;
    }

    /**
     *
     @param source
     @param destination
     */
    public static void copy(double[][] source, double[][] destination) {

        int m = source.length;
        int n = source[0].length;

        if (destination.length != m) {
            throw new IllegalArgumentException("destination.length must equal source.length");
        }
        if (destination[0].length != n) {
            throw new IllegalArgumentException("destination[0].length must equal source[0].length");
        }

        for (int i = 0; i < source.length; ++i) {
            System.arraycopy(source[i], 0, destination[i], 0, n);
        }
    }

    /**
     *
     @param a matrix
     @param row0 beginning index, inclusive
     @param row1 end index, inclusive
     @param col0 beginning index, inclusive
     @param col1 end index, inclusive
     @return sub matrix
     */
    public static double[][] copySubMatrix(double[][] a, int row0, int row1, int col0, int col1) {

        int nr2 = row1 - row0 + 1;
        int nc2 = col1 - col0 + 1;

        double[][] m = new double[nr2][];
        int i, j;
        for (i = 0; i < nr2; ++i) {
            m[i] = new double[nc2];
            System.arraycopy(a[row0 + i], col0, m[i], 0, nc2);
        }

        return m;
    }

    /**
     * copy the section of matrix a from row0 to row1 (inclusive) and
     * col0 to col1 (inclusive) into output matrix out which must
     * be size (row1-row0+1) X (col1-col0+1)
     @param a matrix
     @param row0 beginning index, inclusive
     @param row1 end index, inclusive
     @param col0 beginning index, inclusive
     @param col1 end index, inclusive
     @param out output matrix to hold the copied section
     */
    public static void copySubMatrix(double[][] a, int row0, int row1,
        int col0, int col1, double[][] out) {

        int nr2 = row1 - row0 + 1;
        int nc2 = col1 - col0 + 1;

        if (out.length != nr2 || out[0].length != nc2) {
            throw new IllegalArgumentException("out dimensions must be "
                    + " row1 - row0 + 1 X col1 - col0 + 1");
        }

        int i, j;
        for (i = 0; i < nr2; ++i) {
            System.arraycopy(a[row0 + i], col0, out[i], 0, nc2);
        }
    }

    /**
     *
     @param a matrix
     @param col index of column to extract
     @return one dimensional array holding the column col of a
     */
    public static double[] extractColumn(double[][] a, int col) {

        double[] m = new double[a.length];
        extractColumn(a, col, m);
        return m;
    }

    /**
     *
     @param a matrix a
     @param col index of column to extract
     @param out one dimensional array holding the column col of a
     */
    public static void  extractColumn(double[][] a, int col, double[] out) {

        if (out.length != a.length) {
            throw new IllegalArgumentException("out.length must equal a.length");
        }

        int i;
        for (i = 0; i < a.length; ++i) {
            out[i] = a[i][col];
        }
    }

      /**
     * using cofactors and minors of the matrix, return the determinant.
     * in practice one can use any row as the primary set of cofactors or
     * any column.  this method may be optimized in the future, but for now,
     * uses the first column as the cofactors.
     *
     * e.g.    | 1  -5  2 |         | 3 4 |         | 7 4 |         | 7 3 |
     *         | 7   3  4 |  =  1 * | 1 5 |  +  5 * | 2 5 |  +  2 * | 2 1 |  = 11
        + 135 + 2 = 148
     *         | 2   1  5 |
     * <pre>
     * Note that det(a) = 0 shows that matrix a is a singular matrix and is not
     * invertible.
     * </pre>
     @param a a square matrix
     @return the determinant of matrix a
     */
    public static double determinant(Matrix a) {

        double[][] ma = no.uib.cipr.matrix.Matrices.getArray(a);

        return determinant(ma);
    }

    /**
     * using cofactors and minors of the matrix, return the determinant.
     * in practice one can use any row as the primary set of cofactors or
     * any column.  this method may be optimized in the future, but for now,
     * uses the first column as the cofactors.
     *
     * e.g.    | 1  -5  2 |         | 3 4 |         | 7 4 |         | 7 3 |
     *         | 7   3  4 |  =  1 * | 1 5 |  +  5 * | 2 5 |  +  2 * | 2 1 |  = 11 + 135 + 2 = 148
     *         | 2   1  5 |
     * <pre>
     * Note that det(a) = 0 shows that matrix a is a singular matrix and is not
     * invertible.
     * </pre>
     @param a a square matrix
     @return determinant of a
     */
    public static double determinant(double[][] a) {

        if (a == null || a.length == 0) {
            throw new IllegalArgumentException("matrix a cannot be null or empty");
        }
        if (a.length != a[0].length) {
            throw new IllegalArgumentException("matrix a must be square");
        }
        if (a.length == 1) {
            return roundToMachineTol(a[0][0]);
        } else if (a.length == 2) {
            double s = ( a[0][0]*a[1][1] ) - ( a[0][1]*a[1][0] );
            return roundToMachineTol(s);
        } else {
            double s = 0.0;
            // use 1st row as cofactors and minors
            for (int i = 0; i < a.length; i++) {

                double[][] n = copyExcept(a, i, 0);

                double tmp = a[i][0] * determinant(n);

                if ((i & 1) == 0) {
                    s +=  tmp;
                } else {
                    s -=  tmp;
                }
            }
            return roundToMachineTol(s);
        }
    }


    /**
     * using cofactors and minors of the matrix, return the determinant.
     * in practice one can use any row as the primary set of cofactors or
     * any column.  this method may be optimized in the future, but for now,
     * uses the first column as the cofactors.
     *
     * e.g.    | 1  -5  2 |         | 3 4 |         | 7 4 |         | 7 3 |
     *         | 7   3  4 |  =  1 * | 1 5 |  +  5 * | 2 5 |  +  2 * | 2 1 |  = 11 + 135 + 2 = 148
     *         | 2   1  5 |
     * <pre>
     * Note that det(a) = 0 shows that matrix a is a singular matrix and is not
     * invertible.
     * </pre>
     @param a a square matrix
     @return determinant of a
     */
    public static float determinant(float[][] a) {

        if (a == null || a.length == 0) {
            throw new IllegalArgumentException("matrix a cannot be null or empty");
        }
        if (a.length != a[0].length) {
            throw new IllegalArgumentException("matrix a must be square");
        }
        if (a.length == 1) {
            return roundToMachineTol(a[0][0]);
        } else if (a.length == 2) {
            float s = ( a[0][0]*a[1][1] ) - ( a[0][1]*a[1][0] );
            return roundToMachineTol(s);
        } else {

            float[][] n = new float[a.length - 1][a.length - 1];

            float s = 0.0f;
            // use 1st row as cofactors and minors
            for (int i = 0; i < a.length; i++) {

                copyExcept(a, i, 0, n);

                float tmp = a[i][0] * determinant(n);

                if ((i & 1) == 0) {
                    s +=  tmp;
                } else {
                    s -=  tmp;
                }
            }
            return roundToMachineTol(s);
        }
    }

    private static double roundToMachineTol(double a) {
        return Math.round(a * 1E11)/1E11;
    }
    private static float roundToMachineTol(float a) {
        return (float)(Math.round(a * 1E11)/1E11);
    }

    /**
     * calculate the determinant of a using the diagonal of U from the
     * LU decomposition.
     *
     @param a a square matrix
     @return determinant of a
     */
    public static double determinantFromLU(double[][] a) {

        if (!isSquare(a)) {
            throw new IllegalArgumentException("matrix a must be square");
        }
        switch (a.length) {
            case 0:
                throw new IllegalArgumentException("matrix a length must be larger than 0");
            case 1:
                return a[0][0];
            case 2: {
                double s = ( a[0][0]*a[1][1] ) - ( a[0][1]*a[1][0] );
                return s;
            }
            default: {
                LU lu = LinearEquations.LUDecomposition(a);
                double d = 1;
                int i;
                for (i = 0; i < lu.u.length; ++i) {
                    d *= lu.u[i][i];
                }
                return d;
            }
        }
    }


    /**
     * create copy of matrix m except row and col and store it in out matrix
     @param m matrix
     @param col column
     @param row row
     @return extracted copy of matrix minus row and col
     */
    private static void copyExcept(float[][] m, int col, int row, float[][] out) {

        if (out.length != m.length - 1 || out[0].length != m.length - 1) {
            throw new IllegalArgumentException("out must be [m.length-1 X m.length-1]");
        }

        int nr = 0;
        int nc = 0;

        for (int mCol = 0; mCol < m.length; mCol++) {
            if (mCol == col) {
                continue;
            }

            nr = 0;
            for (int mRow = 0; mRow < m[0].length; mRow++) {
                if (mRow == row) {
                    continue;
                }

                out[nc][nr] = m[mCol][mRow];
                nr++;
            }
            nc++;
        }
    }

    /**
     * create copy of matrix m except row and col
     @param m matrix
     @param col column
     @param row row
     @param out extracted copy of matrix minus row and col
     */
    private static void copyExcept(double[][] m, int col, int row, double[][] out) {

        if (out.length != m.length - 1 || out[0].length != m.length - 1) {
            throw new IllegalArgumentException("out must be [m.length-1 X m.length-1]");
        }

        int nr = 0;
        int nc = 0;

        for (int mCol = 0; mCol < m.length; mCol++) {
            if (mCol == col) {
                continue;
            }

            nr = 0;
            for (int mRow = 0; mRow < m[0].length; mRow++) {
                if (mRow == row) {
                    continue;
                }

                out[nc][nr] = m[mCol][mRow];
                nr++;
            }
            nc++;
        }

    }

    /**
     * create copy of matrix m except row and col
     @param m matrix
     @param col column
     @param row row
     @return extracted copy of matrix minus row and col
     */
    private static double[][] copyExcept(double[][] m, int col, int row) {

        double[][] n = new double[m.length - 1][m.length - 1];

        int nr = 0;
        int nc = 0;

        for (int mCol = 0; mCol < m.length; mCol++) {
            if (mCol == col) {
                continue;
            }

            n[nc] = new double[m.length - 1];

            nr = 0;
            for (int mRow = 0; mRow < m[0].length; mRow++) {
                if (mRow == row) {
                    continue;
                }

                n[nc][nr] = m[mCol][mRow];
                nr++;
            }
            nc++;
        }

        return n;
    }

    /**
     * constructs the 3x3 skew-symmetric matrices for use in cross products,
     * notation is [v]_x.
     * the operator is also called "hat operator".
     * v cross product with w is v X w = [v]_x * w.
     * It’s individual terms are a_j_i = -a_i_j.
       <pre>
       |    0   -v[2]   v[1] |
       |  v[2]    0    -v[0] |
       | -v[1]  v[0]      0  |

       Note that the skew symmetric matrix equals its own negative, i.e. A^T = -A.
       </pre>
     @param v vector
     @return skew-symmetric matrix of vector v
     */
    public static double[][] skewSymmetric(double[] v) {
        if (v.length != 3) {
            throw new IllegalArgumentException("v.length must be 3");
        }

        double[][] out = new double[3][3];
        skewSymmetric(v, out);
        return out;
    }

    /**
     * constructs the 3x3 skew-symmetric matrices for use in cross products,
     * notation is [v]_x.
     * v cross product with w is v X w = [v]_x * w.
     * It’s individual terms are a_j_i = -a_i_j.
       <pre>
       |    0   -v[2]   v[1] |
       |  v[2]    0    -v[0] |
       | -v[1]  v[0]      0  |

       Note that the skew symmetric matrix equals its own negative, i.e. A^T = -A.
       </pre>
     * the operator is also called "hat operator".
     Its opposite is the anti-symmetric matrix and is equal to -1 * skewSymmetric.
     <pre>
     some of its properties:
        let [u]_x be the skew symmetric of u and [v]_x be the skew-symmetric of v.
     where u and v are column vectors.

     -[u]_x * v = [v]_x * u

     -[u]_x * u = 0

     [u]_x = [ -[u]_x ]^T

     -[u]_x * -[v]_x = (u dot v) * I + v * u^T where '*' is matrix multiplication as usual.
     [u]_x * [v]_x = (u dot v) * I + v * u^T

     [u]_x * [v]_x - [v]_x * [u]_x = v * u^T - u*v^T = [u cross v]]_x

     u * v^T * -[w]_x + -[w]_x * v * u^T = -[u cross (v cross w)]_x
           where the later, triple cross product is used in making the Fundamental Matrix
           of photogrammetry, for example.

     </pre>
     @param v direction vector
     @param out output skew-symmetric matrix for v
     */
    public static void skewSymmetric(double[] v, double[][] out) {
        if (v.length != 3) {
            throw new IllegalArgumentException("v.length must be 3");
        }
        if (out.length != 3 || out[0].length != 3) {
            throw new IllegalArgumentException("out must be 3X3");
        }

        /*
        out[0] = new double[]{0,     -v[2],  v[1]};
        out[1] = new double[]{v[2],     0,  -v[0]};
        out[2] = new double[]{-v[1], v[0],    0};
        */

        out[0][0] = 0;
        out[0][1] = -v[2];
        out[0][2] = v[1];
        out[1][0] = v[2];
        out[1][1] = 0;
        out[1][2] = -v[0];
        out[2][0] = -v[1];
        out[2][1] = v[0];
        out[2][2] = 0;

    }

    /**
     * extract the vector from a skew-symmetric matrix.  the operation is called vee operator.
     @param vHat skew-symmetric matrix
     @return extracted vector
     */
    public static double[] extractVectorFromSkewSymmetric(double[][] vHat) {
        double[] out = new double[3];
        out[0] = 0.5*(vHat[2][1] - vHat[1][2]);
        out[1] = 0.5*(vHat[0][2] - vHat[2][0]);
        out[2] = 0.5*(vHat[1][0] - vHat[0][1]);
        return out;
    }

    /**
     *
     @param p0
     @param p1
     @return
     */
    public static double[] crossProduct(double[] p0, double[] p1) {
        if (p0.length != 3 || p1.length != 3) {
            throw new IllegalArgumentException("expecting p0 and p1 lengths to be 3");
        }
        double[] out = new double[3];
        out[0] = p0[1]*p1[2] - p0[2]*p1[1];
        out[1] = p0[2]*p1[0] - p0[0]*p1[2];
        out[2] = p0[0]*p1[1] - p0[1]*p1[0];
        return out;
    }

    /**
     *
     @param p0
     @param p1
     @param eps
     @return
     */
    public static boolean areColinear(double[] p0, double[] p1, double eps) {
        if (p0.length != 3 || p1.length != 3) {
            throw new IllegalArgumentException("expecting p0 and p1 lengths to be 3");
        }
        double[] cp = crossProduct(p0, p1);
        double norm = lPSum(cp, 2);
        return Math.abs(norm) < eps;
    }

    /**
     *
     @param p0
     @param p1
     @param p2
     @param eps
     @return
     */
    public static boolean areColinear(double[] p0, double[] p1, double[] p2, double eps) {
        if (p0.length != 3 || p1.length != 3 || p2.length != 3) {
            throw new IllegalArgumentException("expecting p0, p1 and p2 lengths to be 3");
        }
        //adapted from code by Peter Kovesis for function
        // iscolinear.m.  r =  norm(cross(p2-p1, p3-p1)) < eps
        // https://www.peterkovesi.com/matlabfns/
        double[] p21 = MatrixUtil.subtract(p1, p0);
        double[] p31 = MatrixUtil.subtract(p2, p0);
        double[] cp = crossProduct(p21, p31);
        double norm = lPSum(cp, 2);
        return Math.abs(norm) < eps;
    }

    /**
     * normalize vector v by euclidean, that is the square root of the sum of
     * its squared components.  notation is sometimes ||v||_2.
     @param v array
     @return array v divided by the LP2 norm of v.
     */
    public static double[] normalizeL2(double[] v) {
        return normalizeLP(v, 2);
    }

    /**
     * normalize each column of matrix a by the square root of the sum of
     * its squared components. ||v||_2 for each column...
     @param a matrix
     */
    public static void normalizeColumnsL2(double[][] a) {
        double sum;
        int j;
        for (int i = 0; i < a[0].length; ++i) {
            sum = 0;
            for (j = 0; j < a.length; ++j) {
                sum += (a[j][i]*a[j][i]);
            }
            sum = Math.sqrt(sum);
            for (j = 0; j < a.length; ++j) {
                if (sum > 0.) {
                    a[j][i] /= sum;
                }
            }
        }
    }

    /**
     * normalize each row of matrix a by the square root of the sum of
     * its squared components. ||v||_2 for each row...
     @param a matrix
     */
    public static void normalizeRowsL2(double[][] a) {
        double sum;
        int j;
        for (int i = 0; i < a.length; ++i) {
            sum = 0;
            for (j = 0; j < a[0].length; ++j) {
                sum += (a[i][j]*a[i][j]);
            }
            sum = Math.sqrt(sum);
            for (j = 0; j < a[0].length; ++j) {
                if (sum > 0.) {
                    a[i][j] /= sum;
                }
            }
        }
    }

    /**
     * calculate the condition number as the largest singular value divided
     * by the singular value for i==(rank-1) of A where the singular values are
     * found using the SVD.
     <pre>
     from https://blogs.mathworks.com/cleve/2017/07/17/what-is-the-condition-number-of-a-matrix/
     A condition number for a matrix and computational task measures how
     sensitive the answer is to perturbations in the input data and to roundoff
     errors made during the solution process....If a matrix is singular,
     then its condition number is infinite.
     ...(A large condition number means that the matrix is close to being singular).

     also see Section 9.2 of Strang's "Introduction to Linear Algebra".
     </pre>
     *
     @param a matrix
     @return condition number
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static double conditionNumber(double[][] a) throws NotConvergedException {

        SVDProducts svd = performSVD(a);
        int rm1 = -1;
        for (int i = 0; i < svd.s.length; ++i) {
            if (svd.s[i] > 0) {
                rm1 = i;
            }
        }
        if (rm1 == -1) {
            throw new IllegalStateException("all singular values of a are 0");
        }
        double c = svd.s[0]/svd.s[rm1];

        return c;
    }

    /**
     pre- and post- multiply using the diagonal matrices D_1 and D_2, depending upon the problem.
         example use of a pre-conditioned matrix A:
           solve equation A_hat * x_hat = b_hat

         For pre-multiplication:
                A_hat = D_1 * A * D_2
                b_hat = (D_2)^-1 * b
                x_hat = (D_2)^-1 * x  (this isn't done when solving for x in A*x=b)
         For post-multiplication:
                x = x_hat * D_2
     */
    public static class ConditionedMatrix {
        /**
         * A_hat = D_1 * A * D_2
            where D_1 is diagonal(d1)
            and D_2 is diagonal(d2)
         */
        public double[][] aHat;

        /**
         * This matrix diagonal can be used to pre-condition (pre-multiply) vector b in A*x=b:
         * b_hat = D_1 * b
            where D_1 is diagonal(d1)
         */
        public double[] d1;

        /**
         * This matrix diagonal can be used to post-multiply x_hat in A_hat * x_hat = b_hat.
         * x = x_hat * D_2
         * where D_w is diagonal(d2)
         */
        public double[] d2;

    }

    /**
     * apply Ruiz scaling as "equilibration" for matrix a as a linear system in order to
     * pre-condition it before use.
     * Note: before using this method, you may want to apply permutations following
     * Duff And Koster 1999 and HSL 2000 routine MC64.
     <pre>
     References:
     Ruiz 2001 "A Scaling Algorithm to Equilibrate Both Rows and Columns Norms in Matrices"
     iterative algorithm that has fast linear convergence and assymptotic rate of 1/2.

     Liu, Paul. "An exploration of matrix equilibration." (2015).
     "equilibration" of a linear system applies scaling and permutations
     to an ill-conditioned matrix as a pre-processing step to improve
     the conditioning.
     The Ruiz algorithm and the Liu algorithm both attempt to minimize the condition number of a
     matrix by scaling it by the infinity-norm of each row and column of a matrix to value 1.
     </pre>
     <pre>
     A_hat = the scaled matrix = D_1 * A * D_2

     rewrite:  A_hat =  D_1 * A * D_2
          as:  A_hat * (D_2)^-1 = D_1 * A

     For use in the example solving for x in A * x = b:
        A_hat * x_hat = b_hat
        For pre-multiplication:
            A_hat =  D_1 * A * D_2
            b_hat = D_1 * b
            x_hat = (D_2)^-1 * x  (this isn't done when solving for x in A*x=b)
        For post-multiplication:
            x = x_hat * D_2

     </pre>
     TODO: see LAPack xGEBAL
        can access it in MTJ using LAPack.getInstance()
        see more in com.github.fommil.netlib.
     @param a an m X n matrix.
     @return the scaled matrix a_hat and the 2 matrices used to scale the matrix.  the
       matrices can be used to pre-multiply and post-multiply other data structures in a problem.
     */
    public static ConditionedMatrix preconditionScaling(double[][] a) {

        /*
        rows r_i = a_i^T and of dimension n+1 where i=1,...m
        cols c_j = a_j   and of dimension m?+1 where j=1...n
        D_R is mxm diagonal matrix = diag( sqrt( ||r_i||_inf ) for i=1,...m
        D_C is nxn diagonal matrix = diag( sqrt( ||c_j||_inf ) for j=1,...n
          where ||.||_inf is the infinity norm and is the maximum of absolute values of the internal set
        if any inf-norm is 0, the inf-norm is replaced by value 1
        */

        int m = a.length;
        int n = a[0].length;
        double eps = 1E-8;

        /*
        alg:
          set A_hat(0) = A
              D_1(0) = I_m
              D_2(0) = I_n
          for k=0,1,2,... until convergence do:
            D_R = diag( sqrt( ||r_i(k)||_inf ) for i=1,...m
            D_C = diag( sqrt( ||c_j(k)||_inf ) for j=1,...n
            A_hat(k+1) = (D_R)^-1 * A_hat(k) * (D_C)^-1
            D_1(k+1) = D_1(k) * (D_R)^-1
            D_2(k+1) = D_2(k) * (D_C)^-1

          convergence:
             max_{1<=i<=m} ( | (1 - ||r_i(k)||_inf) | ) <= eps
             and max_{1<=j<=n} ( | (1 - ||c_j(k)||_inf) | ) <= eps

         rewriting to use prev and current:
           aHatPrev = A
           d1Prev = I_m
           d2Prev = I_n
           while (true) {
               dR = diag( sqrt( ||r_i in aHatPrev||_inf ) for i=1,...m
               dC = diag( sqrt( ||c_j in aHatPrev||_inf ) for j=1,...j
               aHatCurr = (dR)^-1 * aHatPrev * (dC)^-1
               d1Curr = d1Prev * (dR)^-1
               d2Curr = d2Prev * (dC)^-1

               if (conv conditions) {
                   break;
               }
               aHatPrev = aHatCurr
               d1Prev = d1Curr;
               d2Prev = d2Curr;
           }
        */
        double[][] aHatCurr = MatrixUtil.zeros(m, n);
        double[][] tmp = MatrixUtil.zeros(m, n);
        double[] d1Curr = new double[m];
        double[] d2Curr = new double[n];
        double[] rInfNorm = new double[m];
        double[] cInfNorm = new double[n];
        double dR;
        double dC;
        double[] dRInv = new double[n];
        double[] dCInv = new double[n];

        // init the prev arrays
        double[][] aHatPrev = MatrixUtil.copy(a);
        double[] d1Prev = new double[m];
        Arrays.fill(d1Prev, 1);
        double[] d2Prev = new double[n];
        Arrays.fill(d2Prev, 1);

        int i;
        double s;
        double t;
        final int maxIter = 100;
        int nIter = 0;

        /*
        System.out.printf("  %d) aHat=\n%s\n", nIter, FormatArray.toString(a, "%.5e"));
        try {
            System.out.printf("    %.3e\n", MatrixUtil.conditionNumber(a));
        } catch (NotConvergedException e) {
            e.printStackTrace();
        }*/

        while (nIter < maxIter) {
            calculateInfNormForRows(aHatCurr, rInfNorm);
            calculateInfNormForColumns(aHatCurr, cInfNorm);
            for (i = 0; i < m; ++i) {
                dR = Math.sqrt(rInfNorm[i]);
                if (dR < 1E-11) {// eq 0
                    dRInv[i] = 1;
                } else {
                    dRInv[i] = 1. / dR;
                }
            }
            for (i = 0; i < n; ++i) {
                dC = Math.sqrt(cInfNorm[i]);
                if (dC < 1E-11) {// eq 0
                    dCInv[i] = 1;
                } else {
                    dCInv[i] = 1. / dC;
                }
            }

            //aHatCurr = (dR)^-1 * aHatPrev * (dC)^-1
            MatrixUtil.multiplyDiagonalByMatrix(dRInv, aHatPrev, tmp);
            MatrixUtil.multiplyByDiagonal(tmp, dCInv, aHatCurr);

            /*
            System.out.printf("  %d)\n", nIter);
            System.out.printf("      rInfNorm=%s\n", FormatArray.toString(rInfNorm, "%.5e"));
            System.out.printf("      cInfNorm=%s\n", FormatArray.toString(cInfNorm, "%.5e"));
            System.out.printf("      dRInv=%s\n", FormatArray.toString(dRInv, "%.5e"));
            System.out.printf("      dCInv=%s\n", FormatArray.toString(dCInv, "%.5e"));
            System.out.printf("      aHat=\n%s", FormatArray.toString(aHatCurr, "%.5e"));
            try {
                System.out.printf("    c.n.=%.3e\n", MatrixUtil.conditionNumber(aHatCurr));
            } catch (NotConvergedException e) {
                e.printStackTrace();
            }
            */

            //d1Curr = d1Prev * (dR)^-1
            for (i = 0; i < m; ++i) {
                d1Curr[i] = d1Prev[i] * dRInv[i];
            }
            //d2Curr = d2Prev * (dC)^-1
            for (i = 0; i < n; ++i) {
                d2Curr[i] = d2Prev[i] * dCInv[i];
            }

            ++nIter;
            // test convergence
            s = Double.NEGATIVE_INFINITY;
            for (i = 0; i < m; ++i) {
                t = Math.abs(1 - rInfNorm[i]);
                if (t > s) {
                    s = t;
                }
            }
            System.out.printf("      s=%.3e\n", s);
            if (s <= eps) {
                // test column convergence
                s = Double.NEGATIVE_INFINITY;
                for (i = 0; i < n; ++i) {
                    t = Math.abs(1 - cInfNorm[i]);
                    if (t > s) {
                        s = t;
                    }
                }
                if (s <= eps) {
                    break;
                }
            }
            //aHatPrev = aHatCurr
            MatrixUtil.copy(aHatCurr, aHatPrev);
            //d1Prev = d1Curr;
            System.arraycopy(d1Curr,0, d1Prev, 0, d1Curr.length);
            //d2Prev = d2Curr;
            System.arraycopy(d2Curr,0, d2Prev, 0, d2Curr.length);
        }

        ConditionedMatrix c = new ConditionedMatrix();
        c.aHat = aHatCurr;
        c.d1 = d1Curr;
        c.d2 = d2Curr;
        return c;
    }

    /**
     * calculate the infinity norm for each column of matrix a where the infinity norm
     * for each col as a set is the maximum absolute value for that col.
     @param a an m x n matrix
     @param output array of size a.length to be populated with the infinity norms for each column of a
     */
    public static void calculateInfNormForColumns(double[][] a, double[] output) {
        int n = a[0].length;
        if (output.length != n) {
            throw new IllegalArgumentException("output.length must equal a[0].length");
        }
        Arrays.fill(output, Double.NEGATIVE_INFINITY);
        double c;
        int i, j;
        for (j = 0; j < n; ++j) {
            for (i = 0; i < a.length; ++i) {
                c = Math.abs(a[i][j]);
                if (c > output[j]) {
                    output[j] = c;
                }
            }
        }
    }

    /**
     * calculate the infinity norm for each row of matrix a where the infinity norm
     * for each row as a set is the maximum absolute value for that row.
     @param a an m x n matrix
     @param output array of size a.length to be populated with the infinity norms for each row of a
     */
    public static void calculateInfNormForRows(double[][] a, double[] output) {
        int m = a.length;
        if (output.length != m) {
            throw new IllegalArgumentException("output.length must equal a.length");
        }
        Arrays.fill(output, Double.NEGATIVE_INFINITY);
        double r;
        int j;
        for (int i = 0; i < m; ++i) {
            for (j = 0; j < a[i].length; ++j) {
                r = Math.abs(a[i][j]);
                if (r > output[i]) {
                    output[i] = r;
                }
            }
        }
    }

    /**
     * summation = the (1/p) power of sum of its (components)^p.
     @param v array v
     @param p the power of the polynomial of v in the sum
     @return he (1/p) power of sum of its (components)^p
     */
    public static double lPSum(double[] v, double p) {
        double sum = 0;
        for (double a : v) {
            sum += Math.pow(a, p);
        }
        sum = Math.pow(sum, 1./p);
        return sum;
    }

    /**
    following the convention used by Matlab
    <pre>
    https://www.mathworks.com/help/matlab/ref/norm.html#bvhji30-3

    For p-norm = 1, the L1-norm is the maximum absolute column sum of the matrix.
    ||X||_1 = max sum for an arg j where (0.lte.j.lte.n-1) sum_(i=0 to n-1) ( |a[i][j] )
    </pre>
    @param a matrix
    @return the maximum absolute column sum of the matrix
    */
    public static int lp1Norm(int[][] a) {
        int maxSum = Integer.MIN_VALUE;
        int i, j, sum;
        for (j = 0; j < a[0].length; ++j) {
            sum = 0;
            for (i = 0; i < a.length; ++i) {
                sum += Math.abs(a[i][j]);
            }
            if (sum > maxSum) {
                maxSum = sum;
            }
        }
        return maxSum;
    }

    /**
     * calculate the Frobenius Norm of matrix a.
     * It's the square root of the sum of squares of each element.
     * It can also be calculated as the square root of the
     * trace of a_conjugate*a, or as
     * the square root of the sums of the squares of the singular
     * values of a.
     * (The later can be compared to the spectral norm which is the largest singular value of a.)
     @param a matrix
     @return the sum of each element of a squared
     */
    public static double frobeniusNorm(double[][] a) {
        double sum = 0;
        int i, j;
        for (i = 0; i < a.length; ++i) {
            for (j = 0; j < a[i].length; ++j) {
                sum += (a[i][j]*a[i][j]);
            }
        }
        return Math.sqrt(sum);
    }

    /**
     * calculate the value of the largest singular value of matrix r.  this is the induced norm
     * which measures what is the maximum of ‖𝐴𝑥‖‖𝑥‖ for any 𝑥≠0 (or, equivalently, the maximum of ‖𝐴𝑥‖ for ‖𝑥‖=1).
     * Note that the frobenius norm can be defined as the square root of the sums of the squares of the singular
     * values of and therefore can be different than the spectral norm.
     * <pre>
     *     reference for documentation:
     *     https://math.stackexchange.com/questions/33083/what-is-the-difference-between-the-frobenius-norm-and-the-2-norm-of-a-matrix
     * </pre>
     @param r matrix
     @return the spectral norm of r, which is the largest singular value of r
     * @throws NotConvergedException thrown by MTJ when SVD could not converge
     */
    public static double spectralNorm(double[][] r) throws NotConvergedException {
        SVDProducts svd = performSVDATransposeA(r);
        double norm = Math.sqrt(svd.s[0]);
        return norm;
    }

    /**
     * Given a symmetric matrix and a nonnegative number eps, find the
     * nearest symmetric positive semidefinite matrices with eigenvalues at least eps.
       a pos def matrix is symmetric and its eigenvalues are all positive, 
       so, there's a unique minimum for quadratic equations, but the nearest
       will be a function of eps.
     * <pre>
     * References:
     * https://nhigham.com/2021/01/26/what-is-the-nearest-positive-semidefinite-matrix/
     * Cheng and Higham, 1998
     *
     * </pre>
     @param a a symmetric matrix
     @param eps a tolerance or error above 0 such as machine precision.  must
     * be greater than or equal to 0.  If it is above 0, this method attempts
     * to return a symmetric positive definite matrix.
     @return a matrix which is the nearest positive semidefinite matrix to a
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static double[][] nearestPositiveSemidefiniteToASymmetric(double[][] a,
        double eps) throws NotConvergedException {

        if (eps < 0) {
            throw new IllegalArgumentException("eps must be .gte. 0");
        }

        if (!isSymmetric(a, eps)) {
            System.err.println("matrix is not symmetric in current form");
        }

        // uses the Frobenius norm as a distance.  notation: ||A||_F

        // spectral decomposition
        //  X = Q * diag(tau_i) * Q^T
        //    where tau_i = {    lambda_i for lambda_i >= eps
        //                  { or eps_i    for lambda_i <  eps
        //        where lambda_i are the eigenvalues of A

        EVD evd = EVD.factorize(new DenseMatrix(a));
        double[][] leftEV = MatrixUtil.convertToRowMajor(evd.getLeftEigenvectors());

        double[] tau = Arrays.copyOf(evd.getRealEigenvalues(),
            evd.getRealEigenvalues().length);

        int i;
        for (i = 0; i < tau.length; ++i) {
            if (tau[i] < eps) {
                tau[i] = eps;
            }
        }

        double[][] aPSD = MatrixUtil.multiplyByDiagonal(leftEV, tau);
        aPSD = MatrixUtil.multiply(aPSD, MatrixUtil.transpose(leftEV));

        /*
        // another solution for X uses polar decomposition B = UH,
        //     where U is orthogonal (U*U^T=I)
        //     and H is symmetric positive definite (H = H^T >= 0)
        //         (with H = Q * diag(lambda_i) * Q^T.
        //         (where Q is svd.v)
        //
        //   since A is symmetric):
        //   answer is X_F = (A + H)/2

        // can use polar decompositon on square matrices:
        double[][] b = MatrixUtil.pointwiseAdd(a, MatrixUtil.transpose(a));
        MatrixUtil.multiply(b, 0.5);

        QH qh = performPolarDecomposition(b);

        //X_F:
        double[][] aPSD2 = MatrixUtil.pointwiseAdd(a, qh.h);
        MatrixUtil.multiply(aPSD2, 0.5);
        */

        return aPSD;
    }

    /**
     * Given a matrix A that is not necessarily symmetric, find the
     * nearest symmetric matrix to A.
     * <pre>
     * References:
     * https://nhigham.com/2021/01/26/what-is-the-nearest-positive-semidefinite-matrix/
     *
     * </pre>
     @param a a square matrix which can be non-symmetric.
     @return a matrix which is the nearest symmetric to a
     */
    public static double[][] nearestSymmetricToA(double[][] a) {
        if (!isSquare(a)) {
            throw new IllegalArgumentException("a must be square");
        }
        //0.5 * (A + A^T)
        int n = a.length;
        double[][] out = new double[n][];
        int i, j;
        for (i = 0; i < a.length; ++i) {
            out[i] = new double[n];
            for (j = 0; j < a[i].length; ++j) {
                out[i][j] = 0.5 * (a[i][j] + a[j][i]);
            }
        }
        return out;
    }

    /**
     * Given a matrix a that is not necessarily symmetric,
     * and a nonnegative number eps, find the
     * nearest symmetric positive definite matrices with eigenvalues at least eps.
     * Note that this method attempts to make it symmetric positive definite by
     * adding a small perturbation of size smallest eigenvalue to the diagonal
     * of the resulting matrix.  The result satisfies the broader, non-negative
     * definition of positive semi-definite matrix.
     * <pre>
     * References:
     * https://nhigham.com/2021/01/26/what-is-the-nearest-positive-semidefinite-matrix/
     *
     * </pre>
     @param a a square matrix which can be non-symmetric.
     @param eps a tolerance or error above 0 such as machine precision.  must
     * be greater than or equal to 0.  If it is above 0, this method attempts
     * to return a symmetric positive definite matrix.
     @return matrix which is the nearest positive semidefinite to a
     * @throws no.uib.cipr.matrix.NotConvergedException  thrown by MTJ when SVD could not converge
     */
    public static double[][] nearestPositiveSemidefiniteToA(double[][] a,
        double eps) throws NotConvergedException {

        if (eps < 0) {
            throw new IllegalArgumentException("eps must be .gte. 0");
        }
        // uses the Frobenius norm as a distance

        if (!isSquare(a)) {
            System.err.println("matrix must be square for nearestPositiveSemidefiniteToA()");
        }

        // uses the Frobenius norm as a distance.  notation: ||A||_F

        //     the symmetric part of A is matrix B = 0.5*(A + A^T)
        //     the skew-symmetric part of A is matrix C = 0.5*(A - A^T)
        //                          note that C = -C^T.

        // another solution for X uses polar decomposition B = UH,
        //     where U is orthogonal (U*U^T=I)
        //     and H is symmetric positive definite (H = H^T >= 0)
        //         (with H = Q * diag(lambda_i) * Q^T.
        //
        //   since A is symmetric):
        //   answer is X_F = (A + H)/2

        // can use polar decomposition on square matrices:

        double[][] b = nearestSymmetricToA(a);

        /* another to method to solve for X_F = PSD
        // [Q,d] = eig(B,'vector'); X_F = Q*(max(d,0).*Q'); X_F = (X_F + X_F')/2;
        // where Q are the right eigenvectors, and D is a vector of corresponding eigenvalues
        EVD evdB = EVD.factorize(new DenseMatrix(b));
        double[][] q = MatrixUtil.convertToRowMajor(evdB.getRightEigenvectors());
        double[] d = evdB.getRealEigenvalues();
        double dMax = MiscMath0.findMax(d);
        double[][] xF = MatrixUtil.transpose(q);
        MatrixUtil.multiply(xF, dMax);
        xF = MatrixUtil.multiply(q, xF);
        xF = MatrixUtil.pointwiseAdd(xF, MatrixUtil.transpose(xF));
        MatrixUtil.multiply(xF, 0.5);
        EVD evdXF = EVD.factorize(new DenseMatrix(xF));
        double[] eigXF = evdXF.getRealEigenvalues();
        // no perturbations needed for xF
         */

        // H is formed from V, S, and V^T of SVD(B)
        QH qh = performPolarDecomposition(b);
        double[][] h = qh.h;

        //X = (B + H)/2.
        double[][] aPSD = MatrixUtil.pointwiseAdd(b, h);
        MatrixUtil.multiply(aPSD, 0.5);

        // check that all eigenvalues are >= eps
        EVD evd = EVD.factorize(new DenseMatrix(aPSD));
        int i;
        double[] eig = evd.getRealEigenvalues();
        boolean ok = true;
        for (i = 0; i < eig.length; ++i) {
            if (eig[i] < eps) {
                ok = false;
                break;
            }
        }
        if (!ok) {
            //System.out.printf("before: evd eigenvalues=%s\n", FormatArray.toString(evd.getRealEigenvalues(), "%.5e"));

            //still needs a small perturbation to make the eigenvalues all >= eps
            // see https://nhigham.com/2021/02/16/diagonally-perturbing-a-symmetric-matrix-to-make-it-positive-definite/

            double[] e = Arrays.copyOf(eig, eig.length);
            Arrays.sort(e);
            double eigMin = Math.max(-e[0], eps);
            double[][] D = MatrixUtil.zeros(a.length, a.length);
            for (i = 0; i < D.length; ++i) {
                D[i][i] = eigMin;
            }
            aPSD = MatrixUtil.pointwiseAdd(aPSD, D);
            //evd = EVD.factorize(new DenseMatrix(aPSD));
            //System.out.printf("after:  evd eigenvalues=%s\n", FormatArray.toString(evd.getRealEigenvalues(), "%.5e"));
        }

        return aPSD;
    }

    /**
     * normalize vector v by power p, that is the (1/p) power of sum of
     * its (components)^p.  notation is sometimes ||v||_p.
     <pre>
     when p = 0, this is the manhattan normalization or taxi-cab normalization,
     when p = 2, this is the euclidean normalization,
     when p = Double.POSITIVE_INFINITY, this is the max.
     </pre>
     @param v array
     @param p
     @return array v normalized by its lp_p sum.
     */
    public static double[] normalizeLP(double[] v, double p) {
        double divisor = 0;
        if (Double.isInfinite(p)) {
            divisor = MiscMath0.findMax(v);
        } else {

            for (double a : v) {
                divisor += Math.pow(a, p);
            }
            divisor = Math.pow(divisor, 1. / p);
        }
        double[] out = Arrays.copyOf(v, v.length);
        if (divisor != 0.) {
            for (int i = 0; i < v.length; ++i) {
                out[i] /= divisor;
            }
        }
        return out;
    }

    /**
     * the outer product of vectors v1 and v2, which is v1 as a single row matrix
     * and v2 as a single column matrix, so is v1 * v2^T.
     @param v1 array 1
     @param v2 array 2
     @return the outer product of v1 and v2 as double array of
     * size v1.length X v2.length.
     */
    public static double[][] outerProduct(double[] v1, double[] v2) {
        int n = v1.length;
        int m = v2.length;
        double[][] out = new double[n][m];
        for (int row = 0; row < n; ++row) {
            for (int col = 0; col < m; ++col) {
                out[row][col] = v1[row] * v2[col];
            }
        }
        return out;
    }

    /**
     * TODO: write a version for sparse graphs.
     * determine the largest eigenvalue using the power method.  note that
     * matrix A must be diagonalize-able, that is, a positive definite matrix.
     * for best results, perform standard normalization on matrix A first
     * because the first initial guess of an eigenvector of a is composed
     * of random values between [0 and 1).
     * The method is implemented from pseudocode in Golub and van Loan
     * "Matrix Computations".
     *
     * Note that this method normalizes the rows of matrix a during operation.
     * If one is using a matrix that expects normalized columns, such as the transition
     * matrix of page rank in MMDS chap 5, then you'll want to transpose the
     * matrix before using this method).
     *
     * calculates lambda in lambda * v = M * v for some constant eigenvalue lambda.
     *
     * NOTE that the number of necessary iterations is dependent upon
     * how close the largest and second largest eigenvalues are and that ratio
     * tends to be near "1" for large matrices and in that case, the power
     * method isn't the right method (consider QR or SVD).
     * TODO:consider implementing the inverse power method also to determine the
     * smallest eigenvalue and its eigenvector
     @param a a positive definite matrix.  if you have a transition matrix of normalized columns,
     you'll want to transpose the matrix before using this method.
     @param nIterations number of iterations to use
     @return the largest eigenvalue of a
     */
    public static double powerMethod(double[][] a, int nIterations) {

        // is diagonalizable.   From wikipedia: "Real symmetric matrices are diagonalizable by orthogonal matrices; i.e.,
        // given a real symmetric matrix A, Q^T * A * Q is diagonal for some orthogonal matrix Q
        if (!isSquare(a)) {
            throw new IllegalArgumentException("a must be a square matrix");
        }

        // v_k = A^k * v_0  = (c_1*(lambda_1)^k * x_1) + ... (c_n*(lambda_n)^k * x_n)

        int nR = a.length;
        double[] v = new double[nR];
        double[] z;
        double norm;
        double eig = 0;
        int row;

        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        //avoid orthogonal first guess at v using randomization.
        for (row = 0; row < nR; ++row) {
            v[row] = rand.nextDouble();
        }

        for (int i = 0; i < nIterations; ++i) {

            z = MatrixUtil.multiplyMatrixByColumnVector(a, v);
            norm = 0;
            for (row = 0; row < nR; ++row) {
                norm += (z[row]*z[row]);
            }
            norm = Math.sqrt(norm);
            eig = norm;
            for (row = 0; row < nR; ++row) {
                v[row] = z[row] / eig;
            }
            //System.out.printf("eig=%.3f\n  v=%s\n  z=%s\n", eig, Arrays.toString(v),
            //    Arrays.toString(z));
        }
        return eig;
    }

    /**
     * determine the largest eigenvalue using the power method.  note that
     * matrix A must be diagonalizable, that is, a positive definite matrix.
     * for best results, perform standard normalization on matrix A first
     * because the first initial guess of an eigenvector of a is composed
     * of random values between [0 and 1).
     * The method is implemented from pseudocode in Golub and van Loan
     * "Matrix Computations".
     * NOTE that the number of necessary iterations is dependent upon
     * how close the largest and second largest eigenvalues are and that ratio
     * tends to be near "1" for large matrices and in that case, the power
     * method isn't the right method (consider QR or SVD).
     * The eigenvalue is returned and the eigenvector is returned by setting it into the elements of x.
     * A * x = lamba * x where lambda is the eigenvalue corresponding to eigenvector x and x is the
     * largest eigenvalue of matrix A, that is, the principal eigenvector.
     @param a a positive definite matrix
     @param tolerance iterations are stopped when the current multiplication vector
     * difference from previous is smaller than tolerance for each item.
     @param x an initialized vector of size a.length that will be filled by
     * this method to hold the vector used to calculate eig = x^T * a * x
     @return the largest eigenvalue of a
     */
    public static double powerMethod(double[][] a, double tolerance, double[] x) {
        int nR = a.length;
        double[] v = new double[nR];
        double[] z;
        double norm, t;
        double eig = 0;
        int row, stop;

        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        //avoid orthogonal first guess at v using randomization.
        for (row = 0; row < nR; ++row) {
            v[row] = rand.nextDouble();
        }

        int nIter = 0;

        while (true) {

            z = MatrixUtil.multiplyMatrixByColumnVector(a, v);
            norm = 0;
            for (row = 0; row < nR; ++row) {
                norm += (z[row]*z[row]);
            }
            norm = Math.sqrt(norm);
            eig = norm;
            stop = 1;
            for (row = 0; row < nR; ++row) {
                t = z[row]/eig;
                if (Math.abs(v[row] - t) > tolerance) {
                    stop = 0;
                }
                v[row] = t;
            }
            if (stop == 1) {
                break;
            }
            //System.out.printf("nIter=%d eig=%.3f\n  v=%s\n  z=%s\n", nIter,
            //    eig, Arrays.toString(v), Arrays.toString(z));
            nIter++;
        }
        System.arraycopy(v, 0, x, 0, v.length);
        return eig;
    }

     /**
     * determine the largest eigenvalue using the power method.  note that
     * array a must be diagonalizable, that is, a positive definite matrix.
     * for best results, perform standard normalization on matrix a first
     * because the first initial guess of an eigenvector of a is composed
     * of random values between [0 and 1).
     * The method is implemented from pseudocode in Golub and van Loan
     * "Matrix Computations".
     * NOTE that the number of necessary iterations is dependent upon
     * how close the largest and second largest eigenvalues are and that ratio
     * tends to be near "1" for large matrices and in that case, the power
     * method isn't the right method (consider QR or SVD).
     @param a a positive definite matrix
     @param tolerance iterations are stopped when the current multiplication vector
     * difference from previous is smaller than tolerance for each item.
     @return the largest eigenvalue of a
     */
    public static double powerMethod(double[][] a, double tolerance) {
        double[] x = new double[a.length];
        return powerMethod(a, tolerance, x);
    }

    /**
     * determine the eigenvalue pairs using the power method.  note that
     * array a must be diagonalizable, that is, a positive definite matrix.
     * for best results, perform standard normalization on matrix a first
     * because the first initial guess of an eigenvector of a is composed
     * of random values between [0 and 1).
     * The method follows "Mining of Massive Datasets" by Leskovec, Rajaraman,
     * and Ullman.  http://www.mmds.org/
     @param a a positive definite matrix
     @param tolerance iterations are stopped when the current multiplication vector
     * difference from previous is smaller than tolerance for each item.
     @return an array of a.length eigenvectors
     */
    public static double[] powerMethodEigenPairs(double[][] a, double tolerance) {

        double[] x = new double[a.length];
        double eig;
        double[] eigs = new double[a.length];
        double[][] x2, a2;

        for (int nr = 0; nr < a.length; ++nr) {

            eigs[nr] = powerMethod(a, tolerance, x);

            x2 = outerProduct(x, x);
            a2 = copy(a);
            for (int i = 0; i < a2.length; ++i) {
                for (int j = 0; j < a2[i].length; ++j) {
                    a2[i][j] -= (eigs[nr] * x2[i][j]);
                }
            }
            a = a2;

            System.out.printf("eig[%d]=%.5f x=%s\n", nr, eigs[nr], Arrays.toString(x));
            System.out.println("  a2=\n");
            for (int i = 0; i < a2.length; ++i) {
                System.out.printf("  %s\n", Arrays.toString(a2[i]));
            }
        }

        return eigs;
    }

    /**
     * calculate (matrix A)^power using its eigen decompostion:
     * A^power = S * (Delta^power) * S^-1
     * where S holds eigenvectors in its columns.  Delta is a diagonal matrix
     * holding the eigenvalues.
     @param a a square matrix
     @param power integer power to apply to matrix a
     @return  (matrix A)^power
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public double[][] powerOf(double[][] a, int power) throws NotConvergedException {
        if (!isSquare(a)) {
            throw new IllegalArgumentException("matrix a must be square");
        }
        double eps = 1E-15;
        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(a);

        int rank = rank(svd);

        if (rank < svd.u.length) {
            // this case has eigenvalues that are 0.
            //  det(A - lambda*I) = 0, so is not diagonalizable.
            throw new IllegalStateException("there are eigenvalues with value 0, so matrix a is not diagonalizable");
        }

        double[] d = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            d[i] = Math.pow(svd.s[i], power);
        }

        double[][] s = svd.u;

        double[][] sInv;
        // check that the eigenvectors are independent.
        // using left divide notation:
        //  x = A\B solves the system of linear equations A*x = B for x.
        //  X = A\B in MTJ is X = A.solve(B, X), that is, inputs are A and B.
        //  S*S^-1 = I
        //  S.solve(I, sInv) to get sInv using MTJ:
        DenseMatrix _sInv = new DenseMatrix(s.length, s.length);
        _sInv = (DenseMatrix) new DenseMatrix(s).solve(
             new DenseMatrix(MatrixUtil.createIdentityMatrix(s.length)),
             _sInv);
        sInv = MatrixUtil.convertToRowMajor(_sInv);

        /*
            // check this condition... S not invertible so check math in A = S^-1 * Delta * S
            //                         where S^-1 is pseudoinverse
            //sInv = svd(s).V * pseudoinverse(svd(s).s) * svd(s).U^T
            sInv = MatrixUtil.pseudoinverseRankDeficient(s, false);
            double[][] chk = MatrixUtil.multiply(s, sInv);
            System.out.printf("check S*S^-1 ~ I:\n%s\n", FormatArray.toString(chk, "%.5e"));
        */

        // A = S^-1 * delta^power * S
        double[][] aP = MatrixUtil.multiply(MatrixUtil.multiplyByDiagonal(s, d), sInv);

        return aP;
    }

    /**
     * calculate the square root of symmetric positive definite matrix A using SVD.
     *
     * <pre>
     *    [U, S, V] = svd(A)
     *    J = V * S^(1/2) * V^T is a symmetric n×n matrix, such that square root of A = JJ.
     *    J is non-negative definite.
     * </pre>
     * from Allan Jepson's lecture on Gilbert Strang's SVD in machine learning
     * http://www.cs.toronto.edu/~jepson/csc420/notes/introSVD.pdf
     * Also see Chap 7.4 of "Introduction to LinearAlgebra" by Strang, the section
     * on Polar Decomposition.
     @param a a square symmetric positive definite matrix.  If the matrix is not
     *          positive definite matrix, use nearestPositiveSemidefiniteToA() first.
     @return 
     * @throws no.uib.cipr.matrix.NotConvergedException thrown by MTJ when SVD could not converge
     */
    public static double[][] squareRoot(double[][] a) throws NotConvergedException {

        if (!MatrixUtil.isPositiveDefinite(a)) {
            throw new IllegalArgumentException("a must be a non-negative definite matrix");
        }

        int m = a.length;
        int n = a[0].length;

        if (m < n) {
            a = MatrixUtil.transpose(a);
            m = a.length;
            n = a[0].length;
        }

        // limit for a number to be significant above 0
        double eps = 1e-16;

        DenseMatrix aMatrix = new DenseMatrix(a);
        SVD svd = SVD.factorize(aMatrix); // U is mxn; S=nxn; V=nxn

        // s is an array of size min(m,n)  which is nxn here
        double[] s = svd.getS();

        double[][] sMatrix = new double[n][n];
        for (int i = 0; i < n; ++i) {
            sMatrix[i] = new double[n];
            sMatrix[i][i] = Math.sqrt(s[i]);
        }

        int rank = rank(svd);

        // the number of distinct non-zero eigenvalues is rank.
        // there are 2^(rank) square roots

        /*
        S has singular values along the diagonal and those are the square roots
        of positive eigenvalues of A.  the number of them is == the rank of A.
        NOTE the singular values are the eigenvalues of A*A^T of A^T*A.

        S = [ sqrt(λ_1)     0        ...
            [     0      sqrt(λ_2)   ...

        */

        /*
        U is mxm orthonormal columns
        S is mxn with non-negative singular values.  rank is number of non-zero entries
        V is  nxn

            [nxn] * [nxn] * [nxn]
        J = V * S^(1/2) * V^T  is [nxn]
        */

        /*
          a00  a01   s00  0
          a10  a11     0  s11
          a20  a21

        a_r0_c0 * s_r0_c0 + 0(=s_r1_c0)    0(=s_r0_c1) + a_r0_c1 * s_r1_c1
        a_r1_c0 * s_r0_c0 + 0(=s_r1_c0)    0(=s_r0_c1) + a_r1_c1 * s_r1_c1
        a_r2_c0 * s_r0_c0 + 0(=s_r1_c0)    0(=s_r0_c1) + a_r2_c1 * s_r1_c1
        =
        a_r0_c0 * s[0]    a_r0_c1 * s[1]
        a_r1_c0 * s[0]    a_r1_c1 * s[1]
        a_r2_c0 * s[0]    a_r2_c1 * s[1]
        */
        DenseMatrix vT = (DenseMatrix) svd.getVt();
        //DenseMatrix u = (DenseMatrix) svd.getU();

        double[][] _vT = Matrices.getArray(vT);

        double[][] j = MatrixUtil.multiply(MatrixUtil.transpose(_vT), sMatrix);
        j = MatrixUtil.multiply(j, _vT);

        // Strang,"Introduction to Linear Algebra", chap 7, Section Polar Decomposition"
        // A = U * Σ * V^T
        //   = (U * V^T) * (V * Σ * V^T)
        //   =  (Q)      *  (H)

        return j;
    }

    /**
     *
     @param a
     @return
     */
    public static boolean isSquare(double[][] a) {
        return (a.length == a[0].length);
    }

    /**
     *  checks for a_ij == a_ji.  no normalization for row factors is made.
     *  that is, there may be a factor which makes the matrix symmetric.
        e.g. [33, 24]
             [48, 57] here, divide by 2 and the matrix is symmetric
       
     @param a matrix
     @param tol
     @return true if a is symmetric
     */
    public static boolean isSymmetric(double[][] a, double tol) {
        if (!isSquare(a)) {
            return false;
        }
        // compare aij with aji

        int n = a.length;
        int i, j;
        for (i = 0; i < n; ++i) {
            for (j = i; j < n; ++j) {
                if (Math.abs(a[i][j] - a[j][i]) > tol) {
                    return false;
                }
            }
        }

        // if false, there may be a factor which makes the matrix symmetric.
        // e.g. [33, 24]
        //      [48, 57] <--- divide by 2 and the matrix is symmetric
        //

        return true;
    }

    /**
     *
     @param a
     @return
     */
    public static boolean isPositiveSymmetric(double[][] a) {
        if (!isSquare(a)) {
            return false;
        }

        double tol = 1e-3;
        int n = a.length;
        int i, j;
        for (i = 0; i < n; ++i) {
            for (j = i; j < n; ++j) {
                if (Math.abs(a[i][j]) < 1e-17 || Math.abs(a[j][i]) < 1e-17) {
                    return false;
                }
                if (Math.abs(a[i][j] - a[j][i]) > tol) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * holds data structures for having solved for the vector x which is the
     * closest to a given vector b in the subspace defined by A which is
     * n columns of linearly independent vectors of length m (they are in
     * real space R^m).
     * x is an approximation so is noted as x^{hat}.
     * the projection p = A*x^{hat}.
     * The matrix projection P = p*b.
     */
    public static class ProjectionResults {
        double[] x;
        double[] p;
        double[][] pMatrix;
    }
    /**
     * solve for the vector x which is the closest to a given vector b in the
     * subspace defined by A which is n columns of linearly independent vectors
     * of length m (they are in real space R^m).
     * x is an approximation so is noted as x^{hat}.
     * the projection p = A*x^{hat}.
     * The matrix projection P = p*b.
     <pre>
     Chap 4 of the book "Introduction to Linear Algebra" by W Gilbert Strang
     </pre>
     @param a subspace defined by A which is n columns of linearly
     * independent vectors of length m (they are in real space R^m).
     @param b  array b
     @return matrix projection P = p*b where p = A*x^{hat}
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static ProjectionResults projection(double[][] a, double[] b) throws NotConvergedException {

        /*
        from Chap 4 of the book "Introduction to Linear Algebra" by W Gilbert Strang,
           x^{hat} = (A^T*A)^-1 * A^T * b
         and p = A * x^{hat}
         and P = p * b
        */

        ProjectionResults pr = new ProjectionResults();

        // A_pseudoinverse = inverse(A^T*A) * A^T
        // [a[0].length][a.length]
        double[][] aPseudoInv = MatrixUtil.pseudoinverseFullColumnRank(a);

        pr.x = MatrixUtil.multiplyMatrixByColumnVector(aPseudoInv, b);

        pr.p = MatrixUtil.multiplyMatrixByColumnVector(a, pr.x);

        pr.pMatrix = MatrixUtil.multiply(a, aPseudoInv);

        return pr;
    }

    /**
     * check whether all eigenvalues are >=  tol
     * @param a
     * @param tol
     * @return
     */
    public static boolean allEigenValuesGEQ(double[][] a, double tol) throws NotConvergedException {
        EVD evd = EVD.factorize(new DenseMatrix(a));
        int i;
        double[] eig = evd.getRealEigenvalues();
        for (i = 0; i < eig.length; ++i) {
            if (eig[i] < tol) {
                return false;
            }
        }
        return true;
    }

    /**
     * A matrix is positive definite if it’s symmetric and all its eigenvalues are positive,
     * which means its pivots are positive.
     * a quick test is that x^T * A * x is strictly positive for every non-zero
     * column vector x of n real numbers where A is n x n.
     * Every positive definite matrix can be factored into L*D*L^T with positive pivots.
     * Also note that if a has n independent columns (r==n), A^T*A is positive definite.
     @param a matrix a
     @return true if a is positive definite
     */
    public static boolean isPositiveDefinite(double[][] a) throws NotConvergedException {

        if (true) {
            // checking that eigenvalues are all >0 within machine precision
            return allEigenValuesGEQ(a, 1E-11);
        }
        if (true) {
            DenseCholesky chol = no.uib.cipr.matrix.DenseCholesky.factorize(new DenseMatrix(a));
            return (chol.getU() != null);
        }

        if (!isSymmetric(a, 1E-7)) {
            System.err.println("matrix is not symmetric in current form");
            //return false;
        }

        /*
        from Strang "Introduction to Linear Algebra", chapter 6.5.

        matrix A is positive definite for every non-zero vector x if x^T*A*x > 0

        when a symmetric nxn matrix has 1 of these 4, it has all 4:
            1) both eigenvalues are positive
               Note: mathworks recommends a tolerance for error, so all eigenvalues > tolerance above 0.
            2) all upper left determinants (the 1x1 and 2x2 ... ) are positive
            3) the pivots are positive a>0 and a*c-b^2>0
            4) the function x^T * A * x is positive except at x = 0
        */

        int n = a.length;
        int i;

        double[] x = new double[n];
        for (i = 0; i < n; ++i) {
            x[i] = i + 1;
        }
        double p = MatrixUtil.innerProduct(MatrixUtil.multiplyRowVectorByMatrix(x, a), x);
        //return p > 0;

        // using rule that left upward rooted determinants > 0
        double[][] sa;
        double det;
        for (i = 0; i < n; ++i) {
            sa = MatrixUtil.copySubMatrix(a, 0, i, 0, i);
            det = MatrixUtil.determinant(sa);
            if (det < 0) {
                return false;
            }
        }
        return true;
    }

        /*
         * e.g.    | 1  -5  2 |         | 3 4 |         | 7 4 |         | 7 3 |
         *         | 7   3  4 |  =  1 * | 1 5 |  +  5 * | 2 5 |  +  2 * | 2 1 |  = 11 + 135 + 2 =
 148
         *         | 2   1  5 |
         *
         *          3 4     7  4    7  3
         *          1 5     2  5    2  1
         *
         *         -5 2     1  2    1 -5
         *          1 5     2  5    2  1
         *
         *         -5 2     1  2    1 -5
         *          3 4     7  4    7  3
    */

    /**
     *
     @param m
     @return
     */

    public static double[][] createCofactor(double[][] m) {

        int ncols = m.length;
        int nrows = m[0].length;

        double[][] cofactor = new double[ncols][nrows];

        for (int i = 0; i < ncols; i++) {

            cofactor[i] = new double[nrows];

            boolean si = ((i & 1) == 1); // sign is -

            for (int j = 0; j < nrows; j++) {

                boolean sj = ((j & 1) == 1); // sign is -

                double[][] n = copyExcept(m, i, j);

                double cfctr = determinant(n);

                if (si ^ sj) { // XOR if either is 1 but not both
                    cfctr = -1*cfctr;
                }

                cofactor[i][j] = cfctr;
            }
         }
        return cofactor;
    }

    /**
     * find the equation for which A * A^(-1) = the identity matrix using cramer's rule.
     *
     * note that for a to be invertible, none of its eigenvalues can be 0.
     * also note that if the number of linearly independent vectors os matrix
     * A is equal to the number of columns of A, one can use the spectral
     * decomposition: A^-1 = Q * (delta)^-1 * Q^-1
     * where Q is a matrix whose columns hold eigenvectors and delta is a diagonal
     * matrix holding the eigenvalues.
     *
     *             1
     * A^(-1) =  ------ C^(T)  where C_ij = cofactor of a_ij
     *            det A
     *
     @param m a square invertible matrix.
     @return inverse of matrix m
     */
    public static double[][] inverse(double[][] m) {

        // create cofactor of matrix:
        double[][] cofactor = createCofactor(m);

        double[][] cofactorTransposed = transpose(cofactor);

        double det = determinant(m);

        multiply(cofactorTransposed, 1./det);

        return cofactorTransposed;
    }

    /**
     * given a as vectors of data of nSamples of nVariables, return the
     * mean of each of the variables.
     * note that the format must be a[nSamples][nVariables],
     * e.g. a[0] = [10, 100, 1000]', a[1] = [9, 101, 999]; for nSamples = 2
     * and nVariables = 3;
     @param a matrix
     @return mean of each column as an array of size a[0].length
     */
    public static double[] columnMeans(double[][] a) {
        int nSamples = a.length;
        int nVariables = a[0].length;

        int i, j;
        double[] mean = new double[nVariables];
        double sum;
        for (j = 0; j < nVariables; ++j) {
            sum = 0;
            for (i = 0; i < nSamples; ++i) {
                sum += (a[i][j]);
            }
            mean[j] = sum/(double)nSamples;
        }
        return mean;
    }

    /**
     * given a, return the mean of each row of a.
     @param a matrix
     @return mean of each row as an array of size a.length
     */
    public static double[] rowMeans(double[][] a) {
        int nSamples = a.length;
        int nVariables = a[0].length;

        int i, j;
        double[] mean = new double[nSamples];
        double sum;
        for (i = 0; i < nSamples; ++i) {
            sum = 0;
            for (j = 0; j < nVariables; ++j) {
                sum += (a[i][j]);
            }
            mean[i] = sum/(double)nVariables;
        }
        return mean;
    }


    /**
     * given a as vectors of data of nSamples of nVariables, return the
     * mean of each of the variables.
     * note that the format must be a[nSamples][nVariables],
     * e.g. a[0] = [10, 100, 1000]', a[1] = [9, 101, 999]; for nSamples = 2
     * and nVariables = 3;
     @param a vectors of data of nSamples of nVariables in format [nSamples X nVariables]
     @return standard deviation of the variables.  length is a[0].length
     */
    public static double[] standardDeviation(double[][] a) {
        int nSamples = a.length;
        int nVars = a[0].length;

        double[] c = columnMeans(a);

        double[] out = new double[nVars];
        int i, d;
        double diff;
        for (i = 0; i < nSamples; ++i) {
            for (d = 0; d < nVars; ++d) {
                diff = a[i][d] - c[d];
                out[d] += (diff*diff);
            }
        }
        for (d = 0; d < nVars; ++d) {
            out[d] = Math.sqrt(out[d]/(nSamples - 1.0));
        }

        return out;
    }

    /**
     * calculate a - v*I where A is square matrix and v is a vector.  I is the identity
     * matrix.
     @param a a square matrix.
     @param v a vector of length of a.length.
     @return the matrix a - v*I.
     */
    public static double[][] aMinusVectorTimesIdentity(double[][] a, double[] v) {
        int n = a.length;
        if (n != a[0].length) {
            throw new IllegalArgumentException("a must be a square matrix");
        }
        if (n != v.length) {
            throw new IllegalArgumentException("v must be same length as a");
        }

        double[][] out = copy(a);

        int i;
        for (i = 0; i < n; ++i) {
            out[i][i] -= v[i];
        }

        return out;
    }

    /**
     * calculate the sum of the diagonal elements of a.
     * Note that the trace of matrix A equals the sum of its eigenvalues.
     * Note: the trace of A is equal to the sum of its eigenvalues.
     @param a a square matrix.
     @return the sum of the diagonal elements of a
     */
    public static double trace(double[][] a) {
        int n = a.length;
        if (n != a[0].length) {
            throw new IllegalArgumentException("a must be a square matrix");
        }

        double sum = 0;

        int i;
        for (i = 0; i < n; ++i) {
            sum += a[i][i];
        }

        return sum;
    }


    /**
     * calculate the sum of the diagonal elements of a.
     * Note that the trace of matrix A equals the sum of its eigenvalues.
     * Note: the trace of A is equal to the sum of its eigenvalues.
     @param a a square matrix.
     @return the sum of the diagonal elements of a
     */
    public static float trace(float[][] a) {
        int n = a.length;
        if (n != a[0].length) {
            throw new IllegalArgumentException("a must be a square matrix");
        }

        float sum = 0f;

        int i;
        for (i = 0; i < n; ++i) {
            sum += a[i][i];
        }

        return sum;
    }
    /**
     * calculate the sum of the diagonal elements of v*I (i.e. sum of all elements of v)
     @param v a vector to be treated as diagonal elements of an identity matrix.
     @return the sum of the elements of v
     */
    public static double trace(double[] v) {
        int n = v.length;

        double sum = 0;

        int i;
        for (i = 0; i < n; ++i) {
            sum += v[i];
        }

        return sum;
    }

    /**
     * pointwise multiplication (a.k.a. Hadamard product, entrywise product, and Schur product)
     @param a matrix a
     @param b matrix b
     @return element-wise multiplication of elements in a and b
     */
    public static double[][] pointwiseMultiplication(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;

        if (b.length != m || b[0].length != n) {
            throw new IllegalArgumentException("a and b must have same dimensions");
        }

        int j;
        double[][] out = new double[m][n];
        for (int i = 0; i < m; ++i) {
            out[i] = new double[n];
            for (j = 0; j < n; ++j) {
                out[i][j] = a[i][j] * b[i][j];
            }
        }
        return out;
    }

    /**
     * pointwise multiplication (a.k.a. Hadamard product, entrywise product, and Schur product)
     @param a matrix a
     @param b matrix b
     @return element-wise multiplication of elements in a and b
     */
    public static int[][] pointwiseMultiplication(int[][] a, int[][] b) {
        int m = a.length;
        int n = a[0].length;

        if (b.length != m || b[0].length != n) {
            throw new IllegalArgumentException("a and b must have same dimensions");
        }

        int j;
        int[][] out = new int[m][n];
        for (int i = 0; i < m; ++i) {
            out[i] = new int[n];
            for (j = 0; j < n; ++j) {
                out[i][j] = a[i][j] * b[i][j];
            }
        }
        return out;
    }

    /**
     * pointwise addition
     @param a matrix a
     @param b matrix b
     @return element-wise addition of elements in a and b
     */
    public static double[][] pointwiseAdd(double[][] a, double[][] b) {
        double[][] out = MatrixUtil.zeros(a.length, a[0].length);
        pointwiseAdd(a, b, out);
        return out;
    }

    /**
     * pointwise addition
     @param a matrix a
     @param b matrix b
     @param out the results of element wise add of a + b. Note that it
     * is safe to provide out as the same object as input argument a or b.
     */
    public static void pointwiseAdd(double[][] a, double[][] b, double[][] out) {
        int m = a.length;
        int n = a[0].length;

        if (b.length != m || b[0].length != n) {
            throw new IllegalArgumentException("a and b must have same dimensions");
        }
        if (out.length != m || out[0].length != n) {
            throw new IllegalArgumentException("a and out must have same dimensions");
        }

        int i, j;
        for (i = 0; i < out.length; ++i) {
            for (j = 0; j < out[i].length; ++j) {
                out[i][j] = a[i][j] + b[i][j];
            }
        }
    }

    /** pointwise subtraction
     @param a matrix a
     @param b matrix b
     @return element-wise subtraction of elements in a and b
     */
    public static double[][] pointwiseSubtract(double[][] a, double[][] b) {
        double[][] out = MatrixUtil.zeros(a.length, a[0].length);
        pointwiseSubtract(a, b, out);
        return out;
    }

    /**
     * pointwise subtraction
     @param a matrix
     @param b matrix
     @param out the results of element wise subtraction, a - b. Note that it
     */
    public static void pointwiseSubtract(double[][] a, double[][] b, double[][] out) {
        int m = a.length;
        int n = a[0].length;

        if (b.length != m || b[0].length != n) {
            throw new IllegalArgumentException("a and b must have same dimensions");
        }
        if (out.length != m || out[0].length != n) {
            throw new IllegalArgumentException("a and out must have same dimensions");
        }

        int i, j;
        for (i = 0; i < out.length; ++i) {
            for (j = 0; j < out[i].length; ++j) {
                out[i][j] = a[i][j] - b[i][j];
            }
        }
    }

    /**
     * point-wise subtraction
     @param a matrix a
     @param b matrix b
     @return the results of point-wise subtraction, a - b.
     */
    public static int[][] pointwiseSubtract(int[][] a, int[][] b) {
        int m = a.length;
        int n = a[0].length;

        if (b.length != m || b[0].length != n) {
            throw new IllegalArgumentException("a and b must have same dimensions");
        }

        int[][] out = new int[m][];

        int i, j;
        for (i = 0; i < out.length; ++i) {
            out[i] = new int[n];
            for (j = 0; j < out[i].length; ++j) {
                out[i][j] = a[i][j] - b[i][j];
            }
        }
        return out;
    }

    /**
     * pointwise subtraction
     @param a array
     @param b array
     @param out the results of element wise subtraction, a - b. Note that it
     * is safe to provide out as the same object as input argument a or b.
     */
    public static void pointwiseSubtract(double[] a, double[] b, double[] out) {
        int m = a.length;

        if (b.length != m) {
            throw new IllegalArgumentException("a and b must have same dimensions");
        }
        if (out.length != m) {
            throw new IllegalArgumentException("a and out must have same dimensions");
        }

        int i;
        for (i = 0; i < out.length; ++i) {
            out[i] = a[i] - b[i];
        }
    }

    /**
     * pointwise multiplication (a.k.a. Hadamard product, entrywise product, and Schur product)
     @param a matrix
     @param b matrix
     @return  element-wise multiplication of elements in a and b
     */
    public static double[] pointwiseMultiplication(double[] a, double[] b) {
        int m = a.length;

        if (b.length != m) {
            throw new IllegalArgumentException("a and b must have same dimensions");
        }

        int j;
        double[] out = new double[m];
        for (int i = 0; i < m; ++i) {
            out[i] = a[i] * b[i];
        }
        return out;
    }

    /**
     * dot product, summation_over_i(a[i]*b[i])
     @param a array
     @param b array
     @return dot product of a and b which is a.k.a.inner product
     */
    public static double dot(double[] a, double[] b) {
        int m = a.length;

        if (b.length != m) {
            throw new IllegalArgumentException("a and b must have same dimensions");
        }

        double out = 0;
        for (int i = 0; i < m; ++i) {
            out += a[i] * b[i];
        }
        return out;
    }

    /**
     * create an array of zeros
     @param nRows number of rows for new matrix
     @param nCols number of columns for new matrix
     @return a new empty double array of size [nRows X nCols] in row major notation.
     */
    public static double[][] zeros(int nRows, int nCols) {
        double[][] out = new double[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            out[i] = new double[nCols];
            // java, by default, initializes with zeroes
        }
        return out;
    }

    /**
     *
     @param nRows
     @return
     */
    public static double[][] createIdentityMatrix(int nRows) {
        double[][] out = new double[nRows][nRows];
        for (int i = 0; i < nRows; ++i) {
            out[i] = new double[nRows];
            out[i][i] = 1;
        }
        return out;
    }

    /**
     * right divide is pointwise division
     @param a matrix a
     @param b matrix b
     @return element-wise division of elements in a and b
     */
    public static double[][] pointwiseDivision(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;

        if (b.length != m || b[0].length != n) {
            throw new IllegalArgumentException("a and b must have same dimensions");
        }

        int j;
        double[][] out = new double[m][n];
        for (int i = 0; i < m; ++i) {
            out[i] = new double[n];
            for (j = 0; j < n; ++j) {
                out[i][j] = a[i][j] / b[i][j];
            }
        }
        return out;
    }

    /**
     * right divide is pointwise division, that is a[i]/b[i] for i = [0, a.length).
     @param a array a
     @param b array b
     @return element-wise division of elements in a and b
     */
    public static double[] pointwiseDivision(double[] a, double[] b) {
        int m = a.length;

        if (b.length != m) {
            throw new IllegalArgumentException("a and b must have same lengths");
        }

        int j;
        double[] out = Arrays.copyOf(a, m);
        for (int i = 0; i < m; ++i) {
            out[i] /= b[i];
        }
        return out;
    }

    /**
     * perform a left-right swap of the columns of a, flipping the matrix
     * vertically.  the method mimics matlab's flipur.
     @param a  matrix a
     */
    public static void flipLR(double[][] a) {
        int idxHi = a[0].length - 1;
        int idxLo = 0;
        int n = idxHi - idxLo + 1;
        int end = idxLo + (n/2);

        double swap;
        int col, row, idx2;
        int count = 0;
        for (col = idxLo; col < end; col++) {
            idx2 = idxHi - count;
            for (row = 0; row < a.length; row++) {
                swap = a[row][col];
                a[row][col] = a[row][idx2];
                a[row][idx2] = swap;
            }
            count++;
        }
    }

    /**
     * perform an up-down swap of the rows of a, flipping the matrix
     * horizontally.  the method mimics matlab's flipud.
     @param a matrix a
     */
    public static void flipUD(double[][] a) {
        int idxHi = a.length - 1;
        int idxLo = 0;
        int n = idxHi - idxLo + 1;

        int end = idxLo + (n/2);
        double[] swap;
        int count = 0;
        int i, idx2;
        for (i = idxLo; i < end; i++) {
            idx2 = idxHi - count;
            swap = a[i];
            a[i] = a[idx2];
            a[idx2] = swap;
            count++;
        }
    }

    /**
     * solves for vector x in the equation L*x=b where L is the lower triangular
     * matrix and b is a vector.
     * runtime complexity is approx (b.length)^2.
     @param lowerTriangular the lower triangular matrix
     @param b vector on the righthand side of the equation L*x=b
     @return x in equation L*x = b
     */
    public static double[] forwardSubstitution(double[][] lowerTriangular, double[] b) {
        double[] outX = new double[b.length];
        forwardSubstitution(lowerTriangular, b, outX);
        return outX;
    }

    /**
     * solves for vector x in the equation L*x=b where L is the lower triangular
     * matrix and b is a vector.
     * runtime complexity is approx (b.length)^2.
     * method follows Golub and Van Loan algorithm 4.1-1.
     @param lowerTriangular the lower triangular matrix
     @param b vector on the righthand side of the equation L*x=b
     @param outX output variable x in equation L*x = b.  length is b.length.
     */
    public static void forwardSubstitution(double[][] lowerTriangular, double[] b,
        double[] outX) {

        int m = b.length;

        if (lowerTriangular[0].length != m) {
            throw new IllegalArgumentException("the number of columns in "
                    + "lowerTriangular must equal the length of b");
        }
        if (outX.length != m) {
            throw new IllegalArgumentException("outX.length must equal the length of b");
        }

        int i, j;
        for (i = 0; i < m; ++i) {
            outX[i] = b[i];
            for (j = 0; j <= (i-1); ++j) {
                outX[i] -= (lowerTriangular[i][j] * outX[j]);
            }
            outX[i] /= lowerTriangular[i][i];
        }
    }

    /**
     * solves for vector x in the equation L*x=b where L is the lower triangular
     * matrix and b is a vector.
     * runtime complexity is approx (b.length)^2.
     * method follows Golub and Van Loan algorithm 4.1-1.
     @param lowerTriangular the lower triangular matrix
     @param b vector on the righthand side of the equation L*x=b
     @return x in equation L*x = b
     */
    public static double[] forwardSubstitution(LowerTriangDenseMatrix lowerTriangular, double[] b) {

        int m = b.length;

        if (lowerTriangular.numColumns() != m) {
            throw new IllegalArgumentException("the number of columns in "
                    + "lowerTriangular must equal the length of b");
        }

        double[] x = new double[b.length];

        int i, j;
        for (i = 0; i < m; ++i) {
            x[i] = b[i];
            for (j = 0; j <= (i-1); ++j) {
                x[i] -= (lowerTriangular.get(i, j) * x[j]);
            }
            x[i] /= lowerTriangular.get(i, i);
        }

        return x;
    }

    /**
     * solves for vector x in the equation U*x = y where
     * U is an upper triangular matrix and y is a vector.
     * runtime complexity is approx (y.length)^2.
     * method follows Golub and Van Loan algorithm 4.1-2.
     @param upperTriangular the upper triangular matrix
     * (a_i_j=0 where i .gt. j)
     * <pre>
     *     0  1  2
        2  *  *  *
        1  *  *
        0  *
           0  1  2
     * </pre>
     @param y vector on righthand side of equation
     @return x in equation U*x = y
     */
    public static double[] backwardSubstitution(double[][] upperTriangular, double[] y) {

        double[] outX = new double[y.length];

        backwardSubstitution(upperTriangular, y, outX);

        return outX;
    }

    /**
     * solves for vector x in the equation U*x = y where
     * U is an upper triangular matrix and y is a vector.
     * runtime complexity is approx (y.length)^2.
     @param upperTriangular the upper triangular matrix
     * (a_i_j=0 where i .gt. j)
     * <pre>
     *     0  1  2
        2  *  *  *
        1  *  *
        0  *
           0  1  2
     * </pre>
     @param y vector on righthand side of equation
     @param outX output variable x in equation U*x = y.  must be length y.length.
     */
    public static void backwardSubstitution(double[][] upperTriangular, double[] y,
        double[] outX) {

        int m = y.length;

        if (upperTriangular[0].length != m) {
            throw new IllegalArgumentException("the number of columns in "
                    + "upperTriangular must equal the length of y");
        }
        if (outX.length != m) {
            throw new IllegalArgumentException("outX.length must equal the length of y");
        }

        int i, j;
        for (i = m-1; i >= 0; i--) {
            outX[i] = y[i];
            for (j = i+1; j < m; ++j) {
                outX[i] -= (outX[j]*upperTriangular[i][j]);
            }
            outX[i] /= upperTriangular[i][i];
        }
    }

    /**
     * solves for vector x in the equation U*x = y where
     * U is an upper triangular matrix and y is a vector.
     * runtime complexity is approx (y.length)^2.
     @param upperTriangular the upper triangular matrix
     @param y vector on righthand side of equation
     @return x in equation U*x = y
     */
    public static double[] backwardSubstitution(UpperTriangDenseMatrix upperTriangular,
        double[] y) {

        int m = y.length;

        if (upperTriangular.numColumns() != m) {
            throw new IllegalArgumentException("the number of columns in "
                    + "upperTriangular must equal the length of y");
        }

        double[] x = new double[m];
        int i, j;
        for (i = m-1; i >= 0; i--) {
            x[i] = y[i];
            for (j = i+1; j < m; ++j) {
                x[i] -= (x[j]*upperTriangular.get(i, j));
            }
            x[i] /= upperTriangular.get(i, i);
        }
        return x;
    }

    /**
     *
     @param a
     @param value
     */
    public static void fill(double[][] a, double value) {
        int i;
        for (i = 0; i < a.length; ++i) {
            Arrays.fill(a[i], value);
        }
    }

    /**
     *
     @param a
     @return
     */
    public static BlockMatrixIsometric transpose(BlockMatrixIsometric a) {

        double[][] _b = MatrixUtil.zeros(a.getA()[0].length, a.getA().length);

        BlockMatrixIsometric b = new BlockMatrixIsometric(_b, a.getBlockSize1(), a.getBlockSize0());

        transpose(a.getA(), b.getA());

        /*
        double[][] block = MatrixUtil.zeros(a.getBlockSize0(), a.getBlockSize1());
        double[][] blockT = MatrixUtil.zeros(a.getBlockSize1(), a.getBlockSize0());

        int nb0 = a.getA().length / a.getBlockSize0();
        int nb1 = a.getA()[0].length / a.getBlockSize1();

        int i, j;

        for (i = 0; i < nb0; ++i) {
            for (j = 0; j < nb1; ++j) {
                a.getBlock(block, i, j);
                transpose(block, blockT);
                b.setBlock(blockT, j, i);
                System.out.printf("(%d,%d) blockT=\n%s\n", j, i, FormatArray.toString(blockT, "%.3f"));
                System.out.printf("b=\n%s\n", FormatArray.toString(b.getA(), "%.3f"));
            }
        }
        */
        
        return b;
    }
    
    /**
     * rewrite matrix a into a vector using the order of all rows of column 0,
     * then all rows of column 1, etc.
     @param a matrix a
     @return matrix a reshaped into an array
     */
    public static double[] reshapeToVector(double[][] a) {
        double[] out = new double[a.length * a[0].length];
        int c = 0;
        int i, j;
        for (j = 0; j < a[0].length; ++j) {
            for (i = 0; i < a.length; ++i) {
                out[c] = a[i][j];
                c++;
            }
        }
        return out;
    }
     
    /**
     * calculate a pre-conditoner matrix based upon the matrix U in
     * the decompositon P*A = L*U where P is a permutation matrix,
     * L is a lower triangular matrix, U is an upper triangular matrix.
     * 
     * see Bjork 1991 Section 4.3, "Algorithms for linear least squares problems",
     * especially the end of the section.
     * 
     * from wikipedia:
     * LU factorization with partial pivoting: It turns out that a proper 
     * permutation in rows (or columns) is sufficient for LU factorization. 
     * LU factorization with partial pivoting (LUP) refers often to LU 
     * factorization with row permutations only.
     * 
     * The choice of a good preconditioner may improve the speed of an iterative
     * method.
     * 
     * Note that if matrix a has full rank, the ideal choice of a preconditioner
     * is instead of this method's results, is to use the Cholesky factor of 
     * A^T*A. 
     * Note that another preconditioner uses the columns of a to form a diagonal
     * matrix:   D^(1/2) = diag( sqrt(d_1), ... sqrt(d_n)) where 
     *   d_j = (||a_j||_2)^2  where a_j is column j of a.
     * <pre>
     * using a preconditioner involves reforming A*x = b into
     * A * (M^-1) * y = b, solve for y
     * and M * x = y, solve for x.
     * 
     * Or the left precondition:
     *   (M^-1) * (A*x - b) = 0
     * 
     * </pre>
     * 
     @param a and m X n matrix
     @return matrix with dimensions of a^T
     */
    public static double[][] calculatePreconditionerFromLUP(double[][] a) {
        
        LUP lup = LinearEquations.LUPDecomposition(a);
        
        return lup.u;
    }
    
    /**
     * calculate a pre-conditoner matrix based upon the columns of matrix a to 
     * form a diagonal matrix:   D^(1/2) = diag( sqrt(d_1), ... sqrt(d_n)) where 
     *   d_j = (||a_j||_2)^2  where a_j is column j of a.
     * 
     * see Bjork 1991 Section 4.3, "Algorithms for linear least squares problems",
     * especially the end of the section.
     * 
     * from wikipedia:
     * LU factorization with partial pivoting: It turns out that a proper 
     * permutation in rows (or columns) is sufficient for LU factorization. 
     * LU factorization with partial pivoting (LUP) refers often to LU 
     * factorization with row permutations only.
     * 
     * The choice of a good preconditioner may improve the speed of an iterative
     * method.
     * 
     * Note that if matrix a has full rank, the ideal choice of a preconditioner
     * is instead of this method's results, is to use the Cholesky factor of 
     * A^T*A. 
     * 
     * <pre>
     * using a preconditioner involves reforming A*x = b into
     * A * (M^-1) * y = b, solve for y
     * and M * x = y, solve for x.
     * 
     * Or the left precondition:
     *   (M^-1) * (A*x - b) = 0
     * 
     * </pre>
     * 
     @param a and m X n matrix
     @return matrix with dimensions of a^T
     */
    public static double[][] calculatePreconditionerFromColumns(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        /*
        D^(1/2) = diag( sqrt(d_1), ... sqrt(d_n)) where 
        d_j = (||a_j||_2)^2  where a_j is column j of a.
        */
        int j, i;
        double[][] diag = MatrixUtil.zeros(n, n);
        
        for (j = 0; j < n; ++j) {
            for (i = 0; i < m; ++i) {
                diag[j][j] += (a[i][j] + a[i][j]);
            }
            diag[j][j] = Math.sqrt(diag[j][j]);
        }
        
        return diag;
    }
    
    /**
     * perform a polar decomposition on square matrix a.
     * A = Q*H where Q is orthogonal and H is a symmetric positive semidefinite matrix.  
     If A is invertible, then H is symmetric positive definite.
     The method follows Strang "Introduction to Linear Algebra" Chapter 7 section G.
     @param a square matrix
     @return polar decomposition of a
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static QH performPolarDecomposition(double[][] a) throws NotConvergedException {
        
        if (!isSquare(a)) {
            throw new IllegalArgumentException("a must be a square matrix");
        }
        
        SVDProducts svd = performSVD(a);
        
        double[][] q = MatrixUtil.multiply(svd.u, svd.vT);
        
        double[][] v = MatrixUtil.transpose(svd.vT);
        
        double[][] h = MatrixUtil.multiplyByDiagonal(v, svd.s);
        h = MatrixUtil.multiply(h, svd.vT);
        
        QH qh = new QH();
        qh.q = q;
        qh.h = h;
        
        return qh;
    }
    
    /**
     * a class to hold the results of a polar decomposition
     */
    public static class QH {
        
        /**
         * Q is an orthogonal matrix
         */
        public double[][] q;
        
        /**
         * H is a symmetric positive semidefinite matrix.  
         * If the decomposed matrix A is invertible, then H is
         * symmetric positive definite.
         */
        public double[][] h;
    }
    
    /**
     * given 2 non-decreasing ordered sequences of numbers, find their intersection.
     * The method is called multiset because the sequences may contain more than
     * one element having the same value... the method is used for multisets as multi-sequences.
     * The runtime complexity is O(max(a.length, b.length)).
     @param orderedA an increasing sequence of numbers (i.e. ascending sorted). 
     @param orderedB an increasing sequence of numbers (i.e. sorted by non-decreasing order)
     @return the intersection of orderedA and orderedB
     */
    public static int[] multisetIntersection(int[] orderedA, int[] orderedB) {
        
        TIntList c = new TIntArrayList();
        
        int i = 0, j = 0;
        while (i < orderedA.length && j < orderedB.length) {
            if (orderedA[i] == orderedB[j]) {
                c.add(orderedA[i]);
                i++; 
                j++;
            } else if (orderedA[i] < orderedB[j]) {
                i++;
            } else {
                j++;
            }
        }
        return c.toArray();
    }
    
    /**
     * given 2 sequences of numbers, find their intersection.
     * The method is called multiset because the sequences may contain more than
     * one element having the same value... the method is used for multisets as multi-sequences.
     * The runtime complexity is O(N*log_2(N)) where N is max(a.length, b.length).
     @param a a sequence of numbers 
     @param b a sequence of numbers
     @return intersection of a and b
     */
    public static int[] multisetUnorderedIntersection(int[] a, int[] b) {
        
        a = Arrays.copyOf(a, a.length);
        b = Arrays.copyOf(b, b.length);
        Arrays.sort(a);
        Arrays.sort(b);
        
        return multisetIntersection(a, b);
    }
    
    /**
     *
     @param a
     @return
     */
    public static TIntObjectMap<TIntSet> copy(TIntObjectMap<TIntSet> a) {
        TIntObjectMap<TIntSet> c = new TIntObjectHashMap<TIntSet>();
        TIntObjectIterator<TIntSet> iter = a.iterator();
        TIntSet set, cSet;
        TIntIterator iter2;
        int i, key;
        for (i = 0; i < a.size(); ++i) {
            iter.advance();
            
            key = iter.key();
            set = iter.value();
            
            cSet = new TIntHashSet();
            
            iter2 = set.iterator();
            
            while (iter2.hasNext()) {
                cSet.add(iter2.next());
            }
            c.put(key, cSet);
        }
        return c;
    }
    
    /**
     *
     @param a
     @return
     */
    public static TObjectDoubleMap<PairInt> copy(TObjectDoubleMap<PairInt> a) {
        
        TObjectDoubleMap<PairInt> out = new TObjectDoubleHashMap<PairInt>();
        
        TObjectDoubleIterator<PairInt> iter = a.iterator();
        PairInt p;
        int i;
        for (i = 0; i < a.size(); ++i) {
            iter.advance();
            p = iter.key();
            out.put(p.copy(), iter.value());
        }
        
        return out;
    }
    
    /**
     * create a symmetric adjacency map from a
     @param a an adjacency map
     @return a sthe ymmetric adjacency map of a
     */
    public static TIntObjectMap<TIntSet> copyToSymmetricMap(TIntObjectMap<TIntSet> a) {
        TIntObjectMap<TIntSet> out = new TIntObjectHashMap<TIntSet>();
        
        TIntObjectIterator<TIntSet> iter = a.iterator();
        TIntSet set, uSet, vSet;
        TIntIterator iter2;
        
        int i, u, v;
        for (i = 0; i < a.size(); ++i) {
            iter.advance();
            
            u = iter.key();
            set = iter.value();
            
            uSet = out.get(u);
            if (uSet == null) {
                uSet = new TIntHashSet();
                out.put(u, uSet);
            }
            
            iter2 = set.iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                
                uSet.add(v);
                
                vSet = out.get(v);
                if (vSet == null) {
                    vSet = new TIntHashSet();
                    out.put(v, vSet);
                }
                vSet.add(u);
            }
        }
        
        return out;
    }
    
    /**
     *
     @param a
     @return
     */
    public static TIntIntMap copy(TIntIntMap a) {
        TIntIntMap c = new TIntIntHashMap();
        TIntIntIterator iter = a.iterator();
        int i;
        for (i = 0; i < a.size(); ++i) {
            iter.advance();
            c.put(iter.key(), iter.value());
        }
        return c;
    }
    
    /**
     *
     @param a
     @return
     */
    public static float[][] convertToFloat(double[][] a) {
        float[][] c = new float[a.length][];
        int i, j;
        for (i = 0; i < a.length; ++i) {
            c[i] = new float[a[0].length];
            for (j = 0; j < a[0].length; ++j) {
                c[i][j] = (float) a[i][j];
            }
        }
        return c;
    }

    /**
     *
     @param a
     @return
     */
    public static float[] convertToFloat(double[] a) {
        float[] c = new float[a.length];
        int i;
        for (i = 0; i < a.length; ++i) {
            c[i] = (float) a[i];
        }
        return c;
    }

    /**
     *
     @param a
     @return
     */
    public static double[][] convertIntToDouble(int[][] a) {
        double[][] c = new double[a.length][];
        int i, j;
        for (i = 0; i < a.length; ++i) {
            c[i] = new double[a[0].length];
            for (j = 0; j < a[0].length; ++j) {
                c[i][j] = a[i][j];
            }
        }
        return c;
    }

    /**
     *
     @param a
     @return
     */
    public static double[][] convertToDouble(float[][] a) {
        double[][] c = new double[a.length][];
        int i, j;
        for (i = 0; i < a.length; ++i) {
            c[i] = new double[a[0].length];
            for (j = 0; j < a[0].length; ++j) {
                c[i][j] = a[i][j];
            }
        }
        return c;
    }

    /**
     *
     @param a
     @return
     */
    public static float[][] convertToFloat(int[][] a) {
        float[][] c = new float[a.length][];
        int i, j;
        for (i = 0; i < a.length; ++i) {
            c[i] = new float[a[0].length];
            for (j = 0; j < a[0].length; ++j) {
                c[i][j] = a[i][j];
            }
        }
        return c;
    }
    
    /**
     * create a permutation matrix given the vector of permuted element indexes.
     * Usage: pre-multiplying, P*A, results in permuting the rows of A.
     * post-multiplying, A*P, results in permuting the columns of A.
     @param assignments the permutation vector.
     @return permutation matrix
     */
    public static int[][] createPermutationMatrix(int[] assignments) {
        int n = assignments.length;
        int[][] p = new int[n][];
        int i;
        for (i = 0; i < n; ++i) {
            p[i] = new int[n];
        }
        for (i = 0; i < n; ++i) {
            p[i][assignments[i]] = 1;
        }
        assert(MatrixUtil.isAPermutationMatrix(p));
        return p;
    }

    /**
     * extract each column of a and append it to an output vector.
     * if a is [mxn], the output vector length will be m*n.
     @param a matrix a
     @return the output vector of stacked columns of a
     */
    public static double[] stack(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        double[] out = new double[m*n];
        int j;
        for (int i = 0; i < n; ++i) {
            for (j = 0; j < m; ++j) {
                out[i*m + j] = a[j][i];
            }
        }

        return out;
    }

    /**
     * Given two matrices A € R^(mxn) and B € R^(kxl), their Kronecker product,
     * denoted by A⨂B, is a new matrix € R^(mk x nl).
     *
     @param a matrix
     @param b matrix
     @return kronecker product of a and b
     */
    public static double[][] kroneckerProduct(double[][] a, double[][] b) {

        int m = a.length;
        int n = a[0].length;
        int k = b.length;
        int l = b[0].length;
        int mk = m * k;
        int nl = n * l;

        double[][] kron = MatrixUtil.zeros(mk, nl);
        int i;
        int j;
        int ii;
        double[][] b2;
        int r = 0;
        int c = 0;
        for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j) {
                // multiply each a[i][j] by entire matrix b2
                b2 = MatrixUtil.copy(b);
                MatrixUtil.multiply(b2, a[i][j]);
                // write the kXl items of b2 into [m*k][n*l] matrix starting at [i*k][j*l]
                for (ii = 0; ii < k; ++ii) {
                    System.arraycopy(b2[ii], 0, kron[i*k + ii], j*l, b2[ii].length);
                }

                // e.g. m=1,n=3,  k=3, l=4
                // B: 0,0  0,1  0,2  0,3    0,4 0,5 0,6 0,7  j=0:[0,l*j] j=1:[0,l*j]
                //    1, .. ..
                //    2, .. ..
                // B: 3, ...
            }
        }
        return kron;
    }

    /**
     * Given two vectors A € R^(m) and B € R^(k), their Kronecker product,
     * denoted by A⨂B, is a new matrix € R^(mk).  It is the outer product of the
     * vectors.
     @param a array a
     @param b array b
     @return the outer product
     */
    public static double[][] kroneckerProduct(double[] a, double[] b) {
        return outerProduct(a, b);
    }

    /**
     * return the nullspace of matrix A using SVD.  This method is more accurate, but slower than
     * using QR decomposition.
     @param a matrix of dimensions [mXn]
     @param tol tolerance as the equivalent to 0 which is approximately machine precision.
     @return nullspace of matrix A, transposed to row vectors for easier use.
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static double[][] nullSpaceUsingSVD(double[][] a, double tol) throws NotConvergedException {
        SVDProducts svd = MatrixUtil.performSVD(a);
        int rank = rank(svd);
        int n = a[0].length;
        // the last n-rank columns of V hold the nullspace
        double[][] nullspace = new double[n-rank][];
        for (int i = rank, c = 0; i < n; ++i, ++c) {
            nullspace[c] = Arrays.copyOf(svd.vT[i], svd.vT[i].length);
        }
        return nullspace;
    }

    /**
     * return the nullspace of matrix A using QR decomposition.
     * This method is less accurate, but faster than using SVD decomposition.
     @param a matrix of dimensions [mXn] where m must be .leq. n.
     @param tol tolerance as the equivalent to 0 which is approximately machine precision.
     @return nullspace of matrix A transposed to row vectors for easier use.
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static double[][] nullSpaceUsingQR(double[][] a, double tol) throws NotConvergedException {

        QR qr = QR.factorize(new DenseMatrix(a));
        UpperTriangDenseMatrix r = qr.getR();
        DenseMatrix q = qr.getQ();
        int rank = 0;
        int i;
        // r is square so can use the number of rows
        for (i = 0; i < r.numRows(); ++i) {
            if (r.get(i, i) > tol) {
                rank++;
            }
        }
        int n = a[0].length;
        double[][] qT = MatrixUtil.transpose(Matrices.getArray(q));
        // the last n-rank columns of q hold the nullspace
        double[][] nullspace = new double[n-rank][];
        int c;
        for (i = rank, c = 0; i < n; ++i, ++c) {
            nullspace[c] = Arrays.copyOf(qT[i], qT[i].length);
        }
        return nullspace;
    }

    /**
     * calculate the triple product as a x (b X c).
     <pre>
     from Boas "Mathematical Methods in the Physical Sciences", eqn 3.8
     a X (b X c) = (a dot c)*B - (a dot b) * c.

     also adapted from MASKS example code triple_product.m.
     "An introduction to 3-D Vision"
     by Y. Ma, S. Soatto, J. Kosecka, S. Sastry (MASKS)
     </pre>
     @param a array of length 3
     @param b array of length 3
     @param c array of length 3
     @return triple product of a, b, c.
     */
    public static double tripleProduct(double[] a, double[] b, double[] c) {

        return c[0]*(a[1]*b[2] - b[1]*a[2])
                + c[1]*(a[2]*b[0] - b[2]*a[0])
                + c[2]*(a[0]*b[1] - b[0]*a[1]);
    }


    public static int rank(SVD svd) {
        int rank = 0;
        TDoubleSet set = new TDoubleHashSet();
        // machine tolerance for 0:
        double tol = 1E-11;
        for (double s : svd.getS()) {
            if (Math.abs(s) > tol) {
                set.add(s);
            }
        }
        return set.size();
    }

    public static int rank(MatrixUtil.SVDProducts svd) {
        int rank = 0;
        TDoubleSet set = new TDoubleHashSet();
        // machine tolerance for 0:
        double tol = 1E-11;
        for (double s : svd.s) {
            if (Math.abs(s) > tol) {
                set.add(s);
            }
        }
        return set.size();
    }

    public static int rank(double[][] a) throws NotConvergedException {
        return rank(SVD.factorize(new DenseMatrix(a)));
    }

    /**
     * calculate for matrix A, the 2 eigenvalues.
     * The method uses the determinant, trace and quadratic formula.
     <pre>
     https://en.wikipedia.org/wiki/Eigenvalue_algorithm
     </pre>
     * @param a
     * @return
     */
    public static double[] eigenvalues2X2(double[][] a) {
        if (a.length != 2 || a[0].length != 2) {
            throw new IllegalArgumentException("a must be [2 X 2] in size");
        }
        double trA = trace(a);
        double detA = determinant(a);
        double gap = Math.sqrt(trA*trA - 4*detA);
        double eig1 = (trA + gap)/2;
        double eig2 = (trA - gap)/2.;
        return new double[]{eig1, eig2};
    }

    /**
     * calculate for matrix A, the 2 eigenvalues.
     * The method uses the determinant, trace and quadratic formula.
     <pre>
     https://en.wikipedia.org/wiki/Eigenvalue_algorithm
     </pre>
     * @param a
     * @return
     */
    public static float[] eigenvalues2X2(float[][] a) {
        if (a.length != 2 || a[0].length != 2) {
            throw new IllegalArgumentException("a must be [2 X 2] in size");
        }
        float trA = trace(a);
        float detA = determinant(a);
        float gap = (float) Math.sqrt(trA*trA - 4*detA);
        float eig1 = (trA + gap)/2;
        float eig2 = (trA - gap)/2.f;
        return new float[]{eig1, eig2};
    }
}
