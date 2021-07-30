package algorithms.matrix;

import algorithms.matrix.LinearEquations.LUP;
import algorithms.misc.Misc0;
import algorithms.util.FormatArray;
import gnu.trove.list.array.TDoubleArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Random;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.LowerSymmDenseMatrix;
import no.uib.cipr.matrix.LowerTriangDenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.QR;
import no.uib.cipr.matrix.SVD;
import no.uib.cipr.matrix.UpperTriangDenseMatrix;

/**
 
 <pre>
   misc notes:
 
   The eigenvalues can be determined in a few ways depending upon the matrix:
       If A is a positive definite matrix, can use the power method.  
          Caveat is that it performs best when the spectral gap is large 
          (diff between largest and 2nd largest eigenvalues).
       If A is symmetric, can diagonalize A.
       If A is not symmetric, can either diagonalize A^T*A and use the square root 
          of those eigenvalues or can use Singular Value Decomposition on A 
          (the singular values are the eigenvalues of A).
   Note that the eigenvectors are the same for the diagonalization of A, the 
   diagonalization of A^T*A, the SVD(A), and the for same operations performed 
   on the CUR-Decompositions of A (=C) or on C^T*C.
 </pre>
 * @author nichole
 */
public class MatrixUtil {
    
    /**
     * multiply the row vector v by matrix m.
     * @param v
     * @param m
     * @return result is size [1][m[0].length]
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
     * @param m two dimensional array in row major format
     * @param n one dimensional array
     * @return vector of length m.length
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
     * multiply matrix m by vector n and return results in given vector out
     * @param m two dimensional array in row major format
     * @param n one dimensional array
     * @param out vector of length m.length to return results in
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
        if (n.toString().equals(out.toString())) {
            throw new IllegalArgumentException("n cannot be the same as out");
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
     * @param m two dimensional array in row major format
     * @param n one dimensional array
     * @return 
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
     * @param m one dimensional array that is input and output for result
     * @param factor factor to multiply m by
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
     * @param a two dimensional array in that is input and output for result
     * @param m factor to multiply m by
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
     * @param m tow dimensional array in row major format
     * @param n two dimensional array in row major format
     * @return multiplication of m by n.  resulting matrix is size mrows X ncols.
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
     * @param m tow dimensional array in row major format
     * @param n two dimensional array in row major format
     * @param out the results of multiplication of m by n.  the matrix should be size mrows X ncols.
     */
    public static void multiply(double[][] m, double[][] n, double[][] out) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        // identity check:
        if (out.toString().equals(m.toString()) || out.toString().equals(n.toString())) {
            throw new IllegalArgumentException("out must be a different object than n and m");
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
    
   
    public static double[][] createATransposedTimesA(double[][] a) {

        if (a == null || a.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int m = a.length;
        int n = a[0].length;
        
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
        double[][] c = new double[n][n];
        for (i = 0; i < n; ++i) {
            c[i] = new double[n];
        }
        double sum;
        for (i = 0; i < n; i++) {
            for (outCol = 0; outCol < n; outCol++) {
                sum = 0;
                for (j = 0; j < m; j++) {
                    sum += (a[j][outCol] * a[j][i]);
                }
                c[outCol][i] = sum;
            }
        }

        return c;
    }
    
    public static void multiply(double[] a, double f) {
        for (int i = 0; i < a.length; ++i) {
            a[i] *= f;
        }
    }
    
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
     * @param m
     * @param n
     * @return 
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
    
    public static void multiply(TDoubleArrayList a, double f) {
        for (int i = 0; i < a.size(); ++i) {
            a.set(i, f * a.get(i));
        }
    }
    
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
     * @param m
     * @param diag
     * @return 
     */
    public static double[][] multiplyByDiagonal(
        double[][] m, double[] diag) {

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
        
        /*
        a b c      p0 0  0
        d e f      0  p1 0
                   0  0  p2        
        */
        
        double[][] c = new double[mrows][mcols];
        
        for (int row = 0; row < mrows; row++) {
            c[row] = new double[mcols];
            for (int mcol = 0; mcol < mcols; mcol++) {
                c[row][mcol] = m[row][mcol] * diag[mcol];
            }            
        }

        return c;
    }
    
    /**
     * perform dot product of m and a diagonalized matrix of diag,
     * and return matrix of size mrows X mcols
     * @param m
     * @param diag
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
     * @param a
     * @param b
     * @return scale result of a^T * b
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
     * @param a
     * @param b
     * @return scalar result of a^T * b
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
    
    public static double[] subtract(double[] m, double[] n) {

        int len = m.length;

        double[] c = new double[len];

        subtract(m, n, c);

        return c;
    }
    
    /**
     * calculate m[i] - s for each i=[0, m.length).
     * @param m
     * @param s
     * @return 
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
     * @param m
     * @param s
     * @return 
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
     * @param m
     * @param s
     * @return 
     */
    public static double[] add(double s, double[] m) {

        int len = m.length;

        double[] c = new double[len];

        for (int i = 0; i < len; ++i) {
            c[i] = s + m[i];
        }

        return c;
    }
    
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
    
    public static void add(int[] m, int n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int len = m.length;
             
        for (int i = 0; i < len; i++) {
            m[i] += n;
        }
    }
    
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
     * transpose matrix m into out matrix which must be size m[0].length X m.length.
     * @param m matrix to transpose
     * @param out output matrix to hold transposed m
     * @return 
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
    
    public static double[][] convertToRowMajor(DenseMatrix a) {
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
    
    public static double[][] convertToRowMajor(UpperTriangDenseMatrix a) {
        int nc = a.numColumns();
        int nr = a.numRows();
        //TODO: read API to see if can make this more effient by only visiting upper right of matrix
        double[][] out = new double[nr][];
        for (int i = 0; i < nr; ++i) {
            out[i] = new double[nc];
            for (int j = 0; j < nc; ++j) {
                out[i][j] = a.get(i, j);
            }
        }
        return out;
    }
    
    public static double[][] convertToRowMajor(LowerSymmDenseMatrix a) {
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
    
    /**
     * calculate the pseudo-inverse of matrix a for cases when the rank of a
     * is less than the width of matrix a.
     * Note that the term rank can be deceptive for cases when the 
     * original matrix a was rank deficient and then is perturbed to become
     * a full rank matrix leading to possibility of larger errors when treated as full rank.
     * the term rank deficient can be replaced by the ter numerically rank deficient.
     * (see Bjork 1991 Section 2, "Algorithms for linear least squares problems"
     * and Chap 6 of Golub & Van Loan).
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
     * @param a and m X n matrix
     * @return matrix with dimensions of a^T
     * @throws NotConvergedException 
     */
    public static double[][] pseudoinverseRankDeficient(double[][] a) throws NotConvergedException {
        int m = a.length;
        int n = a[0].length;
        
        // limit for a number to be significant above 0
        double eps = 1e-15;
        
        // from Gilbert Strang's "Introduction to Linear Algebra":
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
        
        double[][] sInverse = new double[n][];
        for (int i = 0; i < n; ++i) {
            sInverse[i] = new double[m];
        }
        
        for (int i = 0; i < s.length; ++i) {
            double sI = s[i];
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
     * rank matrix or overdetermined, i.e. rank = m or m > n, using LUP decomposition.
     * A_pseudoinverse = inverse(A^T*A) * A^T
     * <pre>
       for case m .lt. n: 
           should calculate A† = A^T * (A * A^T)^(−1)
        
        for case n .lt. m: 
           should calculate A† = (A^T * A)^(−1) * A^T
     </pre>
     * NOTE that (A^T*A) (or (A * A^T)) has to be invertible, that is, 
     * the reduced echelon form of A has linearly independent columns (rank==n).
     * following pseudocode from Cormen et al. Introduction to Algorithms.
     * @param a two dimensional array in row major format with dimensions
     * m x n.  a is a full-rank matrix.
     * a is a non-singular matrix(i.e. has exactly one solution). 
     * 
     * @return matrix of size [a[0].length][a.length]
     * @throws NotConvergedException 
     */
    public static double[][] pseudoinverseFullRank(double[][] a) throws NotConvergedException {
        int m = a.length;
        int n = a[0].length;
        
        // limit for a number to be significant above 0 (precision of computer)
        double eps = 1e-16;
        
        //from cormen et al: A_pseudoinverse = inverse(A^T*A) * A^T
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
        //     x = A\B solves the system of linear equations A*x = B
        //        can substitute A = A^T*A,   X = (A^T*A)^-1,   B = I
        //  X = A\B in MTJ is X = A.solve(B, X), that is, inputs are A and B.        
        DenseMatrix aTAI = (DenseMatrix)aTA.solve(I, placeholder);
        double[][] _aTAI = MatrixUtil.convertToRowMajor(aTAI);
        // _pseudoinverse = inverse(A^T*A) * A^T
        double[][] inv = MatrixUtil.multiply(_aTAI, _aT);
               
        return inv;
    }
    
    /**
      from Strang "Introduction to Linear Algebra":
      <pre>
       an inverse matrix may or may not exist.  
       (1) has to be a square matrix.
           A^-1 x A = I, where I is the identity matrix.
       (2) an inverse matrix has n pivots remaining after elimination,
                where pivot is the leftmost non-zero variable.
       (3) after elimination, next test for possible invertibility is 
             that the determinant is not zero

       if A is invertible, then A * x = b can be solved as x = A^-1 * b
       and (A * B)^-1 = B^-1 * A^-1
      </pre>
     NOTE: because this uses decomposition, each application using it should decide whether
     * to perform the exterior operations at same time to avoid recomputing
     * any matrices.
     * @param a
     * @return true if is invertible
     */
    public static boolean isInvertible(double[][] a) {
        
        int m = a.length;
        int n = a[0].length;
        if (m != n) {
            throw new IllegalArgumentException("a is not a squae matrix");
        }
        
        //rank:
        // -- could use PackCholesky followed by check for positive definiteness
        // -- could use LUP decomposition andcount the L diagonal 1's
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
     * class to hold the results of the Singular Value Decomposition
     */
    public static class SVDProducts {
        
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
     * @param a
     * @return 
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
     * @param a
     * @return 
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
        return out;
    }
    
    /**
     * create matrix A^T*A then perform SVD on it.  NOTE that the singular values
     * returned in S will have the square of values of SVD(A).s.
     * @param a
     * @return
     * @throws NotConvergedException 
     */
    public static SVDProducts performSVDATransposeA(double[][] a) throws NotConvergedException {
        double[][] aTa = MatrixUtil.createATransposedTimesA(a);
        return performSVDATransposeA(new DenseMatrix(aTa));
    }
    
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
     * A = Q*R of an orthogonal matrix Q and an upper triangular matrix R.
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
     * @param a a square or rectangular matrix.
     * @return 
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
     * @param xy
     * @return a matrix for use for canonical transformation of the points.
     * the format of the result is 
     * <pre>
     *  t[0] = new double[]{scale,       0,     -centroidX*scale};
        t[1] = new double[]{0,           scale, -centroidY*scale};
       </pre>
     */
    public static double[][] calculateNormalizationMatrix2X3(double[][] xy) {
        
        int nRows = xy.length;
        
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
    
     public static float[][] copy(float[][] a) {

        float[][] c = new float[a.length][a[0].length];

        for (int i = 0; i < a.length; ++i) {
            c[i] = new float[a[0].length];
            System.arraycopy(a[i], 0, c[i], 0, a[0].length);
        }
        
        return c;
    }

    public static double[][] copy(double[][] a) {
        
        double[][] m = new double[a.length][];
        
        for (int i = 0; i < a.length; ++i) {
            int n0 = a[i].length;
            m[i] = new double[n0];
            System.arraycopy(a[i], 0, m[i], 0, n0);
        }
        
        return m;
    }
    
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
     * @param a
     * @param row0 beginning index, inclusive
     * @param row1 end index, inclusive
     * @param col0 beginning index, inclusive
     * @param col1 end index, inclusive
     * @return 
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
     * @param a
     * @param row0 beginning index, inclusive
     * @param row1 end index, inclusive
     * @param col0 beginning index, inclusive
     * @param col1 end index, inclusive
     * @param out output matrix to hold the copied section
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
     * @param a
     * @param col index of column to extract
     * @return one dimensional array holding the column col of a
     */
    public static double[] extractColumn(double[][] a, int col) {
                
        double[] m = new double[a.length];
        extractColumn(a, col, m);
        return m;
    }
    
    /**
     * 
     * @param a
     * @param col index of column to extract
     * @param out one dimensional array holding the column col of a
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
     */
    public static double determinant(Matrix m) {

        double[][] a = no.uib.cipr.matrix.Matrices.getArray(m);

        return determinant(a);
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
     * 
     */
    public static double determinant(double[][] m) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (m.length != m[0].length) {
            throw new IllegalArgumentException("m must be a square");
        }
        if (m.length == 1) {
            return m[0][0];
        } else if (m.length == 2) {
            double s = ( m[0][0]*m[1][1] ) - ( m[0][1]*m[1][0] );
            return s;
        } else {
            double s = 0.0;
            // use 1st row as cofactors and minors
            for (int i = 0; i < m.length; i++) {

                double[][] n = copyExcept(m, i, 0);
                
                double tmp = m[i][0] * determinant(n);
                                
                if ((i & 1) == 0) {
                    s +=  tmp;
                } else {
                    s -=  tmp;
                }
            }
            return s;
        }
    }
    
    /**
     * create copy of matrix m except row and col
     * @param m
     * @param col
     * @param row
     * @return
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
     * v cross product with w is v X w = [v]_x * w.
     * It’s individual terms are a_j_i = -a_i_j.
       <pre>
       |    0   -v[2]  v[1] |
       |  v[2]    0    -v[0] |
       | -v[1]  v[0]    0  |
       </pre>
     * @param v
     * @return 
     */
    public static double[][] skewSymmetric(double[] v) {
        if (v.length != 3) { 
            throw new IllegalArgumentException("v.length must be 3");
        }
        
        
        double[][] out = new double[3][3];
        out[0] = new double[]{0,     -v[2],  v[1]};
        out[1] = new double[]{v[2],     0,  -v[0]};
        out[2] = new double[]{-v[1], v[0],    0};
        
        return out;
    }
    
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
    
    public static boolean areColinear(double[] p0, double[] p1, double eps) {
        if (p0.length != 3 || p1.length != 3) {
            throw new IllegalArgumentException("expecting p0 and p1 lengths to be 3");
        }
        double[] cp = crossProduct(p0, p1);
        double norm = lPSum(cp, 2);
        return Math.abs(norm) < eps;
    }
    
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
     * @param v
     * @return 
     */
    public static double[] normalizeL2(double[] v) {
        double sum = 0;
        for (double a : v) {
            sum += (a*a);
        }
        sum = Math.sqrt(sum);
        for (int i = 0; i < v.length; ++i) {
            v[i] /= sum;
        }
        return v;
    }
    
    /**
     * normalize each column of matrix a by the square root of the sum of 
     * its squared components. ||v||_2 for each column...
     * @param a matrix
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
                a[j][i] /= sum;
            }
        }
    }
    
    /**
     * normalize each row of matrix a by the square root of the sum of 
     * its squared components. ||v||_2 for each row...
     * @param a matrix
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
                a[i][j] /= sum;
            }
        }
    }
    
    /**
     * calculate the condition number as the largest singular value divided
     * by the singular value for i==(rank-1) of A where the singular values are
     * found using the SVD.
     * 
     * from https://blogs.mathworks.com/cleve/2017/07/17/what-is-the-condition-number-of-a-matrix/
     * A condition number for a matrix and computational task measures how 
     * sensitive the answer is to perturbations in the input data and to roundoff 
     * errors made during the solution process....If a matrix is singular, 
     * then its condition number is infinite.
     * ...(A large condition number means that the matrix is close to being singular).
     * 
     * 
     * @param a
     * @return 
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
     * determine the rank of a using the SVD.
     * @param a
     * @return
     * @throws NotConvergedException 
     */
    public static double rank(double[][] a) throws NotConvergedException {
        
        SVDProducts svd = performSVD(a);
        int rm1 = -1;
        for (int i = 0; i < svd.s.length; ++i) {
            if (svd.s[i] > 0) {
                rm1 = i;
            }
        }
        
        return rm1+1;
    }
    
    /**
     * summation = the (1/p) power of sum of 
     * its (components)^p.
     * @param v
     * @return 
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
     * calculate the Frobenius Norm of matrix a.
     * It's the square root of the sum of squares of each element.
     * It can also be calculated as the square root of the
     * trace of a_conjugate*a, or as
     * the square root of the sums of the squares of the singular
     * values of a.
     * @param a
     * @return 
     */
    public static double frobeniusNorm(double[][] a) {
        double sum = 0;
        int i, j;
        for (i = 0; i < a.length; ++i) {
            for (j = 0; j < a[i].length; ++j) {
                sum += (a[i][j]*a[i][j]);
            }
        }
        return sum;
    }
    
    public static double spectralNorm(double[][] r) throws NotConvergedException {
        SVDProducts svd = performSVDATransposeA(r);
        double norm = Math.sqrt(svd.s[0]);
        return norm;
    }
    
    /**
     * Given a symmetric matrix and a nonnegative number eps, find the
     * nearest symmetric positive semidefinite matrices with eigenvalues at least eps.
     * <pre>
     * References:
     * https://nhigham.com/2021/01/26/what-is-the-nearest-positive-semidefinite-matrix/
     * Cheng and Higham, 1998
     * 
     * </pre>
     * @param a a symmetric matrix
     * @param eps
     * @return 
     */
    public static double[][] nearestPositiveSemidefiniteToASymmetric(double[][] a,
        double eps) throws NotConvergedException {
        
        if (!isSquare(a)) {
            throw new IllegalArgumentException("a must b a square matrix");
        }
        
        // TODO: fix the test for symmetric and use it here
         
        // uses the Frobenius norm as a distinace.  notation: ||A||_F
        
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
        double[][] b = MatrixUtil.elementwiseAdd(a, MatrixUtil.transpose(a));
        MatrixUtil.multiply(b, 0.5);
        
        QH qh = performPolarDecomposition(b);
        
        //X_F:
        double[][] aPSD2 = MatrixUtil.elementwiseAdd(a, qh.h);
        MatrixUtil.multiply(aPSD2, 0.5);
        */
        
        return aPSD;
    }
    
    /**
     * Given a matrix a that is not necessarily symmetric,
     * and a nonnegative number eps, find the
     * nearest symmetric positive semidefinite matrices with eigenvalues at least eps.
     * <pre>
     * References:
     * https://nhigham.com/2021/01/26/what-is-the-nearest-positive-semidefinite-matrix/
     * 
     * </pre>
     * @param a a square matrix which can be non-symmetric.
     * @param eps
     * @return 
     */
    public static double[][] nearestPositiveSemidefiniteToA(double[][] a,
        double eps) throws NotConvergedException {
        
        // uses the Frobenius norm as a distinace
        
        if (!isSquare(a)) {
            throw new IllegalArgumentException("a must b a square matrix");
        }
        
        // TODO: fix the test for symmetric
                
        // uses the Frobenius norm as a distinace.  notation: ||A||_F
        
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
        
        // can use polar decompositon on square matrices:
                
        double[][] b = MatrixUtil.elementwiseAdd(a, MatrixUtil.transpose(a));
        MatrixUtil.multiply(b, 0.5);
        
        // H is formed from V, S, and V^T of SVD(B)
        QH qh = performPolarDecomposition(b);
        double[][] h = qh.h;
        
        //X_F:
        double[][] aPSD = MatrixUtil.elementwiseAdd(a, h);
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
            System.out.printf("before: evd eigenvalues=%s\n", FormatArray.toString(evd.getRealEigenvalues(), "%.5e"));

            //still needs a small perturbation to make the eigenvalues all >= eps
            // see https://nhigham.com/2021/02/16/diagonally-perturbing-a-symmetric-matrix-to-make-it-positive-definite/
        
            double[] e = Arrays.copyOf(eig, eig.length);
            Arrays.sort(e);
            double eigMin = Math.max(-e[0], eps);
            double[][] D = MatrixUtil.zeros(5, 5);
            for (i = 0; i < e.length; ++i) {
                D[i][i] = eigMin;
            }
            aPSD = MatrixUtil.elementwiseAdd(aPSD, D);
            evd = EVD.factorize(new DenseMatrix(aPSD));
            System.out.printf("after:  evd eigenvalues=%s\n", FormatArray.toString(evd.getRealEigenvalues(), "%.5e"));
        }
        
        return aPSD; 
    }
    
    /**
     * normalize vector v by power p, that is the (1/p) power of sum of 
     * its (components)^p.  notation is sometimes ||v||_p.
     * when p = 0, this is the manhattan normalization or taxi-cab normalization,
     * when p = 2, this is the euclidean normalization.
     * @param v
     * @return 
     */
    public static double[] normalizeLP(double[] v, double p) {
        double sum = 0;
        for (double a : v) {
            sum += Math.pow(a, p);
        }
        sum = Math.pow(sum, 1./p);
        for (int i = 0; i < v.length; ++i) {
            v[i] /= sum;
        }
        return v;
    }
    
    /**
     * the outer product of vectors v1 and v2, which is v1 as a single row matrix
     * and v2 as a single column matrix, so is v1 * v2^T.
     * @param v1
     * @param v2
     * @return the outer product of v1 and v2 as double array of 
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
     * TODO:consider implementing the inverse power method also to determine the
     * smallest eigenvalue and its eigenvector
     * @param a a positive definite matrix
     * @param nIterations
     * @return 
     */
    public static double powerMethod(double[][] a, int nIterations) {

        if (!isPositiveDefinite(a)) {
            throw new IllegalArgumentException("a must be positive definite");
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
     * @param a a positive definite matrix
     * @param tolerance iterations are stopped when the current multiplication vector
     * difference from previous is smaller than tolerance for each item.
     * @param x an initialized vector of size a.length that will be filled by
     * this method to hold the vector used to calculate eig = x^T * a * x
     * @return 
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
                v[row] = z[row] / eig;
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
     * @param a a positive definite matrix
     * @param tolerance iterations are stopped when the current multiplication vector
     * difference from previous is smaller than tolerance for each item.
     * @return 
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
     * @param a a positive definite matrix
     * @param tolerance iterations are stopped when the current multiplication vector
     * difference from previous is smaller than tolerance for each item.
     * @return an array of a.length eigenvectors
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
     * @param a a square symmetric non-negative definite matrix.
     * @throws no.uib.cipr.matrix.NotConvergedException
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
        
        int rank = 0;
        for (double sv : s) {
            if (sv > eps) {
                rank++;
            }
        }
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
        S is nxn with non-negative singular values.  rank is number of non-zero entries
        V is  mxn
        
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
        
        return j;
    }
    
    public static boolean isSquare(double[][] a) {
        return (a.length == a[0].length);
    }
    
    /*
     * NOTE: this method is a very rough look at symmetry.  
       matrix echelon row reduction should be performed 
       (a may have 0's) and isn't present yet.
     * @param a
     * @return 
    public static boolean isSymmetric(double[][] a) {
        if (!isSquare(a)) {
            return false;
        }
                
        // need reduced echelon form then can compare entries
        
        double tol = 1e-3;
        double r;
        int n = a.length;
        int i, j;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                if (Math.abs(a[i][j] - a[j][i]) > tol) {
                    return false;
                }
            }
        }
        return true;
    }
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
     * @param a subspace defined by A which is n columns of linearly 
     * independent vectors of length m (they are in real space R^m).
     * @param b 
     * @return 
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
        double[][] aPseudoInv = MatrixUtil.pseudoinverseFullRank(a);
        
        pr.x = MatrixUtil.multiplyMatrixByColumnVector(aPseudoInv, b);
        
        pr.p = MatrixUtil.multiplyMatrixByColumnVector(a, pr.x);
        
        pr.pMatrix = MatrixUtil.multiply(a, aPseudoInv);
        
        return pr;
    }
    
    /**
     * A matrix is positive definite if it’s symmetric and all its eigenvalues are positive
     * 
     * @param a
     * @return 
     */
    public static boolean isPositiveDefinite(double[][] a) {
        if (!isSquare(a)) {
            return false;
        }
        
        /*
        from Strang "Introduction to Linear Algebra", chapter 6.5.
        
        matrix A is positive definite for every non-zero vector x if x^T*A*x > 0 
        
        when a symmetric nxn matrix has 1 of these 4, it has all 4:
            1) both eigenvectors are positive
               Note: mathworks recommends a tolerance for error, so all eigenvalues > tolerance above 0.
            2) all upper left determinants (the 1x1 and 2x2 ... ) are positive
            3) the pivots are positive a>0 and a*c-b^2>0
            4) the function x^T * A * x is positive except at x = 0

        */
        // using rule that left upward rooted determinants > 0
        int n = a.length;
        double[][] sa;
        double det;
        for (int i = 0; i < n; ++i) {
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
     * find the equation for which A * A^(-1) = the identity matrix
     *
     *             1
     * A^(-1) =  ------ C^(T)  where C_ij = cofactor of a_ij
     *            det A
     *
     * @param m
     * @return
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
     * @param a
     * @return mean of each column as an array of size a[0].length
     */
    public static double[] mean(double[][] a) {
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
     * given a as vectors of data of nSamples of nVariables, return the
     * mean of each of the variables. 
     * note that the format must be a[nSamples][nVariables],
     * e.g. a[0] = [10, 100, 1000]', a[1] = [9, 101, 999]; for nSamples = 2
     * and nVariables = 3;
     * @param a
     * @return 
     */
    public static double[] standardDeviation(double[][] a) {
        int nSamples = a.length;
        int nVars = a[0].length;
        
        double[] c = mean(a);
        
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
     * @param a a square matrix.
     * @param v a vector of length of a.length.
     * @return the matrix a - v*I.
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
     * calculate the sum of the diagonal elements of a
     * @param a a square matrix.
     * @return the sum of the diagonal elements of a
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
     * calculate the sum of the diagonal elements of v*I (i.e. sum of all elements of v)
     * @param v a vector to be treated as diagonal elements of an identity matrix.
     * @return the sum of the elements of v
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
     * element-wise multiplication
     * @param a
     * @param b
     * @return 
     */
    public static double[][] elementwiseMultiplication(double[][] a, double[][] b) {
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
     * element-wise addition
     * @param a
     * @param b
     * @return 
     */
    public static double[][] elementwiseAdd(double[][] a, double[][] b) {
        double[][] out = MatrixUtil.zeros(a.length, a[0].length);
        elementwiseAdd(a, b, out);
        return out;
    }
    
    /**
     * element-wise addition
     * @param a
     * @param b
     * @param out the results of element wise add of a + b. Note that it
     * is safe to provide out as the same object as input argument a or b.
     */
    public static void elementwiseAdd(double[][] a, double[][] b, double[][] out) {
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
    
    /** element-wise subtraction
     * @param a
     * @param b
     * @return 
     */
    public static double[][] elementwiseSubtract(double[][] a, double[][] b) {
        double[][] out = MatrixUtil.zeros(a.length, a[0].length);
        elementwiseSubtract(a, b, out);
        return out;
    }
    
    /**
     * element-wise subtraction
     * @param a
     * @param b
     * @param out the results of element wise subtraction, a - b. Note that it
     * is safe to provide out as the same object as input argument a or b.
     */
    public static void elementwiseSubtract(double[][] a, double[][] b, double[][] out) {
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
     * element-wise subtraction
     * @param a
     * @param b
     * @param out the results of element wise subtraction, a - b. Note that it
     * is safe to provide out as the same object as input argument a or b.
     */
    public static void elementwiseSubtract(double[] a, double[] b, double[] out) {
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
     * element-wise multiplication
     * @param a
     * @param b
     * @return 
     */
    public static double[] elementwiseMultiplication(double[] a, double[] b) {
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
     * create an array of zeros
     * @param nRows
     * @param nCols
     * @return 
     */
    public static double[][] zeros(int nRows, int nCols) {
        double[][] out = new double[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            out[i] = new double[nCols];
            // java, by default, initializes with zeroes
        }
        return out;
    }
    
    public static double[][] createIdentityMatrix(int nRows) {
        double[][] out = new double[nRows][nRows];
        for (int i = 0; i < nRows; ++i) {
            out[i] = new double[nRows];
            out[i][i] = 1;
        }
        return out;
    }
    
    /**
     * right divide is element-wise division
     * @param a
     * @param b
     * @return 
     */
    public static double[][] elementwiseDivision(double[][] a, double[][] b) {
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
     * right divide is element-wise division, that is a[i]/b[i] for i = [0, a.length).
     * @param a
     * @param b
     * @return 
     */
    public static double[] elementwiseDivision(double[] a, double[] b) {
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
     * @param a 
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
     * @param a 
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
     * @param lowerTriangular the lower triangular matrix
     * @param b vector on the righthand side of the equation L*x=b
     * @return x in equation L*x = b
     */
    public static double[] forwardSubstitution(double[][] lowerTriangular, double[] b) {
        
        int m = b.length;
        
        if (lowerTriangular[0].length != m) {
            throw new IllegalArgumentException("the number of columns in "
                    + "lowerTriangular must equal the length of b");
        }
        
        double[] x = new double[m];
        
        int i, j;
        for (i = 0; i < m; ++i) {
            x[i] = b[i];
            for (j = 0; j <= (i-1); ++j) {
                x[i] -= (lowerTriangular[i][j]*x[j]);
            }
            x[i] /= lowerTriangular[i][i];
        }
        return x;
    }
    
    /**
     * solves for vector x in the equation L*x=b where L is the lower triangular
     * matrix and b is a vector.
     * runtime complexity is approx (b.length)^2.
     * @param lowerTriangular the lower triangular matrix
     * @param b vector on the righthand side of the equation L*x=b
     * @return x in equation L*x = b
     */
    public static double[] forwardSubstitution(LowerTriangDenseMatrix lowerTriangular, double[] b) {
        
        int m = b.length;
        
        if (lowerTriangular.numColumns() != m) {
            throw new IllegalArgumentException("the number of columns in "
                    + "lowerTriangular must equal the length of b");
        }
        
        double[] x = new double[m];
        
        int i, j;
        for (i = 0; i < m; ++i) {
            x[i] = b[i];
            for (j = 0; j <= (i-1); ++j) {
                x[i] -= (lowerTriangular.get(i, j)*x[j]);
            }
            x[i] /= lowerTriangular.get(i, i);
        }
        return x;
    }
    
    /**
     * solves for vector x in the equation U*x = y where 
     * U is an upper triangular matrix and y is a vector.
     * runtime complexity is approx (y.length)^2.
     * @param upperTriangular the upper triangular matrix
     * (a_i_j=0 where i>j)
     * <pre>
     *     0  1  2
        2  *  *  *
        1  *  *  
        0  *
           0  1  2
     * </pre>
     * @param y vector on righthand side of equation
     * @return x in equation U*x = y
     */
    public static double[] backwardSubstitution(double[][] upperTriangular, double[] y) {
        
        int m = y.length;
        
        if (upperTriangular[0].length != m) {
            throw new IllegalArgumentException("the number of columns in "
                    + "upperTriangular must equal the length of y");
        }
        
        double[] x = new double[m];
        int i, j;
        for (i = m-1; i >= 0; i--) {
            x[i] = y[i];
            for (j = i+1; j < m; ++j) {
                x[i] -= (x[j]*upperTriangular[i][j]);
            }
            x[i] /= upperTriangular[i][i];
        }
        
        return x;
    }
    
    /**
     * solves for vector x in the equation U*x = y where 
     * U is an upper triangular matrix and y is a vector.
     * runtime complexity is approx (y.length)^2.
     * @param upperTriangular the upper triangular matrix
     * @param y vector on righthand side of equation
     * @return x in equation U*x = y
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
            
    public static void fill(double[][] a, double value) {
        int i;
        for (i = 0; i < a.length; ++i) {
            Arrays.fill(a[i], value);
        }
    }
    
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
     * @param a
     * @return 
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
     * @param a and m X n matrix
     * @return matrix with dimensions of a^T
     * @throws NotConvergedException 
     */
    public static double[][] calculatePreconditionerFromLUP(double[][] a) throws NotConvergedException {
        
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
     * @param a and m X n matrix
     * @return matrix with dimensions of a^T
     * @throws NotConvergedException 
     */
    public static double[][] calculatePreconditionerFromColumns(double[][] a) throws NotConvergedException {
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
     * @param a square matrix
     * @return 
     */
    public static QH performPolarDecomposition(double[][] a) throws NotConvergedException {
        
        if (!isSquare(a)) {
            throw new IllegalArgumentException("a must b a square matrix");
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
}
