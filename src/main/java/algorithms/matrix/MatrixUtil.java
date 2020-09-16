package algorithms.matrix;

import algorithms.matrix.LinearEquations.LUP;
import static algorithms.matrix.LinearEquations.LUPDecomposition;
import algorithms.misc.Misc0;
import gnu.trove.list.array.TDoubleArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 * TODO: tidy code for multiply and dot operations
 * 
 * @author nichole
 */
public class MatrixUtil {
    
    /**
     * multiply matrix m by vector n
     * @param m two dimensional array in row major format
     * @param n one dimensional array
     * @return the multiplication of matrix m by n
     */
    public static double[] multiply(double[][] m, double[] n) {

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
     * multiply matrix m by vector n
     * @param m two dimensional array in row major format
     * @param n one dimensional array
     * @return the multiplication of matrix m by n
     */
    public static float[] multiply(float[][] m, float[] n) {

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
        
        Can see example in this repository's related proejct called
        https://github.com/nking/curvature-scale-space-corners-and-transformations.git
        in the test class tests//algorithms/NetlibTest.java
        
        or in MTJ source code
        https://github.com/fommil/matrix-toolkits-java/blob/master/src/main/java/no/uib/cipr/matrix/DenseMatrix.java#L299
        
    }*/
    
    /**
     * multiply matrix m by matrix n
     * @param m tow dimensional array in ro major format
     * @param n two dimensional array in row major format
     * @return multiplication of m by n
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
        double[][] c = new double[mrows][];
        
        for (int mrow = 0; mrow < mrows; mrow++) {
            c[mrow] = new double[ncols];
            for (int ncol = 0; ncol < ncols; ncol++) {
                double sum = 0;                
                for (int mcol = 0; mcol < mcols; mcol++) {
                    sum += (m[mrow][mcol] * n[mcol][ncol]);                    
                }
                c[mrow][ncol] = sum;
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
     * perform dot product and return matrix of size mrows X ncols
     * @param m
     * @param n
     * @return 
     */
    public static DenseMatrix multiply(
        Matrix m, Matrix n) {

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
    
    
    public static double[] multiply(Matrix a, double[] b) {
        
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
    
    public static double multiplyByTranspose(TDoubleArrayList a, 
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
    
    public static double multiplyByTranspose(double[] a, 
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
    
    public static double dot(double[] a, double[] b) {

        if (a.length != b.length) {
            throw new IllegalArgumentException("a.length must == b.length");
        }
        
        double sum = 0;
        for (int i = 0; i < a.length; ++i) {
            sum += (a[i] * b[i]);
        }

        return sum;
    }
    
    public static double dot(int[] a, double[] b) {

        if (a.length != b.length) {
            throw new IllegalArgumentException("a.length must == b.length");
        }
        
        double sum = 0;
        for (int i = 0; i < a.length; ++i) {
            sum += (a[i] * b[i]);
        }

        return sum;
    }
    
    public static double[][] dot(DenseMatrix m1, DenseMatrix m2) {
        
        if (m1 == null) {
            throw new IllegalArgumentException("m1 cannot be null");
        }
        if (m2 == null) {
            throw new IllegalArgumentException("m2 cannot be null");
        }
        int cCols = m2.numColumns();
        int cRows = m1.numRows();
        
        if (m1.numColumns() != m2.numRows()) {
            throw new IllegalArgumentException(
                "the number of columns in m1 != number of rows in m2");
        }
        
        // m1 dot m2
        double[][] m = new double[cRows][cCols];
        for (int i = 0; i < cRows; i++) {
            m[i] = new double[cCols];
        }
        
        /*
        t00  t01  t02       x1  x2  x3  x4
        t10  t11  t12       y1  y2  y3  y4
        t20  t21  t22       1    1   1   1
        
        row=0, col=0:nCols0  times and plus col=0, row=0:nRows1 --> stored in row, row + (cAdd=0)
        row=1, col=0:nCols0  times and plus col=0, row=0:nRows1 --> stored in row, row + (cAdd=0)
                
        row=0, col=0:nCols0  times and plus col=(cAdd=1), row=0:nRows1 --> stored in row, row + (cAdd=0)
        */
        
        for (int colAdd = 0; colAdd < m2.numColumns(); colAdd++) {
            for (int row = 0; row < m1.numRows(); ++row) {
                for (int col = 0; col < m1.numColumns(); col++) {
                    double a = m1.get(row, col);
                    double b = m2.get(col, colAdd);
                    m[row][colAdd] += (a * b);
                }
            }
        }

        return m;
    }
        
    /**
     * apply dot operator to m1 and m2 which are formatted using same as 
     * DenseMatrix, that is row major [row][col].
     * @param m1
     * @param m2
     * @return 
     */
    public static double[][] dot(double[][] m1, double[][] m2) {
        
        if (m1 == null) {
            throw new IllegalArgumentException("m1 cannot be null");
        }
        if (m2 == null) {
            throw new IllegalArgumentException("m2 cannot be null");
        }
        int cCols = m2[0].length;
        int cRows = m1.length;
        
        if (m1[0].length != m2.length) {
            throw new IllegalArgumentException(
                "the number of columns in m1 != number of rows in m2");
        }
        
        // m1 dot m2
        double[][] m = new double[cRows][cCols];
        for (int i = 0; i < cRows; i++) {
            m[i] = new double[cCols];
        }
        
        for (int colAdd = 0; colAdd < m2[0].length; colAdd++) {
            for (int row = 0; row < m1.length; ++row) {
                for (int col = 0; col < m1[0].length; col++) {
                    double a = m1[row][col];
                    double b = m2[col][colAdd];
                    m[row][colAdd] += (a * b);
                }
            }
        }

        return m;
    }
    
    /**
     * calculate the product of vectors a and b as a*b^T.
     * @param a
     * @param b
     * @return outer product in double array of size [a.length][b.length])
     */
    public static double[][] outerProduct(double[] a, double[] b) {
        
        return vTv(a, b);
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
    
    /**
     * calculate the pseudo-inverse of matrix a, using the SVD of a,
     * specifically, V*R*U^T where R is 1/diagonal of S.
     * Note that if A^-1 exists, then the pseudo-inverse of A is equal to the
     * inverse of A.
     * 
     * Following Gilbert Strang's "Introduction to Linear Algebra".
     * 
     * TODO: read "ALTERNATIVE METHODS OF CALCULATION OF THE PSEUDO INVERSE
       OF A NON FULL-RANK MATRIX" by M. A. Murray-Lasso, 2008
       http://www.scielo.org.mx/pdf/jart/v6n3/v6n3a4.pdf
     * @param a
     * @return
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
        DenseMatrix uM = svd.getU();
        double[][] uT = MatrixUtil.convertToRowMajor(uM);
        
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
     * rank matrix, i.e. rank = m, using LUP decomposition
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
     * @return
     * @throws NotConvergedException 
     */
    public static double[][] pseudoinverseFullRank(double[][] a) throws NotConvergedException {
        int m = a.length;
        int n = a[0].length;
        
        /*
        NOTE: 
        for case m < n: 
           should calculate A† = A^T * (A * A^T)^(−1)
        
        for case n < m: 
           should calculate A† = (A^T * A)^(−1) * A^T
        */
        
        // limit for a number to be significant above 0 (precision of computer)
        double eps = 1e-16;
        
        //from cormen et al: A_pseudoinverse = inverse(A^T*A) * A^T
        double[][] aT = MatrixUtil.transpose(a);
        double[][] aTA = MatrixUtil.multiply(aT, a);
        DenseMatrix aTAM = new DenseMatrix(aTA);
        
        //NOTE that (A^T*A) has to be invertible, that is, the reduced echelon form
        // of A has linearly independent columns (no free variables, only pivots.
        // which also means rank==n).
        //
        // could invert (A^T*A) using the cofactor matrix/determinant
        //   or a convenience method from MTJ Matrix.solve
        
        DenseMatrix I = Matrices.identity(aTA[0].length);
        DenseMatrix identity = I.copy();
                
        // A.solve(Matrix B, Matrix X) solves X = A\B.
        // A*(A^-1) = I  ... identity = aTA / I
        DenseMatrix aTAIM = (DenseMatrix)aTAM.solve(I, identity);
        double[][] aTAI = MatrixUtil.convertToRowMajor(aTAIM);
        double[][] inv = MatrixUtil.multiply(aTAI, aT);
               
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
     * @param i
     * @param i0
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
     * the outer product of vectors v1 and v2 as the
     * transpose of v1 times v2.
     * @param v1
     * @param v2
     * @return the outer product of v1 and v2 as double array of 
     * size v1.length X v2.length.
     */
    public static double[][] vTv(double[] v1, double[] v2) {
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
     * @param a
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
            
            z = MatrixUtil.multiply(a, v);
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
     * @param a
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
                        
            z = MatrixUtil.multiply(a, v);
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
     * @param a
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
     * @param a
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
            
            x2 = vTv(x, x);
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
        DenseMatrix u = (DenseMatrix) svd.getU();
        
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
        from Strang "Linear Algebra":
        
        matrix A is positive definite for every non-zero vector x if x^T*A*x > 0 
        
        when a symmetric 2x2 matrix has 1 of these 4, it has all 4:
            1) both eigenvectors are positive
            2) the 1x1 and 2x2 ... determinants are positive a>0 and a*c-b^2>0
            3) the pivots are positive a>0 and a*c-b^2>0
            4) the function x^T*A*x is positive except at (0.0)

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
     * @return 
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
    
}
