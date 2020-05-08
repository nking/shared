package algorithms.matrix;

import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 *
 * @author nichole
 */
public class Misc {
    
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
                "the number of columns in m must equal the number of rows in n");
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
     * calculate the pseudo-inverse of matrix a, using the SVD of a.
     * Note that if A^-1 exists, then the pseudo-inverse of A is equal to the
     * inverse of A.
     * Note also that this method returns valid answer for rows .lte. rank,
     * but for rows .gt. rank, the null-space needs to be solved, 
     * which is not done here.
     * (For an example on solving the null-space, that is AX=0, see
     * the class algorithms.imageProcessing.transform.EpipolarTransformer and
     * method calculateEpipolarProjectionFor7Points in
     * project curvature-scale-space-corners-and-transformations.
     * That method follows the 7-point algorithm in R. Hartley and A. Zisserman, 
     * "Multiple View Geometry in Computer Vision"
     * @param a
     * @return
     * @throws NotConvergedException 
     */
    public static double[][] pseudoinverse(double[][] a) throws NotConvergedException {
        int m = a.length;
        int n = a[0].length;
        
        // limit for a number to be significant above 0
        double eps = 1e-16;
        
        // from Gilbert Strang's "Introduction to Linear Algebra":
        // uses SVD:
        //   A is m X n and n >= m
        //
        //   A_inverse = V * pseudoinverse(S) * U^T
        //        where pseudoinverse(S) is simply an empty matrix with the diagonal
        //        being the reciprocal of each singular value
        //      NOTE: compare number of ops w/ this facoring:
        //          V * (pseudoinverse(S) * U^T)
        //   A_inverse is n X m
        //
        //   V is n X n
        //   pseudoinverse(S) is n X m
        //   U^T is m X m
        
        //NOTE: for the case of n > m or, knowing the rank r and n > r,
        //    there is null-space to solve for (AX = 0).
        //    a solution for a specific case of m=7 and n=9 is present in
        //    the computer vision algorithm for the 7-point solution in
        //    the book of Hartley & Zisserman.
        //    An implementation of that is in this code's sibling project
        //    called curvature-scale-space-corners-and-transformations in
        //    the class algorithms.imageProcessing.transform.EpipolarTransformer
        
        boolean transpose = (n > m);
        if (transpose) {
            a = Misc.transpose(a);
        }
        
        DenseMatrix aMatrix = new DenseMatrix(a);
        SVD svd = SVD.factorize(aMatrix);
        
        //TODO: rewrite to use fewer data structures and multiply in place.
        
        // s is an array of size min(m,n)
        double[] s = svd.getS();
        DenseMatrix vTM = svd.getVt();
        double[][] vT = Misc.convertToRowMajor(vTM);
        double[][] v = Misc.transpose(vT);
        DenseMatrix uM = svd.getU();
        double[][] uT = Misc.convertToRowMajor(uM);
        
        //   A_inverse is n X m
        //   V is n X n
        //   pseudoinverse(S) is n X m
        //   U^T is m X m
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
        
        //   A_inverse is n X m
        //   V is n X n
        //   pseudoinverse(S) is n X m
        //   U^T is m X m
        
        //A_inverse = V * pseudoinverse(S) * U^T
        double[][] inv = Misc.multiply(v, sInverse);
        inv = Misc.multiply(inv, uT);
        
        if (transpose) {
            inv = Misc.transpose(inv);
        }
        
        return inv;
    }
    
    /**
     * calculate the pseudo-inverse of matrix a (dimensions mxn) which is a full
     * rank matrix, i.e. rank = m, using LUP decomposition
     * (inverse(A^T*A) * A^T).
     * following pseudocode from Cormen et al. Introduction to Algorithms.
     * @param a two dimensional array in row major format with dimensions
     * m x n.  a is a full-rank matrix.
     * a is a non-singular matrix(i.e. has exactly one solution). 
     * 
     * NOTE: for a non-full-rank matrix, should implement this:
     * "ALTERNATIVE METHODS OF CALCULATION OF THE PSEUDO INVERSE
       OF A NON FULL-RANK MATRIX" by M. A. Murray-Lasso, 2008
       http://www.scielo.org.mx/pdf/jart/v6n3/v6n3a4.pdf
     
     * @return
     * @throws NotConvergedException 
     */
    public static double[][] pseudoinverse2(double[][] a) throws NotConvergedException {
        int m = a.length;
        int n = a[0].length;
        
        // limit for a number to be significant above 0 (precision of computer)
        double eps = 1e-16;
        
        //from cormen et al: A_pseudoinverse = inverse(A^T*A) * A^T
        double[][] aT = Misc.transpose(a);
        double[][] aTA = Misc.multiply(aT, a);
        DenseMatrix aTAM = new DenseMatrix(aTA);
        
        DenseMatrix I = Matrices.identity(aTA[0].length);
        DenseMatrix identity = I.copy();
                
        // A.solve(Matrix B, Matrix X) solves X = A\B.
        // A*(A^-1) = I  ... identity = aTA / I
        DenseMatrix aTAIM = (DenseMatrix)aTAM.solve(I, identity);
        double[][] aTAI = Misc.convertToRowMajor(aTAIM);
        double[][] inv = Misc.multiply(aTAI, aT);
               
        return inv;
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
    
}
