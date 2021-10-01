package algorithms.matrix;

import algorithms.matrix.MatrixUtil.ProjectionResults;
import algorithms.util.FormatArray;
import java.util.Arrays;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseCholesky;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.LowerSymmDenseMatrix;
import no.uib.cipr.matrix.LowerTriangDenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.QR;
import no.uib.cipr.matrix.UpperTriangDenseMatrix;
    
/**
 *
 * @author nichole
 */
public class MatrixUtilTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public MatrixUtilTest() {
    }
    
    public void testMultiply() throws Exception {

        double[][] a = new double[3][2];
        a[0] = new double[]{1, 2};
        a[1] = new double[]{2, 3};
        a[2] = new double[]{3, 4};

        double[] b = new double[]{4, 3};

        double[] m = MatrixUtil.multiplyMatrixByColumnVector(a, b);

        assertTrue( m[0] ==  10 );
        assertTrue( m[1] ==  17 );
        assertTrue( m[2] ==  24 );
        
        Arrays.fill(m, 0);
        MatrixUtil.multiplyMatrixByColumnVector(a, b, m);
        assertTrue( m[0] ==  10 );
        assertTrue( m[1] ==  17 );
        assertTrue( m[2] ==  24 );
        
        double[][] c = new double[3][];
        for (int i = 0; i < 3; i++) {
            c[i] = new double[3];
            for (int j = 0; j < 3; j++) {
                c[i][j] = 1;
            }
        }
        c[1][2] = 0;
        
           
        a = new double[2][3];
        a[0] = new double[]{1, 2, 3};
        a[1] = new double[]{2, 3, 4};
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0+...  a*p1...  a*p2...
        d*p0+...  d*p1...  d*p2...
        
        1 2 3    1 1 1
        2 3 4    1 1 0
                 1 1 1
        
        (1*1 + 2*1 + 3*1)  (1*1 + 2*1 + 3*1)  (1*1 + 0 + 3*1)
        (2*1 + 3*1 + 4*1)  (2*1 + 3*1 + 4*1)  (2*1 + 0 + 4*1)
        
        6  6 4
        9  9 6
        */
        
        double[][] d = MatrixUtil.multiply(a, c);

        assertTrue(d[0][0] == 6);
        assertTrue(d[0][1] == 6);
        assertTrue(d[0][2] == 4);
        assertTrue(d[1][0] == 9);
        assertTrue(d[1][1] == 9);
        assertTrue(d[1][2] == 6);
     
        MatrixUtil.fill(d, 0);
        MatrixUtil.multiply(a, c, d);
        assertTrue(d[0][0] == 6);
        assertTrue(d[0][1] == 6);
        assertTrue(d[0][2] == 4);
        assertTrue(d[1][0] == 9);
        assertTrue(d[1][1] == 9);
        assertTrue(d[1][2] == 6);
        
        /*
        example:  m is 1 2 3
                       4 5 6
                       7 8 9
        
                  n is 100  101  1
                       200  201  1
        
        multiply m by transpose of n:
        
        1 2 3     100  200
        4 5 6     101  201
        7 8 9      1   1
        
        (1*100 + 2*101 + 3*1)   (1*200 + 2*201 + 3*1)     305   605
        (4*100 + 5*101 + 6*1)   (4*200 + 5*201 + 6*1)  =  911  1811
        (7*100 + 8*101 + 9*1)   (7*200 + 8*201 + 9*1)    1517  3017
        */
        double[][] aa = new double[3][];
        aa[0] = new double[]{1, 2, 3};
        aa[1] = new double[]{4, 5, 6};
        aa[2] = new double[]{7, 8, 9};
        
        double[][] bb = new double[2][];
        bb[0] = new double[]{100, 101, 1};
        bb[1] = new double[]{200, 201, 1};
        
    }

    public void testTranspose() throws Exception {

        /*
        100  101  1    100  200
        200  201  1    101  201
                        1    1
        */
        float[][] bb = new float[2][];
        bb[0] = new float[]{100, 101, 1};
        bb[1] = new float[]{200, 201, 1};

        float[][] expected = new float[3][];
        expected[0] = new float[]{100, 200};
        expected[1] = new float[]{101, 201};
        expected[2] = new float[]{1, 1};

        float[][] cc = MatrixUtil.transpose(bb);
        for (int i = 0; i < cc.length; i++) {
            for (int j = 0; j < cc[i].length; j++) {
                assertTrue(expected[i][j] == cc[i][j]);
            }
        }
        
        float[][] dd = MatrixUtil.transpose(cc);

        for (int i = 0; i < dd.length; i++) {
            for (int j = 0; j < dd[i].length; j++) {
                assertTrue(bb[i][j] == dd[i][j]);
            }
        }
    }
    
    public void testTransposeD() throws Exception {

        /*
        100  101  1    100  200
        200  201  1    101  201
                        1    1
        */
        double[][] bb = new double[2][];
        bb[0] = new double[]{100, 101, 1};
        bb[1] = new double[]{200, 201, 1};

        double[][] expected = new double[3][];
        expected[0] = new double[]{100, 200};
        expected[1] = new double[]{101, 201};
        expected[2] = new double[]{1, 1};

        double[][] cc = MatrixUtil.transpose(bb);
        for (int i = 0; i < cc.length; i++) {
            for (int j = 0; j < cc[i].length; j++) {
                assertTrue(expected[i][j] == cc[i][j]);
            }
        }
        MatrixUtil.fill(cc, 0);
        MatrixUtil.transpose(bb, cc);
        for (int i = 0; i < cc.length; i++) {
            for (int j = 0; j < cc[i].length; j++) {
                assertTrue(expected[i][j] == cc[i][j]);
            }
        }

        double[][] dd = MatrixUtil.transpose(cc);

        for (int i = 0; i < dd.length; i++) {
            for (int j = 0; j < dd[i].length; j++) {
                assertTrue(bb[i][j] == dd[i][j]);
            }
        }
        
        bb = new double[2][];
        bb[0] = new double[]{100, 101, 1};
        bb[1] = new double[]{200, 201, 1};
         MatrixUtil.transpose(bb, cc);
        for (int i = 0; i < cc.length; i++) {
            for (int j = 0; j < cc[i].length; j++) {
                assertTrue(expected[i][j] == cc[i][j]);
            }
        }
    }
    
    public void testTranspose3() {
        int a0 = 4;
        int a1 = 9;
        int b0 = 2;
        int b1 = 3;
        double[][] a = new double[a0][];
        int i, j;
        for (i = 0; i < a0; ++i) {
            a[i] = new double[a1];
        }
        
        double v0 = 0;
        for (i = 0; i < a0; ++i) {
            for (j = 0; j < a1; ++j, v0++) {
                a[i][j] = v0;
            }
            System.out.printf("a[%d]=%s\n", i, FormatArray.toString(a[i], "%.0f"));
        }
        
        /*
         0,   1,  2, |  3,  4,  5,  |  6,  7,  8 
         9,  10, 11, | 12, 13, 14,  | 15, 16, 17 
        ----------------------------------------
         18, 19, 20, | 21, 22, 23,  | 24, 25, 26 
         27, 28, 29, | 30, 31, 32,  | 33, 34, 35 
        
        transposed:
        0  9  | 3  12  |  6 15 |
        1  10 | 4  13  |  7 16 |
        2  11 | 5  14  |  8 17 |
        ------------------------
        18 27 | 21 30  | 24 33 |
        19 28 | 22 31  | 25 34 |
        20 29 | 23 32  | 26 35 |
        
        0  9  | 18 27
        1  10 | 19 28
        2  11 | 20 29 
        --------
        3  12 | ...
        4  13 |
        5  14 |
        --------
        6 15  |
        7 16  |
        8 17  |
        */
        
        double[][] expected = MatrixUtil.transpose(a);
                
        BlockMatrixIsometric mA = new BlockMatrixIsometric(a, b0, b1);
        
        BlockMatrixIsometric mB = MatrixUtil.transpose(mA);
        assertEquals(expected.length, mB.getA().length);
        assertEquals(expected[0].length, mB.getA()[0].length);
        
        double diff;
        double tol = 1e-7;
        for (i = 0; i < expected.length; ++i) {
            for (j = 0; j < expected[0].length; ++j) {
                diff = Math.abs(expected[i][j] - mB.getA()[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        double[] c = MatrixUtil.reshapeToVector(mB.getA());
        assertEquals(expected.length * expected[0].length, c.length);
        for (i = 0; i < c.length; ++i) {
            diff = Math.abs(i - c[i]);
            assertTrue(diff < tol);
        }
    }
    
    public void testPseudoinverse() throws NotConvergedException {
        
        // from example 4, chapter 7 of Strang's Introduction to Lenear Algebra
        double[][] a = new double[2][];
        a[0] = new double[]{2, 2};
        a[1] = new double[]{1, 1};
        
        double eps = 1e-16;
        
        double[][] expected = new double[2][];
        expected[0] = new double[]{0.2, 0.1};
        expected[1] = new double[]{0.2, 0.1};
        
        double[][] inv = MatrixUtil.pseudoinverseRankDeficient(a);
        
        for (int i = 0; i < a[0].length; i++) {
            for (int j = 0; j < a.length; j++) {
                double t1 = expected[i][j];
                double t2 = inv[i][j];
                assertTrue(Math.abs(t2 - t1) < eps);
            }
        }
    }
    
    public void testPseudoinverse2() throws NotConvergedException {
        
        // from Cormen et al Introduction to Algorithms chap 28, near Fig 28.3
        double[][] a = new double[5][];
        a[0] = new double[]{1, -1, 1};
        a[1] = new double[]{1, 1, 1};
        a[2] = new double[]{1, 2, 4};
        a[3] = new double[]{1, 3, 9};
        a[4] = new double[]{1, 5, 25};
                
        double eps = 1e-16;
        eps = 0.01;
        
        double[][] expected = new double[3][];
        expected[0] = new double[]{0.5, 0.3, 0.2, 0.1, -0.1};
        expected[1] = new double[]{-0.388, 0.093, 0.19, 0.193, -0.088};
        expected[2] = new double[]{0.06, -0.036, -0.048, -0.036, 0.06};
        
        //from cormen et al: A_pseudoinverse = inverse(A^T*A) * A^T
        double[][] inv = MatrixUtil.pseudoinverseFullRank(a);
        
        double[] y = new double[]{2, 1, 1, 0, 3};
        
        double[] c = MatrixUtil.multiplyMatrixByColumnVector(inv, y);
        
        for (int i = 0; i < a[0].length; i++) {
            for (int j = 0; j < a.length; j++) {
                double t1 = expected[i][j];
                double t2 = inv[i][j];
                double diff = Math.abs(t2 - t1);
                assertTrue(diff < eps);
            }
        }
        
        a[0] = new double[]{1, -1, 1};
        a[1] = new double[]{1, 1, 1};
        a[2] = new double[]{1, 2, 4};
        a[3] = new double[]{1, 3, 9};
        a[4] = new double[]{1, 5, 25};
        double[][] inv2 = MatrixUtil.pseudoinverseRankDeficient(a);
                
    }
    
    public void testMultiply2() throws Exception {

        double[][] m1 = new double[2][3];
        m1[0] = new double[]{0, 1, 0};
        m1[1] = new double[]{1000, 100, 10};

        double[][] m2 = new double[3][2];
        m2[0] = new double[]{2, 1};
        m2[1] = new double[]{3, 0};
        m2[2] = new double[]{4, 0};
                
        double[][] m = MatrixUtil.multiply(m1, m2);

        assertTrue(m.length == 2);
        assertTrue(m[0].length == 2);
        assertTrue(m[1].length == 2);

        assertTrue(m[0][0] == 3);
        assertTrue(m[1][0] == 2340);
        assertTrue(m[0][1] == 0);
        assertTrue(m[1][1] == 1000);
    }
    
    public void testSubtract() throws Exception {

        double[] a = new double[]{100, 100, 100, 100};
        double[] b = new double[]{1, 2, 3, 4};

        double[] expected = new double[]{99, 98, 97, 96};

        double[] c = MatrixUtil.subtract(a, b);

        assertTrue(Arrays.equals(expected, c));
    }
    
    public void testDot() throws Exception {
       
        double[][] m1 = new double[2][3];
        m1[0] = new double[]{0, 1, 0};
        m1[1] = new double[]{1000, 100, 10};
        
        double[][] m2 = new double[3][2];
        m2[0] = new double[]{2, 1};
        m2[1] = new double[]{3, 0};
        m2[2] = new double[]{4, 0};
         
        /*
        0     1     0     2  1
        1000  100  10     3  0
                          4  0
        
        0*2    + 1*3   + 0*0     0*1    +  1*0   +  0*0
        1000*2 + 100*3 + 10*4    1000*1 +  100*0 + 10*0
        */
       
        double[][] m = MatrixUtil.convertToRowMajor(MatrixUtil.multiply(new DenseMatrix(m1), 
            new DenseMatrix(m2)));
        
        assertTrue(m.length == 2);
        assertTrue(m[0].length == 2);
        assertTrue(m[1].length == 2);
        
        assertTrue(m[0][0] == 3);
        assertTrue(m[1][0] == 2340);
        assertTrue(m[0][1] == 0);
        assertTrue(m[1][1] == 1000);
        
        double[] x = new double[]{1, 2, 3, 4, 5};
        double[] y = new double[]{5, 4, 3, 2, 1};
        double expected = 5*1 + 4*2+ 3*3 + 2*4 + 1*5;
        double dot = MatrixUtil.innerProduct(x, y);
        assertTrue(Math.abs(expected - dot) < 1.e-15);
        
        int[] xInt = new int[]{1, 2, 3, 4, 5};
        dot = MatrixUtil.innerProduct(xInt, y);
        assertTrue(Math.abs(expected - dot) < 1.e-15);
        
        //-----------------------
        double[][] expectedR = new double[3][3];
        expectedR[0] = new double[]{1, 3, 4};
        expectedR[1] = new double[]{0, 1, 9};
        expectedR[2] = new double[]{0, 0, 1};
        UpperTriangDenseMatrix rM = new UpperTriangDenseMatrix(new DenseMatrix(expectedR));
        double diff;
        double tol = 1e-5;
        
        double[][] rT = MatrixUtil.convertToRowMajor(rM);
        assertEquals(expectedR.length, rT.length);
        for (int i = 0; i < expectedR.length; ++i) {
            assertEquals(expectedR[i].length, rT[i].length);
            for (int j = 0; j < expectedR[i].length; ++j) {
                diff = Math.abs(expectedR[i][j] - rT[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        LowerSymmDenseMatrix rLM = new LowerSymmDenseMatrix(new DenseMatrix(expectedR));
        rT = MatrixUtil.convertToRowMajor(rM);
        assertEquals(expectedR.length, rT.length);
        for (int i = 0; i < expectedR.length; ++i) {
            assertEquals(expectedR[i].length, rT[i].length);
            for (int j = 0; j < expectedR[i].length; ++j) {
                diff = Math.abs(expectedR[i][j] - rT[i][j]);
                assertTrue(diff < tol);
            }
        }
    }
    
    public void testAdd() throws Exception {

        double[] a = new double[]{1, 2, 3, 4};
        double[] b = new double[]{100, 100, 100, 100};

        double[] expected = new double[]{101, 102, 103, 104};
        
        double[] c = MatrixUtil.add(a, b);
        
        assertTrue(Arrays.equals(expected, c));
    }
    
    public void testAdd2() throws Exception {

        float[] a = new float[]{1, 2, 3, 4};
        float[] b = new float[]{100, 100, 100, 100};

        float[] expected = new float[]{101, 102, 103, 104};
        
        float[] c = MatrixUtil.add(a, b);
        
        assertTrue(Arrays.equals(expected, c));
    }
    
    public void testAdd3() throws Exception {

        int[] a = new int[]{1, 2, 3, 4};
        int add = -1;
        
        int[] expected = new int[]{0, 1, 2, 3};
        
        MatrixUtil.add(a, add);
        
        assertTrue(Arrays.equals(expected, a));
    }
  
    public void testDeterminant() {

        /**
         | 1  -5  2 |
         | 7   3  4 |
         | 2   1  5 |

                  | 3 4 |         | 7 4 |         | 7 3 |
           =  1 * | 1 5 |  +  5 * | 2 5 |  +  2 * | 2 1 |  = 11 + 135 + 2 = 148

         */

        double[][] m = new double[3][3];
        m[0] = new double[]{1, -5, 2};
        m[1] = new double[]{7, 3, 4};
        m[2] = new double[]{2, 1, 5};

        double det = MatrixUtil.determinant(m);

        assertEquals(148., det);
        
        det = MatrixUtil.determinantFromLU(m);
        assertTrue(Math.abs(148. - det) < 1e-7);

        DenseMatrix _a = new DenseMatrix(m);
        Matrix aT = _a.transpose();
        Matrix a = aT.transpose();
        det = MatrixUtil.determinant(a);
        assertTrue(Math.abs(148. - det) < 1e-7);
    }
   
    public static void testPowerMethod() {
        
        // unit test from Strang "Linear Algebra"
        double[][] a = new double[2][2];
        a[0] = new double[]{0.9, 0.3};
        a[1] = new double[]{0.1, 0.7};
        
        double eEig = 1.0;
        
        double tol = 1.e-2;
        
        int nIter = 10;
        
        double eig = MatrixUtil.powerMethod(a, nIter);
        
        double diff = Math.abs(eig - eEig);
        
        System.out.printf("==> eig=%.10f diff=%.10f\n", eig, diff);
        
        assertTrue(diff < tol);
        
        eig = MatrixUtil.powerMethod(a, 1.e-3);
        diff = Math.abs(eig - eEig);
        System.out.printf("===> eig=%.10f diff=%.10f\n", eig, diff);        
        assertTrue(diff < tol);
    }
    
    public static void testPowerMethod2() {
        
        // unit test from Strang "Linear Algebra"
        double[][] a = new double[2][2];
        a[0] = new double[]{3, 2};
        a[1] = new double[]{2, 6};
        
        double eEig1 = 7.0;
        double eEig2 = 2.0;
                
        double tol = 1.e-2;
        
        double[] eigs = MatrixUtil.powerMethodEigenPairs(a, 1.e-3);
        
        double diff = Math.abs(eigs[0] - eEig1);
        
        assertTrue(diff < tol);
        
        diff = Math.abs(eigs[1] - eEig2);
        
        assertTrue(diff < tol);
        
    }
    
    public static void testSquareRoot() throws NotConvergedException {
        
        double[][] a = new double[2][2];
        a[0] = new double[]{33, 24};
        a[1] = new double[]{48, 57};
        //5, 2    4, 7
        double[][] j = MatrixUtil.squareRoot(a);
        assertEquals(2, j.length);
        double[][] jj = MatrixUtil.multiply(j, j);
        assertEquals(2, jj.length);
    }
    
    public static void testCopySubmatrix() {
        
        double eps = 1e-17;
        
        double[][] a = new double[3][3];
        a[0] = new double[]{0, 1, 2};
        a[1] = new double[]{3, 4, 5};
        a[2] = new double[]{6, 7, 8};
        
        int row0, row1, col0, col1;
        int nre, nce, i, j;
        double[][] c;
        double diff;
        
        row0=0; row1=0; nre=1;
        col0=0; col1=0; nce=1;        
        c = MatrixUtil.copySubMatrix(a, row0, row1, col0, col1);
        assertEquals(nre, c.length);
        assertEquals(nce, c[0].length);
        for (i = 0; i < nre; ++i) {
            for (j = 0; j < nce; ++j) {
                diff = Math.abs(a[row0 + i][col0 + j] - c[i][j]);
                assertTrue(diff <= eps);
            }
        }
        MatrixUtil.fill(c, 0);
        MatrixUtil.copySubMatrix(a, row0, row1, col0, col1, c);
        assertEquals(nre, c.length);
        assertEquals(nce, c[0].length);
        for (i = 0; i < nre; ++i) {
            for (j = 0; j < nce; ++j) {
                diff = Math.abs(a[row0 + i][col0 + j] - c[i][j]);
                assertTrue(diff <= eps);
            }
        }
        
        row0=1; row1=2; nre=2;
        col0=0; col1=2; nce=3;        
        c = MatrixUtil.copySubMatrix(a, row0, row1, col0, col1);
        assertEquals(nre, c.length);
        assertEquals(nce, c[0].length);
        for (i = 0; i < nre; ++i) {
            for (j = 0; j < nce; ++j) {
                diff = Math.abs(a[row0 + i][col0 + j] - c[i][j]);
                assertTrue(diff <= eps);
            }
        }
        
    }
    
     public void testInverse() throws Exception {
        /**
         * e.g.    | 1  0  -1 |
         *         |-2  1   0 |
         *         | 1 -1   2 |
         * 
         * det = 1
         * 
         * cofactor = 
         *            | (2)   -(-4)   (2-1) |
         *            | 1(-1)  (2+1) -(-1)  |
         *            |  1    -(0-2)    1   |
         * 
         *  cofactor^T = | 2  1  1 |
         *               | 4  3  2 |
         *               | 1  1  1 |
         * 
         *   inv = (1/det) * cofactor^T
         */
        System.out.println("testInverse");

        double[][] m = new double[3][3];
        m[0] = new double[]{1, -2, 1};
        m[1] = new double[]{0, 1, -1};
        m[2] = new double[]{-1, 0, 2};
        
        double[][] result = MatrixUtil.inverse(m);

        double[][] e = new double[3][3];
        e[0] = new double[]{2, 4, 1};
        e[1] = new double[]{1, 3, 1};
        e[2] = new double[]{1, 2, 1};

        assertTrue(result.length == e.length);
        assertTrue(result[0].length == e[0].length);
        
        for (int i = 0; i < e.length; i++) {
            double[] a = result[i];
            double[] b = e[i];
            assertTrue(Arrays.equals(a, b));
        }
    }
     
     public void testMean() {
         
         double[][] a = new double[3][2];
         a[0] = new double[]{10, 100};
         a[1] = new double[]{9, 110};
         a[2] = new double[]{11, 90};
         
         double[] expected = new double[]{10, 100};
         double[] mean = MatrixUtil.mean(a);
         
         assertEquals(expected.length, mean.length);
         
         double diff;
         double tol = 1.e-17;
         for (int i = 0; i < expected.length; ++i) {
             diff = expected[i] - mean[i];
             assertTrue(Math.abs(diff) < tol);
         }
         
         double[] expectedStDev = new double[]{1.0, Math.sqrt(200./2.)};
         double[] stdev = MatrixUtil.standardDeviation(a);
         for (int i = 0; i < expected.length; ++i) {
             diff = expectedStDev[i] - stdev[i];
             assertTrue(Math.abs(diff) < tol);
         }
     
     }
     
     public static void testAMinusVectorTimesIdentity() {
         
         double[][] a = new double[3][3];
         a[0] = new double[]{3, 2, 1};
         a[1] = new double[]{0, 0, 0};
         a[2] = new double[]{3, 2, 1};
         
         double[] v = new double[]{1, 1, -2};
         
         double[][] expected = new double[3][3];
         expected[0] = new double[]{2, 2, 1};
         expected[1] = new double[]{0, -1, 0};
         expected[2] = new double[]{3, 2, 3};
         
         double[][] result = MatrixUtil.aMinusVectorTimesIdentity(a, v);
         assertEquals(expected.length, result.length);
         assertEquals(expected[0].length, result[0].length);
         
         assertTrue(Arrays.equals(expected[0], result[0]));
         assertTrue(Arrays.equals(expected[1], result[1]));
         assertTrue(Arrays.equals(expected[2], result[2]));
     }
     
     public static void testTrace() {
         
         double[][] a = new double[3][3];
         a[0] = new double[]{3, 2, 1};
         a[1] = new double[]{0, 0, 0};
         a[2] = new double[]{3, 2, 1};
         
         double[] v = new double[]{1, 1, -1};
         
         double expectedA = 4;
         double expectedV = 1;
         
         double tol = 0.01;
         
         double resultA = MatrixUtil.trace(a);
         double resultV = MatrixUtil.trace(v);
         
         assertTrue(Math.abs(expectedA - resultA) < tol);
         assertTrue(Math.abs(expectedV - resultV) < tol);
     }
     
    
     public void testProjection() throws NotConvergedException {
         // from Chap 4 of "Introduction to Linear Algebra" by String
         double[][] a = new double[3][2];
         a[0] = new double[]{1, 0};
         a[1] = new double[]{1, 1};
         a[2] = new double[]{1, 2};
         
         double[] b = new double[]{6, 0, 0};
         
         ProjectionResults pr = MatrixUtil.projection(a, b);
         
         double[] expectedX = new double[]{5, -3};
         double[] expectedP = new double[]{5, 2, -1};
         double[][] expectedPMatrix = new double[3][3];
         expectedPMatrix[0] = new double[]{5/6., 2/6., -1/6.};
         expectedPMatrix[1] = new double[]{2/6., 2/6., 2/6.};
         expectedPMatrix[2] = new double[]{-1/6., 2/6., 5/6.};
         
         double tol = 0.01;
         
         assertEquals(expectedX.length, pr.x.length);
         assertEquals(expectedP.length, pr.p.length);
         assertEquals(expectedPMatrix.length, pr.pMatrix.length);
         
         int i, j;
         double diff;
         for (i = 0; i < expectedX.length; ++i) {
             diff = expectedX[i] - pr.x[i];
             assertTrue(Math.abs(diff) < tol);
         }
         for (i = 0; i < expectedP.length; ++i) {
             diff = expectedP[i] - pr.p[i];
             assertTrue(Math.abs(diff) < tol);
         }
         for (i = 0; i < expectedPMatrix.length; ++i) {
             for (j = 0; j < expectedPMatrix[i].length; ++j) {
                 diff = expectedPMatrix[i][j] - pr.pMatrix[i][j];
                 assertTrue(Math.abs(diff) < tol);
             }
         }
     }
     
     public void testCrossProduct() {
         double eps = 1e-6;
         //Boas, Section 4 Example from book "Mathematical Methods in the Physical Sciences"
         double[] p0 = new double[]{2, 1, -1};
         double[] p1 = new double[]{1, 3, -2};
         
         double[] expected = new double[]{1, 3, 5};
         
         double[] result = MatrixUtil.crossProduct(p0, p1);
         for (int i = 0; i < result.length; ++i) {
             assertTrue(Math.abs(result[i] - expected[i]) < eps);
         }
         
         // check right hand rule
         p0 = new double[]{0, 0, -1};
         p1 = new double[]{-1, 0, 0};
         expected = new double[]{0, 1, 0};
         result = MatrixUtil.crossProduct(p0, p1);
         for (int i = 0; i < result.length; ++i) {
             //System.out.printf("crossproduct [%d] %.6e\n", i, result[i]);
             assertTrue(Math.abs(result[i] - expected[i]) < eps);
         }
         result = MatrixUtil.normalizeL2(result);
         for (int i = 0; i < result.length; ++i) {
             //System.out.printf("normalized [%d] %.6e\n", i, result[i]);
             assertTrue(Math.abs(result[i] - expected[i]) < eps);
         }
         
         assertFalse(MatrixUtil.areColinear(p0, p1, eps));
         
         p0 = new double[]{0, 0, -1};
         p1 = new double[]{0, 0, -0.35};
         assertTrue(MatrixUtil.areColinear(p0, p1, eps));
     }
     
     public void testCreateATransposedTimesA() {
         double[][] a = new double[3][2];
         a[0] = new double[]{1, 2};
         a[1] = new double[]{4, 9};
         a[2] = new double[]{2, 7};
         
         /*
         1  4  2     1  2
         2  9  7     4  9
                     2  7
         1*1 + 4*4 + 2*2   1+2 + 
         */
         
         double[][] expected = new double[2][2];
         expected[0] = new double[]{1*1 + 4*4 + 2*2, 1*2 + 4*9 + 2*7};
         expected[1] = new double[]{2*1 + 9*4 + 7*2, 2*2 + 9*9 + 7*7};
         
         double[][] aTa = MatrixUtil.createATransposedTimesA(a);
         
         assertEquals(2, aTa.length);
         assertEquals(2, aTa[0].length);
         int i, j;
         double diff;
         double tol = 1e-7;
         for (i = 0; i < aTa.length; ++i) {
             for (j = 0; j < aTa[i].length; ++j) {
                 diff = Math.abs(aTa[i][j] - expected[i][j]);
                 assertTrue(diff < tol);
             }
         }
     }
     
     public void testIdentity() {
         int n = 3;
         double[][] i3 = MatrixUtil.createIdentityMatrix(n);
         assertEquals(n, i3.length);
         assertEquals(3, i3[0].length);
         for (int i = 0; i < n; ++i) {
             assertTrue(Math.abs(i3[i][i] - 1.) < 1.e-7);
         }
     }
     
     public void testExtractColumn() {
         
         double[][] a = new double[3][2];
         a[0] = new double[]{0, 5};
         a[1] = new double[]{9, 3};
         a[2] = new double[]{4, 9};
         
         double[] expected = new double[]{5, 3, 9};
         
         double[] r = MatrixUtil.extractColumn(a, 1);
         assertEquals(expected.length, r.length);
         
         double tol = 1.e-5;
         double diff;
         for (int i = 0; i < r.length; ++i) {
             diff = Math.abs(expected[i] - r[i]);
             assertTrue(diff < tol);
         }
         
         Arrays.fill(r, 0);
         MatrixUtil.extractColumn(a, 1, r);
         for (int i = 0; i < r.length; ++i) {
             diff = Math.abs(expected[i] - r[i]);
             assertTrue(diff < tol);
         }
     }
     
     public void testSkewSymmetric() throws NotConvergedException {
         
         // test s^T = s^-1
         // test a crossProduct b = [a]_x * b
         double[] a = new double[]{1, 7, 5};
         double[] b = new double[]{2, 3, 4};
         
         double[] aXb = MatrixUtil.crossProduct(a, b);
         double[][] s = MatrixUtil.skewSymmetric(a);
         double[] sb = MatrixUtil.multiplyMatrixByColumnVector(s, b);
         
         double diff;
         double tol = 1e-5;
         int i, j;
         
         double[][] sT = MatrixUtil.transpose(s);
         assertEquals(3, sT.length);
         for (i = 0; i < sT.length; ++i) {
             for (j = 0; j < sT[i].length; ++j) {
                 diff = Math.abs(sT[i][j] + s[i][j]);
                 assertTrue(diff < tol);
             }
         }
         
         assertEquals(3, sb.length);
         assertEquals(3, aXb.length);
         for (i = 0; i < sb.length; ++i) {
             diff = Math.abs(sb[i] - aXb[i]);
             assertTrue(diff < tol);
         }
     }
     
     public void testElementwiseAdd() {
         double[][] a = new double[3][3];
         double[][] b = new double[3][3];
         
         a[0] = new double[]{1, 2, 3};
         a[1] = new double[]{2, 2, 3};
         a[2] = new double[]{3, 2, 4};
         
         b[0] = new double[]{10, 0, 30};
         b[1] = new double[]{20, 10, 3};
         b[2] = new double[]{30, 20, 40};
         
         double[][] expected = new double[3][3];
         expected[0] = new double[]{11, 2, 33};
         expected[1] = new double[]{22, 12, 6};
         expected[2] = new double[]{33, 22, 44};
         
         double[][] r = MatrixUtil.elementwiseAdd(a, b);
         
         double eps = 1.e-11;
         double diff;
         int i, j;
         for (i = 0; i < 3; ++i) {
             for (j = 0; j < 3; ++j) {
                 diff = Math.abs(expected[i][j] - r[i][j]);
                 assertTrue(diff < eps);
             }
         }
         
         MatrixUtil.fill(r, 0);
         MatrixUtil.elementwiseAdd(a, b, r);
         for (i = 0; i < 3; ++i) {
             for (j = 0; j < 3; ++j) {
                 diff = Math.abs(expected[i][j] - r[i][j]);
                 assertTrue(diff < eps);
             }
         }
     }
     
     public void testElementwiseSubtract() {
         double[][] a = new double[3][3];
         double[][] b = new double[3][3];
         
         a[0] = new double[]{1, 2, 3};
         a[1] = new double[]{2, 2, 3};
         a[2] = new double[]{3, 2, 4};
         
         b[0] = new double[]{-10, 0, -30};
         b[1] = new double[]{-20, -10, -3};
         b[2] = new double[]{-30, -20, -40};
         
         double[][] expected = new double[3][3];
         expected[0] = new double[]{11, 2, 33};
         expected[1] = new double[]{22, 12, 6};
         expected[2] = new double[]{33, 22, 44};
         
         double[][] r = MatrixUtil.elementwiseSubtract(a, b);
         
         double eps = 1.e-11;
         double diff;
         int i, j;
         for (i = 0; i < 3; ++i) {
             for (j = 0; j < 3; ++j) {
                 diff = Math.abs(expected[i][j] - r[i][j]);
                 assertTrue(diff < eps);
             }
         }
         
         MatrixUtil.fill(r, 0);
         MatrixUtil.elementwiseSubtract(a, b, r);
         for (i = 0; i < 3; ++i) {
             for (j = 0; j < 3; ++j) {
                 diff = Math.abs(expected[i][j] - r[i][j]);
                 assertTrue(diff < eps);
             }
         }
         
         MatrixUtil.fill(r, 0);
         MatrixUtil.elementwiseSubtract(a[0], b[0], r[0]);
         j = 0;
         for (i = 0; i < 3; ++i) {
             diff = Math.abs(expected[0][i] - r[0][i]);
             assertTrue(diff < eps);
         }
     }
     
     public void testFlip() {
         double[][] a = new double[3][3];
         a[0] = new double[]{1, 2, 3};
         a[1] = new double[]{10, 20, 30};
         a[2] = new double[]{100, 200, 300};
         double[][] b = MatrixUtil.copy(a);
         
         double[][] expectedLR = new double[3][3];
         expectedLR[0] = new double[]{3, 2, 1};
         expectedLR[1] = new double[]{30, 20, 10};
         expectedLR[2] = new double[]{300, 200, 100};
         
         double[][] expectedUD = new double[3][3];
         expectedUD[2] = new double[]{1, 2, 3};
         expectedUD[1] = new double[]{10, 20, 30};
         expectedUD[0] = new double[]{100, 200, 300};
         
         double[][] a2 = new double[4][4];
         a2[0] = new double[]{1, 2, 3, 4};
         a2[1] = new double[]{10, 20, 30, 40};
         a2[2] = new double[]{100, 200, 300, 400};
         a2[3] = new double[]{1000, 2000, 3000, 4000};
         double[][] b2 = MatrixUtil.copy(a2);
         
         double[][] expectedLR2 = new double[4][4];
         expectedLR2[0] = new double[]{4, 3, 2, 1};
         expectedLR2[1] = new double[]{40, 30, 20, 10};
         expectedLR2[2] = new double[]{400, 300, 200, 100};
         expectedLR2[3] = new double[]{4000, 3000, 2000, 1000};
         
         double[][] expectedUD2 = new double[4][4];
         expectedUD2[0] = new double[]{1000, 2000, 3000};
         expectedUD2[1] = new double[]{100, 200, 300};
         expectedUD2[2] = new double[]{10, 20, 30, 40};
         expectedUD2[3] = new double[]{1, 2, 3, 4};
         
         MatrixUtil.flipLR(a);
         MatrixUtil.flipUD(b);
         MatrixUtil.flipLR(a2);
         MatrixUtil.flipUD(b2);
         
         double tol = 1.e-5;
         double diff;
         int i, j;
         for (i = 0; i < a.length; ++i) {
             for (j = 0; j < a[i].length; ++j) {
                 diff = Math.abs(expectedLR[i][j] - a[i][j]);
                 assertTrue(diff < tol);
                 diff = Math.abs(expectedUD[i][j] - b[i][j]);
                 assertTrue(diff < tol);
                 diff = Math.abs(expectedLR2[i][j] - a2[i][j]);
                 assertTrue(diff < tol);
                 diff = Math.abs(expectedUD2[i][j] - b2[i][j]);
                 assertTrue(diff < tol);
             }
         }
     }
     
     public void testForwardBackwardSubstitution() {
         /*
         a unit test from https://ece.uwaterloo.ca/~dwharder/NumericalAnalysis/04LinearAlgebra/cholesky/
         
         ---------------------------------------------------------
         
        L = [3 0 0 0; 0.2 4 0 0; -0.1 0.3 2 0; 0.5 -0.4 -0.2 5]
        L =
           3.00000   0.00000   0.00000   0.00000
           0.20000   4.00000   0.00000   0.00000
          -0.10000   0.30000   2.00000   0.00000
           0.50000  -0.40000  -0.20000   5.00000

        >> L * L'
        M = ans =

            9.0000   0.6000  -0.3000   1.5000
            0.6000  16.0400   1.1800  -1.5000
           -0.3000   1.1800   4.1000  -0.5700
            1.5000  -1.5000  -0.5700  25.4500
         
         b = (2.49, 0.566, 0.787, -2.209)T.
         using forward substitution to get y = (0.83, 0.1, 0.42, -0.5)T. 
         using backward substitution to get x = (0.3, 0, 0.2, -0.1)T.
         */
         int i, j;
         double diff;
         double tol = 1.e-7;
         
         /*
         double[][] _m2 = new double[4][];
         _m2[0] = new double[]{9.0000,   0.6000,  -0.3000,   1.5000};
         _m2[1] = new double[]{0.6000,  16.0400,   1.1800,  -1.5000};
         _m2[2] = new double[]{-0.3000,   1.1800,   4.1000,  -0.5700};
         _m2[3] = new double[]{1.5000,  -1.5000,  -0.5700,  25.4500};
         DenseCholesky chol = no.uib.cipr.matrix.DenseCholesky.factorize(new DenseMatrix(_m2));
         */
         
         double[][] _l2 = new double[4][];
         _l2[0] = new double[]{3.00000,  0.00000,   0.00000,   0.00000};
         _l2[1] = new double[]{0.20000,   4.00000,   0.00000,   0.00000};
         _l2[2] = new double[]{-0.10000,   0.30000,   2.00000,   0.00000};
         _l2[3] = new double[]{0.50000,  -0.40000,  -0.20000,   5.00000};
         
         double[][] _m2_r = MatrixUtil.multiply(_l2, MatrixUtil.transpose(_l2));
         System.out.printf("m2e=%s\n", FormatArray.toString(_m2_r, "%.3f"));
         
         double[] _b2 = new double[]{2.49, 0.566, 0.787, -2.209};
         double[] _y2 = new double[]{0.83, 0.1, 0.42, -0.5};
         double[] _x2 = new double[]{0.3, 0, 0.2, -0.1};
         
         double[] _y2_r = MatrixUtil.forwardSubstitution(_l2, _b2);
         System.out.printf("y2e=%s\ny2=%s\n", FormatArray.toString(_y2, "%.3f"),
             FormatArray.toString(_y2_r, "%.3f"));
         assertEquals(_y2.length, _y2_r.length);
         for (i = 0; i < _y2.length; ++i) {
             diff = Math.abs(_y2[i] - _y2_r[i]);
             assertTrue(diff < tol);
         }
         MatrixUtil.forwardSubstitution(_l2, _b2, _y2_r);
         for (i = 0; i < _y2.length; ++i) {
             diff = Math.abs(_y2[i] - _y2_r[i]);
             assertTrue(diff < tol);
         }
         double[] _x2_r = MatrixUtil.backwardSubstitution(MatrixUtil.transpose(_l2), _y2_r);
         System.out.printf("x2e=%s\nx2=%s\n", FormatArray.toString(_x2, "%.3f"),
             FormatArray.toString(_x2_r, "%.3f"));
         assertEquals(_x2.length, _x2_r.length);
         for (i = 0; i < _x2.length; ++i) {
             diff = Math.abs(_x2[i] - _x2_r[i]);
             assertTrue(diff < tol);
         }
         MatrixUtil.backwardSubstitution(MatrixUtil.transpose(_l2), _y2_r, _x2_r);
         for (i = 0; i < _x2.length; ++i) {
             diff = Math.abs(_x2[i] - _x2_r[i]);
             assertTrue(diff < tol);
         }
         
         double[] _b2_r = MatrixUtil.multiplyMatrixByColumnVector(_m2_r, _x2_r);
         System.out.printf("b2e=%s\nb2=%s\n", FormatArray.toString(_b2, "%.3f"),
             FormatArray.toString(_b2_r, "%.3f"));
         assertEquals(_b2.length, _b2_r.length);
         for (i = 0; i < _b2.length; ++i) {
             diff = Math.abs(_b2[i] - _b2_r[i]);
             assertTrue(diff < tol);
         }
         
         //-------------------
         LowerTriangDenseMatrix l2 = new LowerTriangDenseMatrix(new DenseMatrix(_l2));
         UpperTriangDenseMatrix l2T = new UpperTriangDenseMatrix(new DenseMatrix(
             MatrixUtil.transpose(_l2)));
         
         double[] _y2_r2 = MatrixUtil.forwardSubstitution(l2, _b2);
         //System.out.printf("y2e=%s\ny2=%s\n", FormatArray.toString(_y2, "%.3f"),
         //    FormatArray.toString(_y2_r, "%.3f"));
         assertEquals(_y2.length, _y2_r2.length);
         for (i = 0; i < _y2.length; ++i) {
             diff = Math.abs(_y2[i] - _y2_r2[i]);
             assertTrue(diff < tol);
         }
         double[] _x2_r2 = MatrixUtil.backwardSubstitution(l2T, _y2_r2);
         //System.out.printf("x2e=%s\nx2=%s\n", FormatArray.toString(_x2, "%.3f"),
         //    FormatArray.toString(_x2_r, "%.3f"));
         assertEquals(_x2.length, _x2_r2.length);
         for (i = 0; i < _x2.length; ++i) {
             diff = Math.abs(_x2[i] - _x2_r2[i]);
             assertTrue(diff < tol);
         }
         
         double[] _b2_r2 = MatrixUtil.multiplyMatrixByColumnVector(_m2_r, _x2_r2);
         //System.out.printf("b2e=%s\nb2=%s\n", FormatArray.toString(_b2, "%.3f"),
         //    FormatArray.toString(_b2_r, "%.3f"));
         assertEquals(_b2.length, _b2_r2.length);
         for (i = 0; i < _b2.length; ++i) {
             diff = Math.abs(_b2[i] - _b2_r2[i]);
             assertTrue(diff < tol);
         }
         
     }
     
     public void testCopy() {
        double[][] a = new double[2][];
        a[0] = new double[]{0, 1, 2, 3};
        a[1] = new double[]{4, 5, 6, 7};
        double[][] b = new double[2][];
        b[0] = new double[]{0, 10, 20, 30};
        b[1] = new double[]{40, 50, 60, 70};
        
        double[][] ac = MatrixUtil.copy(a);
        int i, j;
        double diff;
        double tol = 1e-11;
        for (i = 0; i < a.length; ++i) {
            for (j = 0; j < a[i].length; ++j) {
                diff = Math.abs(ac[i][j] - a[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        MatrixUtil.copy(b, ac);
        for (i = 0; i < b.length; ++i) {
            for (j = 0; j < b[i].length; ++j) {
                diff = Math.abs(ac[i][j] - b[i][j]);
                assertTrue(diff < tol);
            }
        }
     }
     
     public void testPSD() throws NotConvergedException {
        
         System.out.println("testPSD");
         
        // see https://nhigham.com/2020/07/21/what-is-a-symmetric-positive-definite-matrix/
        // when positive definite, all a[i][i] > 0
        //    and the eigenvalues of A are all positive
        //    ... many more
        // and given a vector x, x^T*A*x > 0
        //   the later happens when the diagonal has positive diagonal elements 
        //   and is diagonally dominant, that is, a[i][i] > summation_over_i_except_for_i==j( |a[i][j]| )
        double[][] a = new double[5][];
        a[0] = new double[]{0, 1, 2, 3, 4};
        a[1] = new double[]{1, 0, 1, 2, 3};
        a[2] = new double[]{2, 1, 0, 1, 2};
        a[3] = new double[]{3, 2, 1, 0, 1};
        a[4] = new double[]{4, 3, 2, 1, 0};
        
        double eps = 1.e-11;
        
        double[][] aPSD = MatrixUtil.nearestPositiveSemidefiniteToASymmetric(a, eps);
        
        double[][] aPSD2 = MatrixUtil.nearestPositiveSemidefiniteToA(a, eps);
        
        System.out.printf("a=\n%s\n", FormatArray.toString(a, "%.5e"));
        System.out.printf("aPSD=\n%s\n", FormatArray.toString(aPSD, "%.5e"));
        System.out.printf("aPSD2=\n%s\n", FormatArray.toString(aPSD2, "%.5e"));
        
        EVD evd1 = EVD.factorize(new DenseMatrix(a));
        EVD evd2 = EVD.factorize(new DenseMatrix(aPSD));
        EVD evd3 = EVD.factorize(new DenseMatrix(aPSD2));
        
        
        System.out.printf("eig(a)=\n%s\n", FormatArray.toString(evd1.getRealEigenvalues(), "%.5e"));
        System.out.printf("eig(aPSD)=\n%s\n", FormatArray.toString(evd2.getRealEigenvalues(), "%.5e"));
        System.out.printf("eig(aPSD2)=\n%s\n", FormatArray.toString(evd3.getRealEigenvalues(), "%.5e"));
        
        
        double[][] aMinusPSD = MatrixUtil.elementwiseSubtract(a, aPSD);
        double dist1 = MatrixUtil.frobeniusNorm(aMinusPSD);
        
        double[][] aMinusPSD2 = MatrixUtil.elementwiseSubtract(a, aPSD2);
        double dist2 = MatrixUtil.frobeniusNorm(aMinusPSD2);
        double[][] aPSDMinusPSD2 = MatrixUtil.elementwiseSubtract(aPSD, aPSD2);
        double dist3 = MatrixUtil.frobeniusNorm(aPSDMinusPSD2);
        
        System.out.printf("dist1=%.7e, dist2=%.7e, dist3=%.7e\n", dist1, dist2, dist3);
        
        
        double[][] g = LinearEquations.choleskyDecompositionViaLDL(aPSD, eps);
        double[][] g2 = LinearEquations.choleskyDecompositionViaLDL(aPSD2, eps);
        EVD evd4 = EVD.factorize(new DenseMatrix(g));
        System.out.printf("chol(aPSD)=\n%s\n", FormatArray.toString(g, "%.5e"));
        System.out.printf("chol(aPSD2)=\n%s\n", FormatArray.toString(g2, "%.5e"));
        
        System.out.printf("eig(g)=\n%s\n", FormatArray.toString(evd4.getRealEigenvalues(), "%.5e"));
        
    }

     public void testRank() throws NotConvergedException {
         double[][] a;
         double eps = 1e-7;
         int r;
         
         a = new double[2][];
         a[0] = new double[]{1, -2};
         a[1] = new double[]{3, -6};
         r = MatrixUtil.rank(a, eps);
         assertEquals(1, r);
         
         a = new double[2][];
         a[0] = new double[]{0, 2, 4};
         a[1] = new double[]{3,-2, 5};
         r = MatrixUtil.rank(a, eps);
         assertEquals(2, r);
         
         a = new double[3][];
         a[0] = new double[]{2, 4, -2, 2};
         a[1] = new double[]{4, 9, -3, 8};
         a[2] = new double[]{-2, -3, 7, 10};
         r = MatrixUtil.rank(a, eps);
         assertEquals(3, r);
         
         
     }
     
     public void testDiagonalization() throws Exception {
         
         System.out.println("testDiagonalization");
         
         double[][] a;
         
         a = new double[2][];
         a[0] = new double[]{0.5, 0.5};
         a[1] = new double[]{0.5, 0.5};
         
         MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(a);
         System.out.printf("svd(a).u=\n%s\n", FormatArray.toString(svd.u, "%.3e"));
         System.out.printf("svd(a).vT=\n%s\n", FormatArray.toString(svd.vT, "%.3e"));
         System.out.printf("svd(a).s=\n%s\n", FormatArray.toString(svd.s, "%.3e"));
         
         double[][] q_uvt = MatrixUtil.multiply(svd.u, svd.vT);
         System.out.printf("Q_uvt = svd(a).u * svd(a).vT=\n%s\n", FormatArray.toString(q_uvt, "%.3e"));
         
         // a must have independent columns, so that's violated.  just looking at results.
         // The columns of Q are the eigenvectors of A.
         // A*Q = Q * diag(eigenvalues of A).
         QR qr = QR.factorize(new DenseMatrix(a));
         System.out.printf("\nqr(a).q=\n%s\n", FormatArray.toString(
             Matrices.getArray(qr.getQ()), "%.3e"));
         System.out.printf("qr(a).r=\n%s\n", FormatArray.toString(
             Matrices.getArray(qr.getR()), "%.3e"));
         
         double[][] s, sInv, sInv2, delta, delta2, a2, a3;
         
         s = MatrixUtil.copy(svd.u);
         DenseMatrix _s = new DenseMatrix(s);
         System.out.printf("s = svd(a).u=\n%s\n", FormatArray.toString(s, "%.3e"));
         
         //x = A\B solves the system of linear equations A*x = B for x.
         //  X = A\B in MTJ is X = A.solve(B, X), that is, inputs are A and B.
         //  S*S^-1 = I
         //  S.solve(I, sInv) to get sInv using MTJ:
         DenseMatrix _sInv = new DenseMatrix(_s.numRows(), _s.numColumns());
             _sInv = (DenseMatrix) _s.solve(
             new DenseMatrix(MatrixUtil.createIdentityMatrix(s.length)), 
             _sInv);
         sInv = MatrixUtil.convertToRowMajor(_sInv);
         
         sInv2 = MatrixUtil.inverse(s);
         System.out.printf("sInv=S\\I:  S.solve(I, sInv)=\n%s\n", 
             FormatArray.toString(sInv, "%.3e"));
         System.out.printf("sInv2=MU.inverse(s)=\n%s\n", 
             FormatArray.toString(sInv2, "%.3e"));
                 
         // Delta = sInv * A * s 
         delta = MatrixUtil.multiply(MatrixUtil.multiply(sInv, a), s);
         System.out.printf("Delta = sInv * A * s=\n%s\n", 
             FormatArray.toString(delta, "%.3e"));
         
         delta2 = MatrixUtil.zeros(svd.s.length, svd.s.length);
         for (int i = 0; i < delta.length; ++i) {
             delta2[i][i] = svd.s[i];
         }
         // A = S * Delta * S^-1
         a2 = MatrixUtil.multiply(MatrixUtil.multiply(s, delta), sInv);
         a3 = MatrixUtil.multiply(MatrixUtil.multiply(s, delta2), sInv);
         System.out.printf(" A = S * Delta * S^-1\n%s\n", FormatArray.toString(a2, "%.3e"));
         System.out.printf(" A = S * Delta2(from svd) * S^-1\n%s\n", FormatArray.toString(a3, "%.3e"));
     }
     
     public void testNormalizeLP() {
         double eps= 1e-7;
         
         double[] a = new double[]{1, 2, 3, 4};
         double p = 3;
         double[] expected = new double[]{
             0.215443469003, 0.430886938006, 0.64633040701, 0.861773876013
         };
         double[] b = MatrixUtil.normalizeLP(a, p);
         int i;
         assertEquals(expected.length, b.length);
         for (i=0; i < b.length; ++i) {
             assertTrue(Math.abs(expected[i]-b[i]) < eps);
         }
         assertTrue(Arrays.equals(new double[]{1, 2, 3, 4}, a));
     }
     public void testNormalizeL2() {
         double eps= 1e-7;
         
         double[] a = new double[]{1, 2, 3, 4};
         double[] expected = new double[]{
             0.182574185835, 0.36514837167, 0.547722557505, 0.73029674334
         };
         double[] b = MatrixUtil.normalizeL2(a);
         int i;
         assertEquals(expected.length, b.length);
         for (i=0; i < b.length; ++i) {
             assertTrue(Math.abs(expected[i]-b[i]) < eps);
         }
         assertTrue(Arrays.equals(new double[]{1, 2, 3, 4}, a));
     }
}
