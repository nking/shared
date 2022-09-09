package algorithms.matrix;

import algorithms.matrix.LinearEquations.LDL;
import algorithms.matrix.LinearEquations.LDM;
import algorithms.matrix.LinearEquations.LU;
import algorithms.matrix.LinearEquations.LUP;
import algorithms.util.FormatArray;
import java.util.Arrays;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

public class LinearEquationsTest extends TestCase { 
    public LinearEquationsTest(String testName) {
        super(testName);
    }
    
    public void testSolveXFromLUDecomposition() {
        /* test from Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms, chap. 28*/
        double[][] a = new double[3][3];
        a[0] = new double[]{1, 2, 0};
        a[1] = new double[]{3, 4, 4};
        a[2] = new double[]{5, 6, 3};

        double[] b = new double[]{3, 7, 8};

        double[] expectedX = new double[]{-1.4, 2.2, 0.6};

        double tol = 1.e-6;
        double diff;

        double[] x = LinearEquations.solveXFromLUDecomposition(a, b);

        assertEquals(expectedX.length, x.length);

        for (int i = 0; i < x.length; ++i) {
            diff = Math.abs(x[i] - expectedX[i]);
            assertTrue(diff < tol);
        }
    }
            //public static double[] solveXFromLUDecomposition(double[][] a, double[] b) {

    public void testLUPSolveX() {

         /* test from Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms, chap. 28*/
         double[][] a = new double[3][3];
         a[0] = new double[]{1, 2, 0};
         a[1] = new double[]{3, 4, 4};
         a[2] = new double[]{5, 6, 3};

         double[] b = new double[]{3, 7, 8};

         double[][] ell = new double[3][3];
         ell[0] = new double[]{1, 0, 0};
         ell[1] = new double[]{0.2, 1, 0};
         ell[2] = new double[]{0.6, 0.5, 1};

         double[][] u = new double[3][3];
         u[0] = new double[]{5,  6,  3};
         u[1] = new double[]{0, 0.8, -0.6};
         u[2] = new double[]{0,  0, 2.5};

         int[] p = new int[]{2, 0, 1};

         double eps = 0.01;
         
         // verify that P*A = L*U
         double[][] permutation = new double[3][3];
         permutation[0] = new double[3];
         permutation[1] = new double[3];
         permutation[2] = new double[3];
         permutation[0][p[0]] = 1;
         permutation[1][p[1]] = 1;
         permutation[2][p[2]] = 1;
         double[][] pa = MatrixUtil.multiply(permutation, a);
         double[][] lu = MatrixUtil.multiply(ell, u);
         for (int i = 0; i < pa.length; ++i) {
             for (int j = 0; j < pa[i].length; ++j) {
                 double t1 = pa[i][j];
                 double t2 = lu[i][j];
                 double diff = Math.abs(t1 - t2);
                 assertTrue(diff < eps);
             }
         }
         
         double[] y = new double[]{8, 1.4, 1.5};

         // assert L*y = permutation * b
         double[] ly = MatrixUtil.multiplyMatrixByColumnVector(ell, y);
         double[] pb = MatrixUtil.multiplyMatrixByColumnVector(permutation, b);
         for (int i = 0; i < pb.length; ++i) {
             double t1 = ly[i];
             double t2 = pb[i];
             double diff = Math.abs(t1 - t2);
             assertTrue(diff < eps);
         }
         
         double[] expectedX = new double[]{-1.4, 2.2, 0.6};
         
         // assert U*x = y
         double[] ux = MatrixUtil.multiplyMatrixByColumnVector(u, expectedX);
         for (int i = 0; i < ux.length; ++i) {
             double t1 = ux[i];
             double t2 = y[i];
             double diff = Math.abs(t1 - t2);
             assertTrue(diff < eps);
         }
         
         double[] resultX = LinearEquations.LUPSolve(ell, u, p, b);
         
         for (int i = 0; i < resultX.length; ++i) {
             double t1 = resultX[i];
             double t2 = expectedX[i];
             double diff = Math.abs(t1 - t2);
             assertTrue(diff < eps);
         }
    }
    
    public void testLUDecompositionOfA() {

         /* test from Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms, chap. 28*/
         
         double[][] a = new double[4][];
         a[0] = new double[]{2, 3, 1, 5};
         a[1] = new double[]{6, 13, 5, 19};
         a[2] = new double[]{2, 19, 10, 23};
         a[3] = new double[]{4, 10, 11, 31};
         
         double[][] ell = new double[4][];
         ell[0] = new double[]{1, 0, 0, 0,};
         ell[1] = new double[]{3, 1, 0, 0};
         ell[2] = new double[]{1, 4, 1, 0};
         ell[3] = new double[]{2, 1, 7, 1};
         
         double[][] u = new double[4][];
         u[0] = new double[]{2, 3, 1, 5};
         u[1] = new double[]{0, 4, 2, 4};
         u[2] = new double[]{0, 0, 1, 2};
         u[3] = new double[]{0, 0, 0, 3};
         
         LU lu = LinearEquations.LUDecomposition(a);
         
         assertNotNull(lu);
         assertEquals(lu.ell.length, a.length);
         assertEquals(lu.ell[0].length, a[0].length);
         assertEquals(lu.u.length, a.length);
         assertEquals(lu.u[0].length, a[0].length);
         
         for (int i = 0; i < a.length; ++i) {
             assertTrue(Arrays.equals(lu.ell[i], ell[i]));
             assertTrue(Arrays.equals(lu.u[i], u[i]));
         }
    }
    
    public void testLUPDecompositionOfA() {

         /* test from Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms, chap. 28*/
         
         double[][] a = new double[4][];
         a[0] = new double[]{2, 0, 2, 0.6};
         a[1] = new double[]{3, 3, 4, -2};
         a[2] = new double[]{5, 5, 4, 2};
         a[3] = new double[]{-1, -2, 3.4, -1};
         
         double[][] ell = new double[4][];
         ell[0] = new double[]{1, 0, 0, 0,};
         ell[1] = new double[]{0.4, 1, 0, 0};
         ell[2] = new double[]{-0.2, 0.5, 1, 0};
         ell[3] = new double[]{0.6, 0, 0.4, 1};
         
         double[][] u = new double[4][];
         u[0] = new double[]{5, 5, 4, 2};
         u[1] = new double[]{0, -2, 0.4, -0.2};
         u[2] = new double[]{0, 0, 4, -0.5};
         u[3] = new double[]{0, 0, 0, -3};
         
         int[] p = new int[]{2, 0, 3, 1};
         
         LUP lup = LinearEquations.LUPDecomposition(a);
         
         assertNotNull(lup);
         assertEquals(lup.ell.length, a.length);
         assertEquals(lup.ell[0].length, a[0].length);
         assertEquals(lup.u.length, a.length);
         assertEquals(lup.u[0].length, a[0].length);
         
         double eps = 1e-7;
         double t1, t2;
         for (int i = 0; i < a.length; ++i) {
             for (int j = 0; j < a[i].length; ++j) {
                 t1 = lup.ell[i][j];
                 t2 = ell[i][j];
                 assertTrue(Math.abs(t1 - t2) < eps);
                 t1 = lup.u[i][j];
                 t2 = u[i][j];
                 assertTrue(Math.abs(t1 - t2) < eps);
             }
         }
         assertTrue(Arrays.equals(lup.p, p));
    }
    
    public void testLeastSquares() throws Exception {
        
        System.out.println("\ntestLeastSquares()");
        
        double[][] xy = new double[5][2];
        xy[0][0] = -1; xy[0][1] = 2;
        xy[1][0] =  1; xy[1][1] = 1;
        xy[2][0] =  2; xy[2][1] = 1;
        xy[3][0] =  3; xy[3][1] = 0;
        xy[4][0] =  5; xy[4][1] = 3;
        
        int polyOrder = 2;
        
        double eps = 0.01;
        
        boolean solveForFullRank = true;
        
        //double[][] tMatrix = MatrixUtil.calculateNormalizationMatrix2X3(xy);
        //xy = MatrixUtil.multiply(xy, tMatrix);
        
        double[] c = LinearEquations.leastSquaresPolynomial(xy, polyOrder, solveForFullRank);
        
        System.out.println(Arrays.toString(c));
        // if used normalization, c[0] should be divided by the scale tMatrix[0][0]
        
        double[] expected = new double[]{1.2, -0.757, 0.214};
        
        assertEquals(expected.length, c.length);
        
        for (int i = 0; i < expected.length; ++i) {
            double t1 = c[i];
            double t2 = expected[i];
            double diff = Math.abs(t1 - t2);
            assertTrue(diff < eps);
        }
        
        printLeastSquaresVariations(xy, solveForFullRank);
    }

    public void testLeastSquares3() throws Exception {
        System.out.println("\ntestLeastSquares3()");
         
        //(0, 0), (0, 0), (0, 12)
        double[][] xy = new double[3][2];
        xy[0][0] =  0; xy[0][1] = 0;
        xy[1][0] =  0; xy[1][1] = 0;
        xy[2][0] =  0; xy[1][1] = 12;
        
        boolean solveForFullRank = false;
        
        printLeastSquaresVariations(xy, solveForFullRank);
    }

    public void testLeastSquares4() throws Exception {
        System.out.println("\ntestLeastSquares4()");
        
        // from Strang's Introduction to Linear Algebra
        //    example 1, chap 4.3
        
        double[][] xy = new double[3][2];
        xy[0][0] = 0; xy[0][1] = 6;
        xy[1][0] =  1; xy[1][1] = 0;
        xy[2][0] =  2; xy[2][1] = 0;
        
        // B.T.W. can see that geometric median is (1, 0)
        
        int polyOrder = 1;
        
        double eps = 1e-17;
        
        // solves using pseudo-inverse = (inverse(A^T*A) * A^T)
        boolean solveForFullRank = true;
        
        //double[][] tMatrix = MatrixUtil.calculateNormalizationMatrix2X3(xy);
        //xy = MatrixUtil.multiply(xy, tMatrix);
        
        double[] c = LinearEquations.leastSquaresPolynomial(xy, polyOrder, solveForFullRank);
        
        // if used normalization, c[0] should be divided by the scale tMatrix[0][0]
              
        /*
        c0 + c1*x0 - err0 = y0
        c0 + c1*x1 - err1 = y1
        c0 + c1*x2 - err2 = y2
        A*X - Y = err
        */
       
        double[] expected = new double[]{5, -3};
        
        assertEquals(expected.length, c.length);
        
        for (int i = 0; i < expected.length; ++i) {
            double t1 = c[i];
            double t2 = expected[i];
            double diff = Math.abs(t1 - t2);
            assertTrue(diff < eps);
        }
         
        printLeastSquaresVariations(xy, solveForFullRank);
    }
    
    private void printLeastSquaresVariations(double[][] xy, boolean solveForFullRank) throws NotConvergedException {
        
        System.out.printf("(x,y)=");
        for (int i = 0; i < xy.length; ++i) {
            System.out.printf(" (%.2f,%.2f)", xy[i][0], xy[i][1]);
        }
        System.out.printf("\n");
        
        int[] polyOrders = new int[]{0, 1, 2};
        
        for (int polyOrder : polyOrders) {
        
            double[] c = LinearEquations.leastSquaresPolynomial(xy, polyOrder, solveForFullRank);
        
            // if used standard unit normalization, c[0] should be divided by the scale tMatrix[0][0]
            // create matrix A
            double[][] a = new double[xy.length][];
            double x;
            for (int i = 0; i < xy.length; ++i) {
                a[i] = new double[polyOrder + 1];
                a[i][0] = 1;
            }
            double[] y = new double[xy.length];
            for (int i = 0; i < xy.length; ++i) {
                x = xy[i][0];
                y[i] = xy[i][1];
                for (int k = 0; k < polyOrder; ++k) {
                    a[i][k + 1] = x * a[i][k];
                }
            }
            double[] ac = MatrixUtil.multiplyMatrixByColumnVector(a, c);
            double[] err = MatrixUtil.subtract(ac, y);
        
            System.out.printf("polyOrder=%d  c=%s  \n    err=%s\n", polyOrder,
                    Arrays.toString(c), Arrays.toString(err));
        }
        
        double[][] yx = new double[xy.length][];
        for (int i = 0; i < xy.length; ++i) {
            yx[i] = new double[]{xy[i][1], xy[i][0]};
        }

        for (int polyOrder : polyOrders) {

            double[] c = LinearEquations.leastSquaresPolynomial(yx, polyOrder, solveForFullRank);

            double[][] a = new double[xy.length][];
            for (int i = 0; i < xy.length; ++i) {
                a[i] = new double[polyOrder + 1];
                a[i][0] = 1;
            }
            double[] y_yx = new double[xy.length];
            double x_yx;
            for (int i = 0; i < xy.length; ++i) {
                x_yx = xy[i][1];
                y_yx[i] = xy[i][0];
                for (int k = 0; k < polyOrder; ++k) {
                    a[i][k + 1] = x_yx * a[i][k];
                }
            }
            double[] ac = MatrixUtil.multiplyMatrixByColumnVector(a, c);
            double[] err = MatrixUtil.subtract(ac, y_yx);

            System.out.printf("polyOrder=%d  c_yx=%s  \n    err=%s\n", polyOrder,
                    Arrays.toString(c), Arrays.toString(err));

        }
        
        System.out.flush();
    }
    
    public void testLDM() throws NotConvergedException {
        /*
        from Golub and van Loan, "Matrix Computations" Example 5.1-1
        */
        
        double[][] a = new double[3][];
        a[0] = new double[]{10, 10, 20};
        a[1] = new double[]{20, 25, 40};
        a[2] = new double[]{30, 50, 61};
        
        double[][] ellE = new double[3][];
        ellE[0] = new double[]{1, 0, 0};
        ellE[1] = new double[]{2, 1, 0};
        ellE[2] = new double[]{3, 4, 1};
        
        double[] dE = new double[]{10, 5, 1};
        
        double[][] mTE = new double[3][];
        mTE[0] = new double[]{1, 1, 2};
        mTE[1] = new double[]{0, 1, 0};
        mTE[2] = new double[]{0, 0, 1};
                
        LDM ldm = LinearEquations.LDMDecomposition(a);
        System.out.printf("A=\n%s\n", FormatArray.toString(a, "%.1f"));
        System.out.printf("L=\n%s\n", FormatArray.toString(ldm.ell, "%.1f"));
        System.out.printf("M^T=\n%s\n", FormatArray.toString(MatrixUtil.transpose(ldm.m), "%.1f"));
        System.out.printf("d=\n%s\n", FormatArray.toString(ldm.d, "%.1f"));

        double diff;
        double tol = 1e-7;
        int i, j;
        
        assertEquals(ellE.length, ldm.ell.length);
        assertEquals(ellE[0].length, ldm.ell[0].length);
        
        assertEquals(3, ldm.d.length);
        
        assertEquals(mTE.length, ldm.m[0].length);
        assertEquals(mTE[0].length, ldm.m.length);
     
        for (i = 0; i < ellE.length; ++i) {
            diff = Math.abs(dE[i] - ldm.d[i]);
            assertTrue(diff < tol);
            for (j = 0; j < ellE[i].length; ++j) {
                diff = Math.abs(ellE[i][j] - ldm.ell[i][j]);
                assertTrue(diff < tol);
                diff = Math.abs(mTE[i][j] - ldm.m[j][i]);
                assertTrue(diff < tol);
            }
        }
        
        double[][] aa = MatrixUtil.multiplyByDiagonal(ldm.ell, ldm.d);
        MatrixUtil.multiply(aa, MatrixUtil.transpose(ldm.m));
        
        System.out.printf("L*D*M^T=\n%s\n", FormatArray.toString(aa, "%.1f"));
        
    }
    
    public void testLDL() {
        /*
        from Golub and van Loan, "Matrix Computations" Example 5.1-2
        */
        
        double[][] a = new double[3][];
        a[0] = new double[]{10, 20, 30};
        a[1] = new double[]{20, 45, 80};
        a[2] = new double[]{30, 80, 171};
        
        double[][] ellE = new double[3][];
        ellE[0] = new double[]{1, 0, 0};
        ellE[1] = new double[]{2, 1, 0};
        ellE[2] = new double[]{3, 4, 1};
        
        double[] dE = new double[]{10, 5, 1};
        
        System.out.printf("A=\n%s\n", FormatArray.toString(a, "%.1f"));
        
        LDL ldl = LinearEquations.LDLDecomposition(a);
        System.out.printf("L=%s\n", FormatArray.toString(ldl.getL(), "%.1f"));
        System.out.printf("d=%s\n", FormatArray.toString(ldl.getD(), "%.1f"));
        System.out.printf("A=\n%s\n", FormatArray.toString(a, "%.1f"));
        
        double diff;
        double tol = 1e-7;
        int i, j;
        
        assertEquals(ellE.length, ldl.getL().length);
        assertEquals(ellE[0].length, ldl.getL()[0].length);
        
        assertEquals(3, ldl.getD().length);
        
        for (i = 0; i < ellE.length; ++i) {
            diff = Math.abs(dE[i] - ldl.getD()[i]);
            assertTrue(diff < tol);
            for (j = 0; j < ellE[i].length; ++j) {
                diff = Math.abs(ellE[i][j] - ldl.getL()[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        double[][] aa = MatrixUtil.multiplyByDiagonal(ldl.getL(), ldl.getD());
        MatrixUtil.multiply(aa, MatrixUtil.transpose(ldl.getL()));
        
        System.out.printf("L*D*L^T=\n%s\n", FormatArray.toString(aa, "%.1f"));
        /*
        for (i = 0; i < a.length; ++i) {
            for (j = 0; j < ellE[i].length; ++j) {
                diff = Math.abs(a[i][j] - aa[i][j]);
                assertTrue(diff < tol);
            }
        }*/
    }
    
    public void testCholeskyDecomposition() {
        // example 5.2-1 from Golub & Van Loan, "Matrix Computations"
        
        double eps = 1e-7;
        
        double[][] a  = new double[2][];
        a[0] = new double[]{2, -2};
        a[1] = new double[]{-2, 5};
        
        double[][] G = LinearEquations.choleskyDecompositionViaLDL(a, eps);
        
        double[][] expectedG = new double[2][];
        expectedG[0] = new double[]{Math.sqrt(2), 0};
        expectedG[1] = new double[]{-Math.sqrt(2), Math.sqrt(3)};
        
        double tol = 1e-3;
        double diff;
        int i, j;
        
        assertEquals(expectedG.length, G.length);
        assertEquals(expectedG[0].length, G[0].length);
        for (i = 0; i < G.length; ++i) {
            for (j = 0; j < G[0].length; ++j) {
                diff = Math.abs(G[i][j] - expectedG[i][j]);
                assertTrue(diff < tol);
            }
        }
        
    }

    public void testGaussianEliminationViaLU() throws NotConvergedException {
        double[][] a = new double[3][];
        a[0] = new double[]{1, 3, 0, 2, -1};
        // multiply a[0] by 2 to make sure the results are rref
        MatrixUtil.multiply(a[0], 2);
        a[1] = new double[]{0, 0, 1, 4, -3};
        a[2] = new double[]{1, 3, 1, 6, -4};

        double[][] expectedR = new double[3][];
        expectedR[0] = new double[]{1, 3, 0, 2, -1};
        expectedR[1] = new double[]{0, 0, 1, 4, -3};
        expectedR[2] = new double[]{0, 0, 0, 0, 0};

        double[][] r;

        r = LinearEquations.gaussianEliminationViaLU(a, 1E-7);

        assertEquals(r.length, expectedR.length);
        assertEquals(r[0].length, expectedR[0].length);
        int i;
        int j;
        for (i = 0; i < expectedR.length; ++i) {
            for (j = 0; j < expectedR[i].length; ++j) {
                assertTrue(Math.abs(expectedR[i][j] - r[i][j]) < 1E-7);
            }
        }

        // ========================
        a = new double[2][];
        a[0] = new double[]{1, 3, -1};
        // multiply a[0] by 2 to make sure the results are rref
        MatrixUtil.multiply(a[0], 2);
        a[1] = new double[]{0, 1, 7};

        expectedR = new double[2][];
        expectedR[0] = new double[]{1, 0, -22};
        expectedR[1] = new double[]{0, 1, 7};

        r = LinearEquations.gaussianEliminationViaLU(a, 1E-7);

        assertEquals(r.length, expectedR.length);
        assertEquals(r[0].length, expectedR[0].length);

        for (i = 0; i < expectedR.length; ++i) {
            for (j = 0; j < expectedR[i].length; ++j) {
                assertTrue(Math.abs(expectedR[i][j] - r[i][j]) < 1E-7);
            }
        }
    }
}
