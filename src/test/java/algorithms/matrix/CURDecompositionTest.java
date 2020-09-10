package algorithms.matrix;

import algorithms.matrix.CURDecomposition.CUR;
import algorithms.matrix.CURDecomposition.PDFs;
import algorithms.matrix.CURDecomposition.SelectedFromA;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * @author nichole
 */
public class CURDecompositionTest extends TestCase {
    
    public CURDecompositionTest(String testName) {
        super(testName);
    }
    
    public void test0() throws NotConvergedException {
        
        // from Chap 11 of book "Mining of Massive 
        // Datasets" by Jure Leskovec, Anand Rajaraman, Jeff Ullman
        // http://www.mmds.org/
        
        double[][] a = new double[7][5];
        a[0] = new double[]{1, 1, 1, 0, 0};
        a[1] = new double[]{3, 3, 3, 0, 0};
        a[2] = new double[]{4, 4, 4, 0, 0};
        a[3] = new double[]{5, 5, 5, 0, 0};
        a[4] = new double[]{0, 0, 0, 4, 4};
        a[5] = new double[]{0, 0, 0, 5, 5};
        a[6] = new double[]{0, 0, 0, 2, 2};
        
        int k = 2;
        
        double[] expectedRPDF = new double[]{0.012, 0.111, 0.198, 0.309, 0.132, 0.206, 0.033};
        double[] expectedCPDF = new double[]{0.21, 0.21, 0.21, 0.185, 0.185};
        
        PDFs pdfs = CURDecomposition._calculatePDFs(a);
        
        double tol= 0.01;
        double diff;
        for (int i = 0; i < a.length; ++i) {
            diff = expectedRPDF[i] - pdfs.rowPDF[i];
            assertTrue(Math.abs(diff) < tol);
        }
        
        for (int i = 0; i < a[0].length; ++i) {
            diff = expectedCPDF[i] - pdfs.colPDF[i];
            assertTrue(Math.abs(diff) < tol);
        }
        
        CURDecomposition.CDFs cdfs = CURDecomposition._calculateCDFs(a, k);
        
        // --- change the selected columns and rows to match the text example:        
        
        // test col select 2 and 4:
        //The scaled columns are then  [1.54, 4.63, 6.17, 7.72, 0, 0, 0]^T 
        //  and [0, 0, 0, 0, 6.58, 8.22, 3.29]^T
        cdfs.colsSelected = new int[]{1, 3};
        double[][] expectedC = new double[2][a.length];
        expectedC[0] = new double[]{1.54, 4.63, 6.17, 7.72, 0, 0, 0};
        expectedC[1] = new double[]{0, 0, 0, 0, 6.58, 8.22, 3.29};
        expectedC = MatrixUtil.transpose(expectedC);
        
        // test row select 5 and 3:
        //a[3] = new double[]{5, 5, 5, 0, 0};
        //a[5] = new double[]{0, 0, 0, 5, 5};
        //   [ 0 0 0 7.79 7.79] and [6.36 6.36 6.36 0 0 ]
        cdfs.rowsSelected = new int[]{5, 3};
        double[][] expectedR = new double[2][a[0].length];
        expectedR[0] = new double[]{0, 0, 0, 7.79, 7.79};
        expectedR[1] = new double[]{6.36, 6.36, 6.36, 0, 0};
        
        
        SelectedFromA r = CURDecomposition._calculateR(a, cdfs.rowsSelected, cdfs.pdfs.rowPDF);
        assertEquals(expectedR.length, r.r.length);
        assertEquals(expectedR[0].length, r.r[0].length);
        for (int i = 0; i < r.r.length; ++i) {
            for (int j = 0; j < r.r[i].length; ++j) {
                diff = expectedR[i][j] - r.r[i][j];
                assertTrue(Math.abs(diff) < tol);
            }
        }
        
        SelectedFromA c = CURDecomposition._calculateR(MatrixUtil.transpose(a), 
            cdfs.colsSelected, cdfs.pdfs.colPDF);
        c.r = MatrixUtil.transpose(c.r);
        assertEquals(expectedC.length, c.r.length);
        assertEquals(expectedC[0].length, c.r[0].length);
        for (int i = 0; i < c.r.length; ++i) {
            for (int j = 0; j < c.r[i].length; ++j) {
                diff = expectedC[i][j] - c.r[i][j];
                assertTrue(Math.abs(diff) < tol);
            }
        }
        
        double[][] expectedU = new double[2][2];
        expectedU[0] = new double[]{0, 1./25.};
        expectedU[1] = new double[]{1./25., 0.};
        
        double[][] u = CURDecomposition._calculateU(a, r.indexesUnique, 
            c.indexesUnique);
        assertEquals(expectedU.length, u.length);
        assertEquals(expectedU[0].length, u[0].length);
        for (int i = 0; i < u.length; ++i) {
            for (int j = 0; j < u[i].length; ++j) {
                diff = expectedU[i][j] - u[i][j];
                assertTrue(Math.abs(diff) < tol);
            }
        }
        
        //NOTE: for some reason, the authors multiply all of matrix R by sqrt(2)
        //  in Figure 11.13.
        
        double[][] expectedCUR = new double[7][5];
        expectedCUR[0] = new double[]{0.3929077125588594, 0.3929077125588594,
            0.3929077125588594, 0, 0
        };
        expectedCUR[1] = new double[]{1.1787231376765783, 1.1787231376765783,
            1.1787231376765783, 0, 0
        };
        expectedCUR[2] = new double[]{1.5716308502354377, 1.5716308502354377,
            1.5716308502354377, 0, 0
        };
        expectedCUR[3] = new double[]{1.9645385627942973, 1.9645385627942973,
            1.9645385627942973, 0, 0
        };
        expectedCUR[4] = new double[]{
            0, 0, 0, 2.04915592378911, 2.04915592378911
        };
        expectedCUR[5] = new double[]{
            0, 0, 0, 2.5614449047363874, 2.5614449047363874
        };
        expectedCUR[6] = new double[]{
            0, 0, 0, 1.024577961894555, 1.024577961894555
        };
        
        double[][] result = MatrixUtil.multiply(c.r, u);
        result = MatrixUtil.multiply(result, r.r);
        
        assertEquals(expectedCUR.length, result.length);
        assertEquals(expectedCUR[0].length, result[0].length);
        for (int i = 0; i < result.length; ++i) {
            for (int j = 0; j < result[i].length; ++j) {
                diff = expectedCUR[i][j] - result[i][j];
                assertTrue(Math.abs(diff) < tol);
            }
        }
    }
    
    public void test1() {
        
        // from Chap 11 of book "Mining of Massive 
        // Datasets" by Jure Leskovec, Anand Rajaraman, Jeff Ullman
        // http://www.mmds.org/
        
        double[][] a = new double[7][5];
        a[0] = new double[]{1, 1, 1, 0, 0};
        a[1] = new double[]{3, 3, 3, 0, 0};
        a[2] = new double[]{4, 4, 4, 0, 0};
        a[3] = new double[]{5, 5, 5, 0, 0};
        a[4] = new double[]{0, 0, 0, 4, 4};
        a[5] = new double[]{0, 0, 0, 5, 5};
        a[6] = new double[]{0, 0, 0, 2, 2};
        
        int k = 2;
        
        CUR cur;
        
        int nTests = 100;
        for (int i = 0; i < nTests; ++i) {
            try {
                cur = CURDecomposition.calculateDecomposition(a, k);
                
                //TODO: check that ||CUR-A||_frob_norm <= (1+eps)*||A-A_k||_frob_norm
                //   with probability 98%
                
            } catch (Throwable ex) {
                Logger.getLogger(CURDecompositionTest.class.getName())
                    .log(Level.SEVERE, null, ex);
            }
        }
        
    }
}
