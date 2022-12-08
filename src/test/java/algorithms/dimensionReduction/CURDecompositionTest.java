package algorithms.dimensionReduction;

import algorithms.dimensionReduction.CURDecomposition.CUR;
import algorithms.dimensionReduction.CURDecomposition.PDFs;
import algorithms.dimensionReduction.CURDecomposition.SelectedFromA;
import algorithms.matrix.MatrixUtil;
import algorithms.matrix.MatrixUtil.SVDProducts;
import algorithms.util.FormatArray;
import java.util.Arrays;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVectorSub;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;
import no.uib.cipr.matrix.sparse.ArpackSym;

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
        
        System.out.printf("cdfs.rowsSelected=%s\n", Arrays.toString(cdfs.rowsSelected));
        System.out.printf("cdfs.pdfs.rowPDF=%s\n", FormatArray.toString(cdfs.pdfs.rowPDF, "%.4e"));
        
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
        //  in Figure 11.13... normalization possibly
        
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
        
        {
            
            System.out.printf("A dimensions are: %d x %d\n", a.length, a[0].length);
            System.out.printf("CUR dimensions are: %d x %d\n", result.length, result[0].length);
            
            double[][] aT = MatrixUtil.transpose(a);
            double[][] aTa = MatrixUtil.multiply(aT, a);
            double[][] aaT = MatrixUtil.multiply(a, aT);
            
            double[][] cT = MatrixUtil.transpose(result);
            double[][] cTc = MatrixUtil.multiply(cT, result);
            double[][] ccT = MatrixUtil.multiply(result, cT);
            
            EVD evd1aa = EVD.factorize(new DenseMatrix(aTa));
            EVD evd2aa = EVD.factorize(new DenseMatrix(aaT));
            EVD evd1cc = EVD.factorize(new DenseMatrix(cTc));
            EVD evd2cc = EVD.factorize(new DenseMatrix(ccT));
            
            SVD svd1aa = SVD.factorize(new DenseMatrix(aTa)); 
            SVD svd2aa = SVD.factorize(new DenseMatrix(aaT));
            
            SVD svd1cc = SVD.factorize(new DenseMatrix(cTc)); 
            SVD svd2cc = SVD.factorize(new DenseMatrix(ccT));
            
            SVD svda = SVD.factorize(new DenseMatrix(a));            
            SVD svdc = SVD.factorize(new DenseMatrix(result));
            
            System.out.printf("\nEVD(A^TA):\n  eigenvalues=\n    %s\n  leftEV=\n%s\n  rightEV=\n%s\n  real evs=\n%s\n",
                FormatArray.toString(evd1aa.getRealEigenvalues(), "%.4f"),
                evd1aa.getLeftEigenvectors().toString(),
                evd1aa.getRightEigenvectors().toString(),
                FormatArray.toString(evd1aa.getRealEigenvalues(), "%.4f")
            );
            System.out.printf("\nEVD(AA^T):\n  eigenvalues=\n    %s\n  leftEV=\n%s\n  rightEV=\n%s\n  real evs=\n%s\n",
                FormatArray.toString(evd2aa.getRealEigenvalues(), "%.4f"),
                evd2aa.getLeftEigenvectors().toString(),
                evd2aa.getRightEigenvectors().toString(),
                FormatArray.toString(evd2aa.getRealEigenvalues(), "%.4f")
            );
            System.out.printf("\nSVD(A^TA):\n  singularvalues=\n    %s\n  U=\n%s\n  V^T=\n%s\n  svs=\n%s\n",
                FormatArray.toString(svd1aa.getS(), "%.4f"),
                svd1aa.getU().toString(),
                svd1aa.getVt().toString(),
                FormatArray.toString(svd1aa.getS(), "%.4f")
            );
            System.out.printf("\nSVD(AA^T):\n  singularvalues=\n    %s\n  U=\n%s\n  V^T=\n%s\n  svs=\n%s\n",
                FormatArray.toString(svd2aa.getS(), "%.4f"),
                svd2aa.getU().toString(),
                svd2aa.getVt().toString(),
                FormatArray.toString(svd2aa.getS(), "%.4f")
            );
            System.out.printf("\nSVD(A):\n  singularvalues=\n    %s\n  U=\n%s\n  V^T=\n%s\n  svs=\n%s\n",
                FormatArray.toString(svda.getS(), "%.4f"),
                svda.getU().toString(),
                svda.getVt().toString(),
                FormatArray.toString(svda.getS(), "%.4f")
            );
            //----------------
            System.out.printf("\nEVD(C^TC):\n  eigenvalues=\n    %s\n  leftEV=\n%s\n  rightEV=\n%s\n  real evs=\n%s\n",
                FormatArray.toString(evd1cc.getRealEigenvalues(), "%.4f"),
                evd1cc.getLeftEigenvectors().toString(),
                evd1cc.getRightEigenvectors().toString(),
                FormatArray.toString(evd1cc.getRealEigenvalues(), "%.4f")
            );
            System.out.printf("\nEVD(CC^T):\n  eigenvalues=\n    %s\n  leftEV=\n%s\n  rightEV=\n%s\n  real evs=\n%s\n",
                FormatArray.toString(evd2cc.getRealEigenvalues(), "%.4f"),
                evd2cc.getLeftEigenvectors().toString(),
                evd2cc.getRightEigenvectors().toString(),
                FormatArray.toString(evd1cc.getRealEigenvalues(), "%.4f")
            );
            System.out.printf("\nSVD(C^TC):\n  singularvalues=\n    %s\n  U=\n%s\n  V^T=\n%s\n  svs=\n%s\n",
                FormatArray.toString(svd1cc.getS(), "%.4f"),
                svd1cc.getU().toString(),
                svd1cc.getVt().toString(),
                FormatArray.toString(svd1cc.getS(), "%.4f")
            );
            System.out.printf("\nSVD(CC^T):\n  singularvalues=\n    %s\n  U=\n%s\n  V^T=\n%s\n  svs=\n%s\n",
                FormatArray.toString(svd2cc.getS(), "%.4f"),
                svd2cc.getU().toString(),
                svd2cc.getVt().toString(),
                FormatArray.toString(svd2cc.getS(), "%.4f")
            );
            System.out.printf("\nSVD(C):\n  singularvalues=\n    %s\n  U=\n%s\n  V^T=\n%s\n  svs=\n%s\n",
                FormatArray.toString(svdc.getS(), "%.4f"),
                svdc.getU().toString(),
                svdc.getVt().toString(),
                FormatArray.toString(svdc.getS(), "%.4f")
            );
            System.out.printf("\nSVD(C).vT_0_0 = %.3f\n", svdc.getVt().get(0, 0));
            System.out.printf("\nSVD(C).vT_1_0 = %.3f\n", svdc.getVt().get(1, 0));
            
            /*
            sum of the r eigenvalues = sum of r diagonal elements of A (= trace of A).
            that is, each item i to rank r is the index of the eigenvalues above zero
            so only adding a_i_i for |lambda_i| > 0 and only adding |lambda_i| > 0.
            */
            System.out.printf("A^TA=\n%s\n", new DenseMatrix(aTa).toString());
            System.out.printf("AA^T=\n%s\n", new DenseMatrix(aaT).toString());
            System.out.printf("C^TC=\n%s\n", new DenseMatrix(cTc).toString());
            System.out.printf("CC^T=\n%s\n", new DenseMatrix(ccT).toString());
            System.out.printf("A=\n%s\n", new DenseMatrix(a).toString());
            System.out.printf("C=\n%s\n", new DenseMatrix(result).toString());
            
            // testing det(A - lambda*I) = 0
            double detATA = MatrixUtil.determinant(
                MatrixUtil.aMinusVectorTimesIdentity(aTa, svd1aa.getS()));
            double detCTC = MatrixUtil.determinant(
                MatrixUtil.aMinusVectorTimesIdentity(cTc, svd1cc.getS()));
            
            System.out.printf("det(A^TA - lambda*I) = %.4f\n", detATA);
            System.out.printf("det(C^TC - lambda*I) = %.4f\n", detCTC);
            
            System.out.printf("trace(A^TA) = %.4f\n", MatrixUtil.trace(aTa));
            System.out.printf("trace(eigeinvalues of A^TA) = %.4f\n", MatrixUtil.trace(svd1aa.getS()));
            
            System.out.printf("trace(C^TC) = %.4f\n", MatrixUtil.trace(cTc));
            System.out.printf("trace(eigeinvalues of C^TC) = %.4f\n", MatrixUtil.trace(svd1cc.getS()));
            
            //-------
            SVD svd1aT = SVD.factorize(new DenseMatrix(aT));
            System.out.printf("\nSVD(A^T):\n  singularvalues=\n    %s\n  U=\n%s\n  V^T=\n%s\n   svs=\n%s\n",
                FormatArray.toString(svd1aT.getS(), "%.4f"),
                svd1aT.getU().toString(),
                svd1aT.getVt().toString(),
                FormatArray.toString(svd1aT.getS(), "%.4f")
            );
            
            try {
            CUR curA = CURDecomposition.calculateDecomposition(a, k);
            SVDProducts svdCURA = curA.getApproximateSVD();
            
            System.out.printf("\nCUR(A):\n%s\n", curA.toString());
            System.out.printf("\nCUR.approxSVD(A):\n%s\n", svdCURA.toString());
            } catch (Throwable t) {}
            
            System.out.flush();
        }
    }
    
    /**
     * 
     @param a
     @return 
     */
    public static double trace(double[][] a, int n) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += a[i][i];
        }
        return sum;
    }
    
    
    public void est1() {
        
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
    
    private double find2ndSmallestEigenValue(Map<Double, DenseVectorSub> eig) {
        
        int n = eig.size();
        if (n < 2) {
            return -1;
        }
        
        double min1 = Double.MAX_VALUE;
        double min2 = Double.MAX_VALUE;
        
        for (Map.Entry<Double, DenseVectorSub> entry : eig.entrySet()) {
            double v = entry.getKey().doubleValue();
            if (min1 == Double.MAX_VALUE) {
                min1 = v;
            } else if (min2 == Double.MAX_VALUE) {
                if (v <= min1) {
                    min2 = min1;
                    min1 = v;
                } else {
                    min2 = v;
                }
            } else {
                if (v <= min1) {
                    min2 = min1;
                    min1 = v;
                } else if (v < min2) {
                    min2 = v;
                }
            }
        }
        
        return min2;
    }     
}
