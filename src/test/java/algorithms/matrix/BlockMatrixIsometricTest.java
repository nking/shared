package algorithms.matrix;

import algorithms.util.FormatArray;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BlockMatrixIsometricTest extends TestCase {
    
    public BlockMatrixIsometricTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
        int a0 = 4;
        int a1 = 9;
        int b0 = 2;
        int b1 = 3;
        double[][] a = new double[a0][];
        int i, j;
        for (i = 0; i < a0; ++i) {
            a[i] = new double[a1];
        }
        
        double v0;
        for (i = 0; i < a0; ++i) {
            v0 = Math.pow(10, i/b0);
            for (j = 0; j < a1; ++j) {
                a[i][j] = v0 + (j/b1);
            }
            System.out.printf("a[%d]=%s\n", i, FormatArray.toString(a[i], "%.0f"));
        }
        
        double[][] aOrig = MatrixUtil.copy(a);
        
        BlockMatrixIsometric mA = new BlockMatrixIsometric(a, b0, b1);

        double[][] _a = mA.getA();
        assertEquals(Arrays.hashCode(a), Arrays.hashCode(_a));
        
        double[][] b = new double[b0][b1];
        for (i = 0; i < b0; ++i) {
            Arrays.fill(b[i], 9);
        }
        
        mA.setBlock(b, 0, 1);
        
        double[][] _b = mA.getBlock(0, 1);
        
        double diff;
        double tol = 1.e-7;
        // assert a( block 0, block 1) is not equal to a0 (block 0, block 1)
        // assert _b == b
        int ii = 0;
        int jj;
        for (i = b0*0; i < b0*(0+1); ++i, ++ii) {
            jj = 0;
            for (j = b1*1; j < b1*(1+1); ++j, ++jj) {
                System.out.println("i=" + i + " j=" + j + " ii=" + ii + " jj=" + jj);
                diff = Math.abs(aOrig[i][j] - _b[ii][jj]);
                assertTrue(diff > tol);
                diff = Math.abs(a[i][j] - _b[ii][jj]);
                assertTrue(diff < tol);
                diff = Math.abs(b[ii][jj] - _b[ii][jj]);
                assertTrue(diff < tol);
            }
        }
        
        _b = MatrixUtil.zeros(b0, b1);
        mA.getBlock(_b, 0, 1);
        
        ii = 0;
        for (i = b0*0; i < b0*(0+1); ++i, ++ii) {
            jj = 0;
            for (j = b1*1; j < b1*(1+1); ++j, ++jj) {
                System.out.println("i=" + i + " j=" + j + " ii=" + ii + " jj=" + jj);
                diff = Math.abs(aOrig[i][j] - _b[ii][jj]);
                assertTrue(diff > tol);
                diff = Math.abs(a[i][j] - _b[ii][jj]);
                assertTrue(diff < tol);
                diff = Math.abs(b[ii][jj] - _b[ii][jj]);
                assertTrue(diff < tol);
            }
        }
        
        assertEquals(b0, mA.getBlockSize0());
                
        assertEquals(b1, mA.getBlockSize1());
        
    }
}
