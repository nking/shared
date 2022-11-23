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
            //System.out.printf("a[%d]=%s%n", i, FormatArray.toString(a[i], "%.0f"));
        }
        
        double[][] aOrig = MatrixUtil.copy(a);
        
        BlockMatrixIsometric mA = new BlockMatrixIsometric(a, b0, b1);

        double[][] _a = mA.getA();
        assertTrue(Arrays.equals(a, _a));
        
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
                //System.out.println("i=" + i + " j=" + j + " ii=" + ii + " jj=" + jj);
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
                //System.out.println("i=" + i + " j=" + j + " ii=" + ii + " jj=" + jj);
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
        
        /*
        a[0]=1, 1, 1, 9, 9, 9, 3, 3, 3 
        a[1]=1, 1, 1, 9, 9, 9, 3, 3, 3 
        a[2]=10, 10, 10, 11, 11, 11, 12, 12, 12 
        a[3]=10, 10, 10, 11, 11, 11, 12, 12, 12
        int a0 = 4;
        int a1 = 9;
        int b0 = 2;
        int b1 = 3;
        */
        double[][] before = mA.getBlock(1, 1);
        double[][] c = new double[b0][];
        c[0] = new double[]{1, 2, 1};
        c[1] = new double[]{3, 4, 5};
        mA.subtractFromBlock(c, 1, 1);
        double[][] ec = new double[b0][];
        ec[0] = new double[]{10, 9, 10};
        ec[1] = new double[]{8, 7, 6};
        
        _b = MatrixUtil.zeros(b0, b1);
        mA.getBlock(_b, 1, 1);
        
        for (i = 0; i < b0; ++i) {
            for (j = 0; j < b1; ++j) {
                diff = Math.abs(ec[i][j] - _b[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        mA.addToBlock(c, 1, 1);
        _b = MatrixUtil.zeros(b0, b1);
        mA.getBlock(_b, 1, 1);
        for (i = 0; i < b0; ++i) {
            for (j = 0; j < b1; ++j) {
                diff = Math.abs(before[i][j] - _b[i][j]);
                assertTrue(diff < tol);
            }
        }
        
    }
    
    public void testCopySet() {
        /*
        a = | 0  1  2  3  |
            | 4  5  6  7  |
        b = | 0  10 20 30 |
            | 40 50 60 70 |
        blockSize0 = 1, blockSize2 = 2
        */
        double[][] a = new double[2][];
        a[0] = new double[]{0, 1, 2, 3};
        a[1] = new double[]{4, 5, 6, 7};
        double[][] b = new double[2][];
        b[0] = new double[]{0, 10, 20, 30};
        b[1] = new double[]{40, 50, 60, 70};
        
        BlockMatrixIsometric ba = new BlockMatrixIsometric(a, 1, 2);
        BlockMatrixIsometric bb = new BlockMatrixIsometric(b, 1, 2);
        
        BlockMatrixIsometric ba2 = ba.copy();
        assertEquals(ba.getBlockSize0(), ba2.getBlockSize0());
        assertEquals(ba.getBlockSize1(), ba2.getBlockSize1());
        assertEquals(ba.getA().length, ba2.getA().length);
        assertEquals(ba.getA()[0].length, ba2.getA()[0].length);
        
        int i, j;
        double diff;
        double tol = 1e-11;
        for (i = 0; i < a.length; ++i) {
            for (j = 0; j < a[i].length; ++j) {
                diff = Math.abs(ba.getA()[i][j] - ba2.getA()[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        ba.set(bb);
        assertEquals(bb.getBlockSize0(), ba.getBlockSize0());
        assertEquals(bb.getBlockSize1(), ba.getBlockSize1());
        assertEquals(bb.getA().length, ba.getA().length);
        assertEquals(bb.getA()[0].length, ba.getA()[0].length);
        for (i = 0; i < a.length; ++i) {
            for (j = 0; j < a[i].length; ++j) {
                diff = Math.abs(ba.getA()[i][j] - bb.getA()[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        ba.reset();
        for (i = 0; i < a.length; ++i) {
            for (j = 0; j < a[i].length; ++j) {
                diff = Math.abs(ba.getA()[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        //test setColumnBlock and setRowBlock
        a = new double[2][];
        a[0] = new double[]{0, 1, 2, 3};
        a[1] = new double[]{4, 5, 6, 7};
        ba = new BlockMatrixIsometric(a, 1, 2);
        
        double[] tmp = new double[2];
        double[] row = new double[]{100, 200};
        
        ba.getRowBlock(tmp, 1, 1);
        assertTrue(Math.abs(tmp[0] - 6.) < 1e-15);
        assertTrue(Math.abs(tmp[1] - 7.) < 1e-15);
        
        ba.setRowBlock(row, 1, 1);
        ba.getRowBlock(tmp, 1, 1);
        assertTrue(Math.abs(tmp[0] - 100.) < 1e-15);
        assertTrue(Math.abs(tmp[1] - 200.) < 1e-15);
        
        tmp[0] = -1; tmp[1] = 2;
        ba.addToRowBlock(tmp, 0, 1);
        ba.getRowBlock(tmp, 0, 1);
        assertTrue(Math.abs(tmp[0] - (2-1.)) < 1e-15);
        assertTrue(Math.abs(tmp[1] - (3+2.)) < 1e-15);
        
        tmp[0] = -1; tmp[1] = 2;
        ba.subtractFromRowBlock(tmp, 1, 0);
        ba.getRowBlock(tmp, 1, 0);
        assertTrue(Math.abs(tmp[0] - (4 - -1.)) < 1e-15);
        assertTrue(Math.abs(tmp[1] - (5 - 2.)) < 1e-15);
        
        // test column methods
        a = new double[4][];
        a[0] = new double[]{0, 1};
        a[1] = new double[]{4, 5};
        a[2] = new double[]{2, 3};
        a[3] = new double[]{6, 7};
        ba = new BlockMatrixIsometric(a, 2, 1); 
        
        double[] col = new double[]{100, 200};
        
        ba.getColumnBlock(tmp, 0, 1);
        assertTrue(Math.abs(tmp[0] - 1.) < 1e-15);
        assertTrue(Math.abs(tmp[1] - 5.) < 1e-15);
        
        ba.setColumnBlock(col, 0, 1);
        ba.getColumnBlock(tmp, 0, 1);
        assertTrue(Math.abs(tmp[0] - 100.) < 1e-15);
        assertTrue(Math.abs(tmp[1] - 200.) < 1e-15);
        
        tmp[0] = -1; tmp[1] = 2;
        ba.addToColumnBlock(tmp, 1, 0);
        ba.getColumnBlock(tmp, 1, 0);
        assertTrue(Math.abs(tmp[0] - (2-1.)) < 1e-15);
        assertTrue(Math.abs(tmp[1] - (6+2.)) < 1e-15);
        
        tmp[0] = -1; tmp[1] = 2;
        ba.subtractFromColumnBlock(tmp, 1, 1);
        ba.getColumnBlock(tmp, 1, 1);
        assertTrue(Math.abs(tmp[0] - (3 - -1.)) < 1e-15);
        assertTrue(Math.abs(tmp[1] - (7 - 2.)) < 1e-15);
    }
    
}
