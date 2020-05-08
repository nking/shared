package algorithms.matrix;

import java.util.Arrays;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;
    
/**
 *
 * @author nichole
 */
public class MiscTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public MiscTest() {
    }
    
    public void testMultiply() throws Exception {

        double[][] a = new double[3][2];
        a[0] = new double[]{1, 2};
        a[1] = new double[]{2, 3};
        a[2] = new double[]{3, 4};

        double[] b = new double[]{4, 3};

        double[] m = Misc.multiply(a, b);

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
        
        double[][] d = Misc.multiply(a, c);

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

        float[][] cc = Misc.transpose(bb);
                for (int i = 0; i < cc.length; i++) {
            for (int j = 0; j < cc[i].length; j++) {
                assertTrue(expected[i][j] == cc[i][j]);
            }
        }

        float[][] dd = Misc.transpose(cc);

        for (int i = 0; i < dd.length; i++) {
            for (int j = 0; j < dd[i].length; j++) {
                assertTrue(bb[i][j] == dd[i][j]);
            }
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
        
        double[][] inv = Misc.pseudoinverse(a);
        
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
        double[][] inv = Misc.pseudoinverse2(a);
        
        double[] y = new double[]{2, 1, 1, 0, 3};
        
        double[] c = Misc.multiply(inv, y);
        
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
        double[][] inv2 = Misc.pseudoinverse(a);
                
    }
    
    public void testMultiply2() throws Exception {

        double[][] m1 = new double[2][3];
        m1[0] = new double[]{0, 1, 0};
        m1[1] = new double[]{1000, 100, 10};

        double[][] m2 = new double[3][2];
        m2[0] = new double[]{2, 1};
        m2[1] = new double[]{3, 0};
        m2[2] = new double[]{4, 0};
                
        double[][] m = Misc.multiply(m1, m2);

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

        double[] c = Misc.subtract(a, b);

        assertTrue(Arrays.equals(expected, c));
    }
}
