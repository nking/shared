package algorithms.correlation;

import algorithms.matrix.MatrixUtil;
import java.util.logging.Logger;
import static org.junit.Assert.assertTrue;
import junit.framework.TestCase;
    
/**
 *
 * @author nichole
 */
public class BruteForceTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public BruteForceTest() {
    }
    
    public void testCovariance() {
        
        // unit test from 
        // https://jamesmccaffrey.wordpress.com/2017/11/03/example-of-calculating-a-covariance-matrix/

        double[][] a = new double[5][3];
        a[0] = new double[]{64.0,   580.0, 29.0};
        a[1] = new double[]{66.0,   570.0, 33.0};
        a[2] = new double[]{68.0,   590.0, 37.0};
        a[3] = new double[]{69.0,   660.0, 46.0};
        a[4] = new double[]{73.0,   600.0, 55.0};
        
        double[] expectedMean = new double[]{68.0,   600.0, 40.0};
        double[][] expectedCOV = new double[3][3];
        expectedCOV[0] = new double[]{11.50,    50.00,   34.75};
        expectedCOV[1] = new double[]{50.00,  1250.00,  205.00};
        expectedCOV[2] = new double[]{34.75,   205.00,  110.00};
        
        double[][] cov = BruteForce.covariance(a);

        double eps = 1.e-7;
        double diff;
        
        assertEquals(expectedCOV.length, cov.length);
        for (int i = 0; i < cov.length; ++i) {
            assertEquals(expectedCOV[i].length, cov[i].length);
            for (int j = 0; j < cov[i].length; ++j) {
                diff = Math.abs(expectedCOV[i][j] - cov[i][j]);
                assertTrue(diff < eps);
            };
        }
    }
    
    public void testCovariance2() {
        
        // unit test from 
        // https://jamesmccaffrey.wordpress.com/2017/11/03/example-of-calculating-a-covariance-matrix/

        int i, j;
        double eps = 1.e-7;
        double diff;
        
        double[][] a = new double[5][3];
        a[0] = new double[]{64.0,   580.0, 29.0};
        a[1] = new double[]{66.0,   570.0, 33.0};
        a[2] = new double[]{68.0,   590.0, 37.0};
        a[3] = new double[]{69.0,   660.0, 46.0};
        a[4] = new double[]{73.0,   600.0, 55.0};
        
        double[] mean = MatrixUtil.columnMeans(a);
        double[] stdev = MatrixUtil.standardDeviation(a);
        double[] expectedMean = new double[]{68.0,   600.0, 40.0};
        for (i = 0; i < a.length; ++i) {
            for (j = 0; j < a[i].length; ++j) {
                a[i][j] -= mean[j];
            }
        }
        
        //a = MatrixUtil.transpose(a);
        
        double[][] expectedCOV = new double[3][3];
        expectedCOV[0] = new double[]{11.50,    50.00,   34.75};
        expectedCOV[1] = new double[]{50.00,  1250.00,  205.00};
        expectedCOV[2] = new double[]{34.75,   205.00,  110.00};
        
        diff = Math.abs(expectedCOV[0][0] - (stdev[0]*stdev[0]));
        assertTrue(diff < eps);
        diff = Math.abs(expectedCOV[1][1] - (stdev[1]*stdev[1]));
        assertTrue(diff < eps);
        diff = Math.abs(expectedCOV[2][2] - (stdev[2]*stdev[2]));
        assertTrue(diff < eps);
                
        
        double[][] cov = BruteForce.covariance2(a);
        
        assertEquals(expectedCOV.length, cov.length);
        for (i = 0; i < cov.length; ++i) {
            assertEquals(expectedCOV[i].length, cov[i].length);
            for (j = 0; j < cov[i].length; ++j) {
                diff = Math.abs(expectedCOV[i][j] - cov[i][j]);
                assertTrue(diff < eps);
            };
        }
    }
    
}
