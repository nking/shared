package algorithms.correlation;

import algorithms.matrix.*;
import java.util.Arrays;
import java.util.Iterator;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.sparse.FlexCompColMatrix;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
    
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
    
}
