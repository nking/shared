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
public class DistanceTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public DistanceTest() {
    }
    
    public void testSort() {
        
        //double[][] _sortCheck(double[] x, double[] y)
        
        double eps = 1.e-17;
        
        double[] x = new double[]{1., 1.1, 1.1, 2., 3.};
        double[] y = new double[]{1., 1.1, 2.1, 3., 4.};
        
        double[] eX = Arrays.copyOf(x, x.length);
        double[] eY = Arrays.copyOf(y, y.length);
        
        double[][] sorted = Distance._sortCheck(x, y);
        double diff;
        for (int i = 0; i < x.length; ++i) {
            diff = Math.abs(sorted[0][i] - eX[i]);
            assertTrue(diff < eps);
            diff = Math.abs(sorted[1][i] - eY[i]);
            assertTrue(diff < eps);
        }
        
        x = new double[]{2, 1., 1.1, 3, 1.1};
        y = new double[]{3, 1., 2.1, 4, 1.1};
        sorted = Distance._sortCheck(x, y);
        for (int i = 0; i < x.length; ++i) {
            diff = Math.abs(sorted[0][i] - eX[i]);
            assertTrue(diff < eps);
            diff = Math.abs(sorted[1][i] - eY[i]);
            assertTrue(diff < eps);
        }
    }
}
