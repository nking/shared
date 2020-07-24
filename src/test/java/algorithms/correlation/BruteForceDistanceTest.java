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
public class BruteForceDistanceTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public BruteForceDistanceTest() {
    }
  
    public void testCorrelation0() {
        
        System.out.println("testCorrelation0");
        
        // unit tests from https://blogs.sas.com/content/iml/2018/04/04/distance-correlation.html
        /*
        https://blogs.sas.com/content/iml/2018/04/04/distance-correlation.html
        x = {1,0,-1, 0};
        y = {0,1, 0,-1};
        results = DistCorr(x, y);
        dCov=0.25
        dCor=0.4472136
        dVarX=0.55917
        dVarY=0.55917
        */
        
        double[][] a, b;
        double dCor, eDCor;
        int i, j;
        
        a = new double[4][1];
        b = new double[4][1];
        a[0] = new double[]{1}; a[1] = new double[]{0}; a[2] = new double[]{-1};a[3] = new double[]{0};
        b[0] = new double[]{0}; b[1] = new double[]{1}; b[2] = new double[]{0};b[3] = new double[]{-1};
        
        dCor = BruteForceDistance.correlation1(a, b);

        System.out.println("dCor=" + dCor);
        
    }
    
    public void testCorrelation1() {
        
        System.out.println("testCorrelation1");
        
        // unit tests from https://blogs.sas.com/content/iml/2018/04/04/distance-correlation.html
        /*

x = do(-1, 1, 0.1)`; 
y = x##2; 
dcor=0.4913789
        
        */
        double[][] a, b;
        double dCor, eDCor;
        int i, j;
        
        a = new double[4][1];
        b = new double[4][1];
        a[0] = new double[]{-1}; a[1] = new double[]{1}; a[2] = new double[]{0};a[3] = new double[]{1};
        b[0] = new double[]{1}; b[1] = new double[]{1}; b[2] = new double[]{0};b[3] = new double[]{1};
        
        dCor = BruteForceDistance.correlation1(a, b);

        System.out.println("dCor=" + dCor);
        
    }
}
