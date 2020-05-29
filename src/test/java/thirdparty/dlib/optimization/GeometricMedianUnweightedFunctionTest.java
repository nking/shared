package thirdparty.dlib.optimization;

import java.util.Arrays;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GeometricMedianUnweightedFunctionTest extends TestCase {
    
    public GeometricMedianUnweightedFunctionTest(String testName) {
        super(testName);
    }
    
    public void est0() {
        /*
        from https://www.geeksforgeeks.org/geometric-median/
     -- Input: (1, 1), (3, 3)
        Output: Geometric Median = (2, 2) with minimum distance = 2.82843
        */
        double[] data = new double[]{1, 1, 3, 3};
        int nDimensions = 2;
        
        GeometricMedianUnweightedFunction f 
            = new GeometricMedianUnweightedFunction(data, nDimensions);
        
        double[] init = f.calculateCentroid();
        
        double[] expected = new double[]{2, 2};
    }
    
    public void test1() {
        
        /*
        from https://www.geeksforgeeks.org/geometric-median/
     -- Input: (0, 0), (0, 0), (0, 12)
        Output: Geometric Median = (0, 0) with minimum distance = 12
        */
        double[] data = new double[]{0, 0,  0, 0,  0, 12};
        int nDimensions = 2;
        
        GeometricMedianUnweightedFunction f 
            = new GeometricMedianUnweightedFunction(data, nDimensions);
        
        // gets stuck around centroid
        double[] init = 
            new double[]{0, 0};
            //new double[]{0, 12};
            //f.calculateCentroid();
        
        double[] expected = new double[]{0, 0};
        
        LBFGSSearchStrategy searchStrategy = new LBFGSSearchStrategy(5);
        ObjectiveDeltaStopStrategy stopStrategy 
            = new ObjectiveDeltaStopStrategy(AbstractGeometricMedianFunction.eps, 100);
        
        double acceptableMinF = -12;
        
        LBFGSOptimization opt = new LBFGSOptimization();
        double min = opt.findMin(searchStrategy, stopStrategy, f, init, 
            acceptableMinF);
       
        System.out.println("min=" + min + " \n   coeffs=" +
            Arrays.toString(init));
        System.out.flush();
       
        for (int i = 0; i < init.length; ++i) {
            assertTrue(Math.abs(init[i] - expected[i]) < 0.1);
        }
    }
}
