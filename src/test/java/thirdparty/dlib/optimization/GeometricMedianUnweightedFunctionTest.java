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
    
    public void test0() {
        
        
        System.out.println("test0");
        
        /*
        from https://www.geeksforgeeks.org/geometric-median/
     -- Input: (1, 1), (3, 3)
        Output: Geometric Median = (1,1) or (2, 2) or (3,3) with minimum distance = 2.82843
        */
        double[] data = new double[]{1, 1, 3, 3};
        int nDimensions = 2;
                
        double[] init;
        double[][] expected = new double[3][2];
        expected[0] = new double[]{1, 1};
        expected[1] = new double[]{2, 2};
        expected[2] = new double[]{3, 3};
        double expectedDist = Math.sqrt(8);
        
        GeometricMedianUnweightedFunction f 
            = new GeometricMedianUnweightedFunction(data, nDimensions);
        
        for (int ii = 0; ii < 3; ++ii) {
            switch(ii) {
                case 0:
                    init = new double[]{1, 1};            
                    break;
                case 1:
                    init = new double[]{1, 3};
                    break;
                default:
                    init = f.calculateCentroid();;
                    break;
            }
            
            System.out.printf("\nbegin w/ init=(%.3f, %.3f); true=(%s), (%s), (%s)\n",
                init[0], init[1], 
                AbstractGeometricMedianFunction.toString(expected[0]),
                AbstractGeometricMedianFunction.toString(expected[1]),
                AbstractGeometricMedianFunction.toString(expected[2]));
           
            LBFGSSearchStrategy searchStrategy = new LBFGSSearchStrategy(5);
            ObjectiveDeltaStopStrategy stopStrategy 
                = new ObjectiveDeltaStopStrategy(1e-5, 
                100);

            double acceptableMinF = -1;

            LBFGSOptimization opt = new LBFGSOptimization();
            double min = opt.findMin(searchStrategy, stopStrategy, f, init, 
                acceptableMinF);

            System.out.println("min=" + min + " \n   coeffs=" +
                Arrays.toString(init));
            System.out.flush();
            
            assertTrue(Math.abs(min - expectedDist) < 0.001);

            /*
            boolean matched = false;
            int c = 0;
            for (int soln = 0; soln < expected.length; ++soln) {
                c = 0;
                double[] s = expected[soln];    
                for (int i = 0; i < init.length; ++i) {
                    double diff = Math.abs(init[i] - s[i]);
                    if (diff < 0.001) {
                        c++;
                    }
                }
                if (c == init.length) {
                    matched = true;
                    break;
                }
            }
            assertTrue(matched);
            */
        }
        
    }
    
    public void est1() {
        
        System.out.println("test1");
        
        /*
        from https://www.geeksforgeeks.org/geometric-median/
     -- Input: (0, 0), (0, 0), (0, 12)
        Output: Geometric Median = (0, 0) with minimum distance = 12
        */
        double[] data = new double[]{0, 0,  0, 0,  0, 12};
        int nDimensions = 2;
                
        double[] init;
        double[] expected = new double[]{0, 0};
        GeometricMedianUnweightedFunction f 
            = new GeometricMedianUnweightedFunction(data, nDimensions);
        
        for (int ii = 0; ii < 2; ++ii) {
            switch(ii) {
                case 0:
                    init = new double[]{0, 0};            
                    break;
                case 1:
                    init = new double[]{0, 12};
                    break;
                default:
                    init = f.calculateCentroid();;
                    break;
            }
            
            System.out.printf("\nbegin w/ init=(%.3f, %.3f); true=(%.3f, %.3f)\n",
                init[0], init[1], expected[0], expected[1]);
            
            LBFGSSearchStrategy searchStrategy = new LBFGSSearchStrategy(5);
            ObjectiveDeltaStopStrategy stopStrategy 
                = new ObjectiveDeltaStopStrategy(AbstractGeometricMedianFunction.eps, 
                100);

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
}
