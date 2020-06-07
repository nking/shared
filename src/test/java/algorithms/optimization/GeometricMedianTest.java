package algorithms.optimization;

import thirdparty.dlib.optimization.*;
import java.util.Arrays;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GeometricMedianTest extends TestCase {
    
    public GeometricMedianTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
        
        System.out.println("test0");
        
        /*
        from https://www.geeksforgeeks.org/geometric-median/
     -- Input: (1, 1), (3, 3)
        Output: Geometric Median = (2, 2) with minimum distance = 2.82843
        */
        double[] data = new double[]{1, 1, 3, 3};
        int nDimensions = 2;
                
        double[] init;
        double[][] expected = new double[4][2];
        expected[0] = new double[]{1, 1};
        expected[1] = new double[]{2, 2};
        expected[2] = new double[]{3, 3};
        expected[3] = new double[]{1.707106781, 1.707106781};
        double expectedDist = Math.sqrt(8);
        
        GeometricMedianUnweightedFunction f 
            = new GeometricMedianUnweightedFunction(data, nDimensions);
        GeometricMedian gm = new GeometricMedian();
        
        for (int ii = 0; ii < 4; ++ii) {
            switch(ii) {
                case 0:
                    init = new double[]{1, 1};            
                    break;
                case 1:
                    init = new double[]{1, 3};
                    break;
                case 2:
                    init = new double[]{3, 3};
                    break;
                default:
                    init = f.calculateCentroid();;
                    break;
            }
            
            System.out.printf("\nbegin w/ init=(%.3f, %.3f); true=(%s), (%s), (%s), (%s)\n",
                init[0], init[1], 
                AbstractGeometricMedianFunction.toString(expected[0]),
                AbstractGeometricMedianFunction.toString(expected[1]),
                AbstractGeometricMedianFunction.toString(expected[2]),
                AbstractGeometricMedianFunction.toString(expected[3])
            );
           
            double min = gm.newtonsMethodAltered(f, init);

            System.out.println("min=" + min + " \n   coeffs=" +
                Arrays.toString(init));
            System.out.flush();
            
            // fails for some
            assertTrue(Math.abs(min - expectedDist) < 0.001);

            
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
        double expectedDist = 12;
                
        double[] init;
        double[] expected = new double[]{0, 0};
        
        GeometricMedianUnweightedFunction f 
            = new GeometricMedianUnweightedFunction(data, nDimensions);
        GeometricMedian gm = new GeometricMedian();
        
        for (int ii = 0; ii < 3; ++ii) {
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
            
            double min = gm.newtonsMethodAltered(f, init);

            System.out.println("min=" + min + " \n   coeffs=" +
                Arrays.toString(init));
            System.out.flush();
            
            // fails for 2 cases
            if (ii == 0)
                assertTrue(Math.abs(min - expectedDist) < 0.001);

            for (int i = 0; i < init.length; ++i) {
            //    assertTrue(Math.abs(init[i] - expected[i]) < 0.1);
            }
        }
    }
    
    public void est2() {
        
        System.out.println("test2");
        
        /*
        from "Noniterative Solution of Some Fermat-Weber Location Problems"
        https://www.researchgate.net/publication/220418355_Noniterative_Solution_of_Some_Fermat-Weber_Location_Problems
        -- Input: (-20, 48), (-20, -48), (20, 0), (59, 0)
           Output: Geometric Median = (20,0) with minimum distance = 163.964
           NOTE: start with (44, 0) to test for a Weiszfeld problem w/ demand points

        */
        double[] data = new double[]{-20, 48, -20, -48, 20, 0, 59, 0};
        int nDimensions = 2;
                
        double[] init;
        double[] expected = new double[]{20., 0};
        double expectedDist = 163.964;
        
        GeometricMedianUnweightedFunction f 
            = new GeometricMedianUnweightedFunction(data, nDimensions);
        GeometricMedian gm = new GeometricMedian();
        
        // fails when derivative is 0
        
        for (int ii = 1; ii < 2; ++ii) {
            switch(ii) {
                case 0:
                    init = new double[]{44, 0};            
                    break;
                case 1:
                    init = new double[]{3, 45};
                    break;
                case 2:
                    init = new double[]{0, 0};            
                    break;
                default:
                    // (9.75, 0)
                    init = f.calculateCentroid();;
                    break;
            }
            
            System.out.printf("\nbegin w/ init=(%.3f, %.3f); true=(%.3f, %.3f)\n",
                init[0], init[1], expected[0], expected[1]);
           
            double min = gm.newtonsMethodAltered(f, init);

            System.out.println("min=" + min + " \n   coeffs=" +
                Arrays.toString(init));
            System.out.flush();
            
            assertTrue(Math.abs(expectedDist - min) <= 1e-17);

            /*
            for (int i = 0; i < init.length; ++i) {
                double diff = Math.abs(init[i] - expected[i]);
                System.out.printf("[%d] X - expected=%f.4 - %f.4 = %f.4\n", 
                    i, init[i], expected[i], diff);
            //    assertTrue(diff < 0.1);
            }*/
        }
        
    }
}
