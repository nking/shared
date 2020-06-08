package algorithms.optimization;

import algorithms.misc.Standardization;
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
        /*
        standardized: Mean = new double[]{2, 2};   -1 -1  1 1
        standardized: StDev = new double[]{1.414, 1.414};
        standardized: new double[]{-0.707, -0.707, 0.707, 0.707};
        */
        
        double[] data0 = new double[]{1, 1, 3, 3};
        int nDimensions = 2;
        double expectedDist = Math.sqrt(8);  // in standaridized=2
        
        double tol = 0.01;
                
        double[] init;
        // the geometric-median in original data frames:
        double[] expected = new double[]{2, 2};
        // the geometric-median in standardized coordinate frames:
        double[] standardizedExpected = new double[]{0., 0};
        
        double[] standardizedMean = new double[nDimensions];
        double[] standardizedStDev = new double[nDimensions];
        double[] data = Standardization.standardUnitNormalization(data0, 
            nDimensions, standardizedMean, standardizedStDev);
        
        GeometricMedianUnweightedFunction f0 
            = new GeometricMedianUnweightedFunction(data0, nDimensions);
        GeometricMedianUnweightedFunction f 
            = new GeometricMedianUnweightedFunction(data, nDimensions);
        GeometricMedian gm = new GeometricMedian();
        
        for (int ii = 0; ii < 3; ++ii) {
            switch(ii) {
                case 0:
                    init = new double[]{0, 0};            
                    break;
                case 1:
                    //init = new double[]{0, 12};
                    init = new double[]{0, 1.1547};
                    break;
                default:
                    init = f.calculateCentroid();;
                    break;
            }
            
            System.out.printf("\nbegin w/ init=(%.3f, %.3f); true=(%.3f, %.3f)\n",
                init[0], init[1], 
                standardizedExpected[0], standardizedExpected[1]);
                //expected[0], expected[1]);
            
            double min = gm.newtonsMethod2(f, init);

            System.out.println("in standardized units: min=" + min + " \n   coeffs=" +
                Arrays.toString(init));
            System.out.flush();
            
            // de-normalize the geometric-median "init" and recalculate min
            //    in the natural coordinates:
            init = Standardization.standardUnitDenormalization(init, nDimensions, 
                standardizedMean, standardizedStDev);
        
            min = f0.f(init);
            
            System.out.println("in data units: min=" + min + " \n   coeffs=" +
                Arrays.toString(init));
            System.out.flush();
            
            assertTrue(Math.abs(min - expectedDist) <= tol);

            for (int i = 0; i < init.length; ++i) {
                assertTrue(Math.abs(init[i] - expected[i]) < tol);
            }
        }
        
    }
    
    public void est1() {
        
        System.out.println("test1");
        
        /*
        from https://www.geeksforgeeks.org/geometric-median/
     -- Input: (0, 0), (0, 0), (0, 12)
        Output: Geometric Median = (0, 0) with minimum distance = 12
        */
        /*
        standardized: Mean = new double[]{0, 4};
        standardized: StDev = new double[]{0, 6.928};
        standardized: new double[]{0, -0.5774, 0, -0.577, 0, 1.1547};
        */
        double[] data0 = new double[]{0, 0,  0, 0,  0, 12};
        int nDimensions = 2;
        double expectedDist = 12;
        
        double tol = 0.01;
                
        double[] init;
        double[] expected = new double[]{0, 0};
        double[] standardizedExpected = new double[]{0., -0.5773};
        // expected in standardized coords = (0-0)/0  , (0-4)/6.9282 = (0, -0.57735)
        
        double[] standardizedMean = new double[nDimensions];
        double[] standardizedStDev = new double[nDimensions];
        double[] data = Standardization.standardUnitNormalization(data0, 
            nDimensions, standardizedMean, standardizedStDev);
        
        GeometricMedianUnweightedFunction f0 
            = new GeometricMedianUnweightedFunction(data0, nDimensions);
        GeometricMedianUnweightedFunction f 
            = new GeometricMedianUnweightedFunction(data, nDimensions);
        GeometricMedian gm = new GeometricMedian();
        
        for (int ii = 0; ii < 3; ++ii) {
            switch(ii) {
                case 0:
                    init = new double[]{0, 0};            
                    break;
                case 1:
                    //init = new double[]{0, 12};
                    init = new double[]{0, 1.1547};
                    break;
                default:
                    init = f.calculateCentroid();;
                    break;
            }
            
            System.out.printf("\nbegin w/ init=(%.3f, %.3f); true=(%.3f, %.3f)\n",
                init[0], init[1], 
                standardizedExpected[0], standardizedExpected[1]);
                //expected[0], expected[1]);
            
            double min = gm.newtonsMethod2(f, init);

            System.out.println("in standardized units: min=" + min + " \n   coeffs=" +
                Arrays.toString(init));
            System.out.flush();
            
            // de-normalize the geometric-median "init" and recalculate min
            //    in the natural coordinates:
            init = Standardization.standardUnitDenormalization(init, nDimensions, 
                standardizedMean, standardizedStDev);
        
            min = f0.f(init);
            
            System.out.println("in data units: min=" + min + " \n   coeffs=" +
                Arrays.toString(init));
            System.out.flush();
            
            assertTrue(Math.abs(min - expectedDist) <= tol);

            for (int i = 0; i < init.length; ++i) {
                assertTrue(Math.abs(init[i] - expected[i]) < tol);
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
           
            double min = gm.newtonsMethod2(f, init);

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
