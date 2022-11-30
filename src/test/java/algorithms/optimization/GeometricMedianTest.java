package algorithms.optimization;

import algorithms.statistics.Standardization;
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
        double expectedDistStandardized = 2.;
        
        double tol = 0.01;
                
        double[] init;
        // the geometric-median in original data reference frames:
        double[] expected = new double[]{2, 2};
        // the geometric-median in standardized coordinate reference frames:
        double[] standardizedExpected = new double[]{0., 0};
        
        double[] standardizedMean = new double[nDimensions];
        double[] standardizedStDev = new double[nDimensions];
        double[] data = Standardization.standardUnitNormalization(data0, 
            nDimensions, standardizedMean, standardizedStDev);
        
        AbstractGeometricMedianFunction f0, f;
        
        for (int type = 0; type < 2; ++type) {
        
            if (type == 0) {
                f0 = new GeometricMedianUnweightedFunction(data0, nDimensions);
                f = new GeometricMedianUnweightedFunction(data, nDimensions);
            } else {
                // consider using 1/n.
                double[] eta = new double[data0.length/nDimensions];
                Arrays.fill(eta, 1.);
                f0 = new GeometricMedianWeightedFunction(data0, nDimensions, eta);
                f = new GeometricMedianWeightedFunction(data, nDimensions, eta);
            }
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

                System.out.printf("\ntype=%d, ii=%d) begin w/ init=(%.3f, %.3f); standardizedExpected=(%.3f, %.3f)\n",
                    type, ii, init[0], init[1], 
                    standardizedExpected[0], standardizedExpected[1]);

                double min;
                if (type == 0) {
                    min = gm.newtonsMethod2(f, init);
                } else {
                    min = gm.newtonsThenVardiZhang(
                        (GeometricMedianWeightedFunction)f, init);
                }

                System.out.printf("in standardized units: min=%.4e \n   coeffs=%s\n   expectedMin==%.4e\n",
                    min, Arrays.toString(init), expectedDistStandardized);
                System.out.flush();

                // de-normalize the geometric-median "init" and recalculate min
                //    in the natural coordinates:
                init = Standardization.standardUnitDenormalization(init, nDimensions, 
                    standardizedMean, standardizedStDev);

                min = f0.f(init);

                System.out.printf("in data units: min=%.4e   coeffs=%s\n   expected=%.4e,   %s\n", 
                    min, Arrays.toString(init), expectedDist,
                    AbstractGeometricMedianFunction.toString(expected));
                System.out.flush();

                assertTrue(Math.abs(min - expectedDist) <= tol);

                for (int i = 0; i < init.length; ++i) {
                    double diff = Math.abs(init[i] - expected[i]);
                    assertTrue(diff < tol);
                }
            }
        }
        
    }
    
    public void test1() {
        
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
        double expectedDistStandardized = 1.7321;
        
        double tol = 0.01;
                
        double[] init;
        double[] expected = new double[]{0, 0};
        double[] standardizedExpected = new double[]{0., -0.5773};
        // expected in standardized coords = (0-0)/0  , (0-4)/6.9282 = (0, -0.57735)
        
        double[] standardizedMean = new double[nDimensions];
        double[] standardizedStDev = new double[nDimensions];
        double[] data = Standardization.standardUnitNormalization(data0, 
            nDimensions, standardizedMean, standardizedStDev);
        
        AbstractGeometricMedianFunction f0, f;
        
        for (int type = 0; type < 2; ++type) {
        
            if (type == 0) {
                f0 = new GeometricMedianUnweightedFunction(data0, nDimensions);
                f = new GeometricMedianUnweightedFunction(data, nDimensions);
            } else {
                // consider using 1/n.
                double[] eta = new double[data0.length/nDimensions];
                Arrays.fill(eta, 1.0);
                f0 = new GeometricMedianWeightedFunction(data0, nDimensions, eta);
                f = new GeometricMedianWeightedFunction(data, nDimensions, eta);
            }
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

                System.out.printf("\ntype=%d, ii=%d) begin w/ init=(%.3f, %.3f); standardizedExpected=(%.3f, %.3f)\n",
                    type, ii, init[0], init[1], 
                    standardizedExpected[0], standardizedExpected[1]);

                double min;
                if (type == 0) {
                    min = gm.newtonsMethod2(f, init);
                } else {
                    min = gm.newtonsThenVardiZhang(
                        (GeometricMedianWeightedFunction)f, init);
                }

                System.out.printf("in standardized units: min=%.4e \n   coeffs=%s\n   expectedMin==%.4e\n",
                    min, Arrays.toString(init), expectedDistStandardized);
                System.out.flush();

                // de-normalize the geometric-median "init" and recalculate min
                //    in the natural coordinates:
                init = Standardization.standardUnitDenormalization(init, nDimensions, 
                    standardizedMean, standardizedStDev);

                min = f0.f(init);

                System.out.printf("in data units: min=%.4e   coeffs=%s\n  expected=%.4e,   %s\n", 
                    min, Arrays.toString(init), expectedDist,
                    AbstractGeometricMedianFunction.toString(expected));
                System.out.flush();

                assertTrue(Math.abs(min - expectedDist) <= tol);

                for (int i = 0; i < init.length; ++i) {
                    double diff = Math.abs(init[i] - expected[i]);
                    assertTrue(diff < tol);
                }
            }
        }
    }
    
    public void test2() {
        
        System.out.println("test2");
        
        /*
        from "Noniterative Solution of Some Fermat-Weber Location Problems"
        https://www.researchgate.net/publication/220418355_Noniterative_Solution_of_Some_Fermat-Weber_Location_Problems
        -- Input: (-20, 48), (-20, -48), (20, 0), (59, 0)
           Output: Geometric Median = (20,0) with minimum distance = 163.964
           NOTE: start with (44, 0) to test for a Weiszfeld problem w/ demand points
        */
        /*
        standardized: Mean = new double[]{9.75, 0}; 
        standardized: StDev = new double[]{37.863, 39.192};
        standardized: new double[]{-0.7857, 1.225, -0.786, -1.225, 0.271, 0, 1.301, 0};
        */
        double[] data0 = new double[]{-20, 48, -20, -48, 20, 0, 59, 0}; //
        int nDimensions = 2;
                
        double[] init;
        // the geometric-median in original data frames:
        double[] expected = new double[]{20., 0};
        // the geometric-median in standardized coordinate frames:
        double[] standardizedExpected = new double[]{0.2707, 0};
        
        double expectedDist = 163.964; // in standardized ~ 4.265
        
        double tol = 0.01;
        
        double[] standardizedMean = new double[nDimensions];
        double[] standardizedStDev = new double[nDimensions];
        //-0.7857, 1.225, -0.786, -1.225, 0.271, 0, 1.301, 
        double[] data = Standardization.standardUnitNormalization(data0, 
            nDimensions, standardizedMean, standardizedStDev);
        
        AbstractGeometricMedianFunction f0, f;
        
        for (int type = 1; type < 2; ++type) {
        
            if (type == 0) {
                // fails for unweighted when derivative becomes 0 and point is
                //   not the minimum, but is a "demand point"
                f0 = new GeometricMedianUnweightedFunction(data0, nDimensions);
                f = new GeometricMedianUnweightedFunction(data, nDimensions);
            } else {
                // consider using 1/n instead of 1.
                double[] eta = new double[data0.length/nDimensions];
                Arrays.fill(eta, 1.0);
                //double[] eta = new double[]{13., 13., 5., 5.};
                f0 = new GeometricMedianWeightedFunction(data0, nDimensions, eta);
                f = new GeometricMedianWeightedFunction(data, nDimensions, eta);
            }
            
            GeometricMedian gm = new GeometricMedian();

            for (int ii = 0; ii < 6; ++ii) {
                switch(ii) {
                    case 0:
                        //init = new double[]{0, 0};
                        init = new double[]{-0.2575, 0};
                        break;
                    case 1:
                        //init = new double[]{44, 0};
                        init = new double[]{0.905, 0};
                        break;
                    case 2:
                        //init = new double[]{3, 45};
                        init = new double[]{-0.178, 1.148};
                        break;
                    case 3:
                        //init = new double[]{0, 0}; 
                        init = new double[]{-0.258, 0.};
                        break;
                    case 4:
                        //init = new double[]{20, 0};
                        init = new double[]{0.2707, 0.};
                        break;
                    default:
                        // (9.75, 0)
                        init = f.calculateCentroid();;
                        break;
                }

                System.out.printf("\n%d) begin w/ init=(%.3f, %.3f); true=(%.3f, %.3f)\n",
                    type, init[0], init[1], 
                    standardizedExpected[0], standardizedExpected[1]);
                    //expected[0], expected[1]);

                double min;
                if (type == 0) {
                    min = gm.newtonsMethod2(f, init);
                } else {
                    min = gm.newtonsThenVardiZhang(
                        (GeometricMedianWeightedFunction)f, init);
                }

                System.out.println("in standardized units: min=" + min + " \n   coeffs=" +
                    Arrays.toString(init));
                System.out.flush();

                boolean b = gm.verify(f, init, tol);
                assertTrue(b);

                // de-normalize the geometric-median "init" and recalculate min
                //    in the natural coordinates:
                init = Standardization.standardUnitDenormalization(init, nDimensions, 
                    standardizedMean, standardizedStDev);

                min = f0.f(init);

                System.out.printf("in data units: min=%.4e   coeffs=%s\n  expected=%.4e,   %s", 
                    min, Arrays.toString(init), expectedDist,
                    AbstractGeometricMedianFunction.toString(expected));
                System.out.flush();

                assertTrue(Math.abs(min - expectedDist) <= tol * expectedDist);

                for (int i = 0; i < init.length; ++i) {
                    double diff = Math.abs(init[i] - expected[i]);
                    double tol2 = 0.03 * Math.abs(standardizedStDev[i]);
                    assertTrue(diff <= tol2);
                }
            }
        }
    }
    
    public void test3() {
        
        System.out.println("test3");
        
        // testing the weighted with multiplicities that are > 1
        
        /*
        from "Noniterative Solution of Some Fermat-Weber Location Problems"
        https://www.researchgate.net/publication/220418355_Noniterative_Solution_of_Some_Fermat-Weber_Location_Problems
        -- Input: (-20, 48), (-20, -48), (20, 0), (59, 0)
           Output: Geometric Median = (20,0) with minimum distance = 163.964
           NOTE: start with (44, 0) to test for a Weiszfeld problem w/ demand points
        */
        /*
        standardized: Mean = new double[]{9.75, 0}; 
        standardized: StDev = new double[]{37.8626905189, 39.192};
        standardized: new double[]{
        //   -0.7857339,1.2247397428, -0.7857339,-1.2247397428, 
        //       0.2707150458545328,0, 1.30075278,0};
        */
        double[] data0 = new double[]{-20, 48, -20, -48, 20, 0, 59, 0}; //
        int nDimensions = 2;
                
        double[] init;
        // the geometric-median in original data frames:
        double[] expected = new double[]{0., 0};
        // the geometric-median in standardized coordinate frames:
        double[] standardizedExpected = new double[]{-0.25750943, 0}; // weights not applied to standardization
        
        double expectedDist = 1747.014152324044;
        double expectedDistStandardized = 45.11109783739542; // weighted
        
        /* ==== found this soln which is close, but not better ====
        expectedDistStandardized = 45.10843371253345;
        expectedDist = 1747.0935;
        expected = new double[]{-0.6609909, 0};
        standardizedExpected = new double[]{-0.2749670125987123, 0};*/
        
        double tol = 0.01;
        
        double[] standardizedMean = new double[nDimensions];
        double[] standardizedStDev = new double[nDimensions];
      
        double[] data = Standardization.standardUnitNormalization(data0, 
            nDimensions, standardizedMean, standardizedStDev);
        
        GeometricMedianWeightedFunction f0, f;

        double[] eta = new double[]{13., 13., 5., 5.};
        f0 = new GeometricMedianWeightedFunction(data0, nDimensions, eta);
        f = new GeometricMedianWeightedFunction(data, nDimensions, eta);

        GeometricMedian gm = new GeometricMedian();

        for (int ii = 0; ii < 6; ++ii) {
            switch (ii) {
                case 0:
                    //init = new double[]{0, 0};
                    init = new double[]{-0.2575, 0};
                    break;
                case 1:
                    //init = new double[]{44, 0};
                    init = new double[]{0.905, 0};
                    break;
                case 2:
                    //init = new double[]{3, 45};
                    init = new double[]{-0.178, 1.148};
                    break;
                case 3:
                    //init = new double[]{0, 0}; 
                    init = new double[]{-0.258, 0.};
                    break;
                case 4:
                    //init = new double[]{20, 0};
                    init = new double[]{0.2707, 0.};
                    break;
                default:
                    // (9.75, 0)
                    init = f.calculateCentroid();
                    ;
                    break;
            }

            System.out.printf("\nii=%d) begin w/ init=(%.3f, %.3f); standardizedExpected=(%.3f, %.3f)\n",
                    ii, init[0], init[1],
                    standardizedExpected[0], standardizedExpected[1]);
            //expected[0], expected[1]);

            double min = gm.newtonsThenVardiZhang(f, init);

            System.out.printf("in standardized units: min=%.4e \n   coeffs=%s\n   expectedMin==%.4e\n",
                min, Arrays.toString(init), expectedDistStandardized);
            System.out.flush();
                
            boolean b = gm.verify(f, init, tol);
            assertTrue(b);

            // de-normalize the geometric-median "init" and recalculate min
            //    in the natural coordinates:
            init = Standardization.standardUnitDenormalization(init, nDimensions,
                    standardizedMean, standardizedStDev);
            
            min = f0.f(init);

            System.out.printf("in data units: min=%.5e   coeffs=%s  expected=%.5e,   %s\n",
                    min, Arrays.toString(init), expectedDist,
                    AbstractGeometricMedianFunction.toString(expected));
            System.out.flush();

            assertTrue(Math.abs(min - expectedDist) <= tol * expectedDist);

            /*for (int i = 0; i < init.length; ++i) {
                double diff = Math.abs(init[i] - expected[i]);
                System.out.printf("diff=%.4e\n", diff);
                assertTrue(diff <= tol);
            }*/
        }

    }
}
