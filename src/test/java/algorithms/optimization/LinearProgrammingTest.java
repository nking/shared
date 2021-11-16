package algorithms.optimization;

import algorithms.matrix.MatrixUtil;
import algorithms.optimization.LinearProgramming.SlackForm;
import algorithms.optimization.LinearProgramming.StandardForm;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LinearProgrammingTest extends TestCase {
    
    public LinearProgrammingTest(String testName) {
        super(testName);
    }
    
    public void testPivot() {
        
        /*
        following Sect 29.3 of Cormen et al. "Introduction to Algorithms"
        
        except changing to "zero based" integer indexes on x, etc.
        
        Example of a linear program in Standard Form:
          maximize 3*x1 + x2 + 2*x3
          subject to:
                   x1 + x2 + 3*x3 .leq. 30
                   2*x1 + 2*x2 + 5*x3 .leq. 24
                   4*x1 + x2 + 2*x3 .leq. 36
                   x1, x2, x3 .geq. 0
        
        converted to Slack Form:
             maximize 3*x0 + x1 + 2*x2
             subject to:
                  x3 = 30 -   x0 -   x1 - 3*x2
                  x4 = 24 - 2*x0 - 2*x1 - 5*x2
                  x5 = 36 - 4*x0 -   x1 - 2*x2
                  x0, x1, x2, x3, x4, x5 .geq. 0
        */
        
        int[] bIndices = new int[]{3, 4, 5};
        int[] nIndices = new int[]{0, 1, 2};
        double[][] a = new double[3][];
        a[0] = new double[]{1, 1, 3};
        a[1] = new double[]{2, 2, 5};
        a[2] = new double[]{4, 1, 2};
        double[] b = new double[]{30,24,36};
        double[] c = new double[]{3, 1, 2};
        double v = 0;
        SlackForm slackForm = new SlackForm(nIndices, bIndices, a, b, c, v);
        
        double[] x = slackForm.computeBasicSolution();
        double[] expectedX = new double[]{0,0,0,30,24,36};
        assertEquals(expectedX.length, x.length);
        
        int i;
        double diff, tol=1e-7;
        for (i = 0; i < expectedX.length; ++i) {
            diff = Math.abs(expectedX[i] - x[i]);
            assertTrue(diff < tol);
        }
        
        // entering variable is the nonbasic.  use x1.  its index is nIndices is eIdx = 0
        // the leaving variable is the basic and that x6 (smallest soln when set x1 to 0).
        //     lIdx = 2 in bIndices for leaving variable x6.
        int eIdx = 0;
        int lIdx = 2;
        LinearProgramming lp = new LinearProgramming();
        SlackForm slackForm2 = lp.pivot(slackForm, lIdx, eIdx);
        
        x = slackForm2.computeBasicSolution();
        double z = slackForm2.evaluateObjective();
        
        //System.out.printf("pivoted eIdx=%d, lIdx=%d slackForm=\n%s\n", 
        //    eIdx, lIdx, slackForm2.toString());
        
        double expectedV = 27;
        double expectedZ = 27;
        double[] expectedB = new double[]{9.000, 21.000, 6.000};
        double[] expectedC = new double[]{0.250, 0.500, -0.750};
        int[] expectedBIndices = new int[]{0, 3, 4};
        int[] expectedNIndices = new int[]{1, 2, 5};
        double[][] expectedA = new double[3][];
        expectedA[0] = new double[]{0.250, 0.500, 0.250};
        expectedA[1] = new double[]{0.750, 2.500, -0.250};
        expectedA[2] = new double[]{1.500, 4.000, -0.500};
        expectedX = new double[]{9,0,0,21,6,0};
        
        assertExpected(tol, slackForm2, expectedV, expectedZ,
            expectedB, expectedC,
            expectedA, expectedBIndices, expectedNIndices, expectedX);
        
        //==========
        //result of entering var = x3 (eIdx=2) and leaving var=x5 (lIdx=2)
        slackForm = new SlackForm(expectedNIndices, expectedBIndices, expectedA, 
            expectedB, expectedC, expectedV);
        eIdx = 1;//x3
        lIdx = 2;//x5
        slackForm2 = lp.pivot(slackForm, lIdx, eIdx);
        x = slackForm2.computeBasicSolution();
        z = slackForm2.evaluateObjective();
        
        //System.out.printf("pivoted eIdx=%d, lIdx=%d obj=%.4f"
        //    + " slackForm=\n%s\n", 
        //    eIdx, lIdx, slackForm2.evaluateObjective(), slackForm2.toString());
        
        expectedV = 111/4.;//27.75
        expectedZ = 111/4.;
        expectedB = new double[]{8.25, 1.5, 17.25};
        expectedC = new double[]{1./16, -1./8, -11./16};
        expectedBIndices = new int[]{1-1, 3-1, 4-1};
        expectedNIndices = new int[]{2-1, 5-1, 6-1};
        expectedA = new double[3][];
        expectedA[0] = new double[]{1./16, -1./8, 5./16};
        expectedA[1] = new double[]{3./8, 1./4, -1./8};
        expectedA[2] = new double[]{-3./16, -5./8, 1./16};
        expectedX = new double[]{8.25, 0, 1.5, 17.25, 0, 0};
        
        assertExpected(tol, slackForm2, expectedV, expectedZ,
            expectedB, expectedC,
            expectedA, expectedBIndices, expectedNIndices, expectedX);
        
        //==========
        //result of entering var = x2 (eIdx=1) and leaving var=x3 (lIdx=1)
        slackForm = new SlackForm(expectedNIndices, expectedBIndices, expectedA, 
            expectedB, expectedC, expectedV);
        eIdx = 0;
        lIdx = 1;
        slackForm2 = lp.pivot(slackForm, lIdx, eIdx);
        x = slackForm2.computeBasicSolution();
        z = slackForm2.evaluateObjective();
        
        //System.out.printf("pivoted eIdx=%d, lIdx=%d slackForm=\n%s\n", 
        //    eIdx, lIdx, slackForm2.toString());
        
        expectedV = 28;
        expectedZ = 28;
        expectedB = new double[]{8, 4, 18};
        expectedC = new double[]{-1./6, -1./6, -2./3};
        expectedBIndices = new int[]{1-1, 2-1, 4-1};
        expectedNIndices = new int[]{3-1, 5-1, 6-1};
        expectedA = new double[3][];
        expectedA[0] = new double[]{-1./6, -1./6, 1./3};
        expectedA[1] = new double[]{8./3, 2./3, -1./3};
        expectedA[2] = new double[]{1./2, -1./2, 0};
        expectedX = new double[]{8, 4, 0, 18,0,0};
        
        assertExpected(tol, slackForm2, expectedV, expectedZ,
            expectedB, expectedC,
            expectedA, expectedBIndices, expectedNIndices, expectedX);
    }
    
    public void testSolveUsingSimplexMethod() {
        /*
        following Sect 29.3 of Cormen et al. "Introduction to Algorithms"
        
        except changing to "zero based" integer indexes on x, etc.
        
        Example of a linear program in Standard Form:
          maximize 3*x1 + x2 + 2*x3
          subject to:
                     x1 + x2 + 3*x3 .leq. 30
                   2*x1 + 2*x2 + 5*x3 .leq. 24
                   4*x1 + x2 + 2*x3 .leq. 36
                   x1, x2, x3 .geq. 0
        */
        double[] b = new double[]{30, 24, 36}; /*m*/
        double[] c = new double[]{3, 1, 2}; /*n*/
        double[][] a = new double[3][]; /*mxn*/
        a[0] = new double[]{1, 1, 3};
        a[1] = new double[]{2, 2, 5};
        a[2] = new double[]{4, 1, 2};
        
        double v = 0;
        
        StandardForm standForm = new LinearProgramming.StandardForm(
            a, b, c, v);
        
        long seed = System.nanoTime();
        LinearProgramming lp = new LinearProgramming(seed);       
        
        SlackForm soln = lp.solveUsingSimplexMethod(standForm);
        
        //System.out.printf("result=\n%s\n", soln.toString());
        
        
        double[] expectedX = new double[]{8, 4, 0, 18, 0, 0};
        double expectedZ = 28;
        
        assertNotNull(soln.x);
        assertEquals(expectedX.length, soln.x.length);
        
        int i;
        double diff, tol=1e-7;
        for (i = 0; i < expectedX.length; ++i) {
            diff = Math.abs(expectedX[i] - soln.x[i]);
            assertTrue(diff < tol);
        }
        assertTrue(Math.abs(expectedZ - soln.evaluateObjective()) < tol);
        
        /*
        xHat1,xHat2,xHat3,xHat4,xHat5,xHat6=(8, 4, 0, 18,0,0))

         Since we know iterations are finished:
          The original variables x1, x2, x3 are (8, 4, 0),
          so z = 3*8 + 1*4 + 2*0 = 28 is optimal.

          The slack variables tell us the amount of slack in each of the original inequalities
          in the Standard Form:
           constraints in Standard Form:
              x1 + x2 + 3*x3 .leq. 30  <== 8 + 4 * 0 .leq. 30  (slack var x4=18 which is 30-12)
              2*x1 + 2*x2 + 5*x3 .leq. 24
              4*x1 + x2 + 2*x3 .leq. 36
              x1, x2, x3 .geq. 0
        */
    }
    
    public void testAuxiliary() {
        /*
        Standard Form:
           maximize  2*x1 - 3*x2 + 3*x3
               subject to x1 +   x2 -   x3 .leq.  7
                          x1 - 2*x2 + 2*x3 .leq.  4
                          -x1 -  x2 +   x3 .leq. -7
                          x1, x2, and x3 .geq. 0
        */
        double[][] a = new double[3][];
        a[0] = new double[]{-1, -1, 1};
        a[1] = new double[]{-1, 2, -2};
        a[2] = new double[]{1, 1, -1};
        double[] b = new double[]{7, 4, -7};
        double[] c = new double[]{2, -3, 3};
        double v = 0;
        StandardForm standForm = new StandardForm(a, b, c, v);
        
        LinearProgramming lp = new LinearProgramming();
        SlackForm initSlackForm = lp.createAuxiliarySlackForm(standForm);
        double[] xBasicSoln = initSlackForm.computeBasicSolution();
        double eval = initSlackForm.evaluateObjective();
        double vHat = eval;
        //System.out.printf("createAuxiliarySlackForm=\n%s\n", initSlackForm.toString());
        //System.out.printf("eval=%.3f\n", eval);
        assertTrue(Math.abs(eval - 0) < 1e-11);
    }
    
    public void testAuxiliary2() {
        /*
        Standard Form:
           maximize  2*x1 - x2
               subject to 2*x1 -   x2  .leq.  2
                            x1 - 5*x2  .leq.  -4
                          x1, x2 .geq. 0
        */
        double[][] a = new double[2][];
        a[0] = new double[]{2, -1};
        a[1] = new double[]{1, -5};
        double[] b = new double[]{2, -4};
        double[] c = new double[]{2, -1};
        double v = 0;
        StandardForm standForm = new StandardForm(a, b, c, v);
        
        LinearProgramming lp = new LinearProgramming();
        SlackForm initSlackForm = lp.createAuxiliarySlackForm(standForm);
        double[] xBasicSoln = initSlackForm.computeBasicSolution();
        double eval = initSlackForm.evaluateObjective();
        boolean isFeasible = initSlackForm.isFeasible();
        System.out.printf("createAuxiliarySlackForm=\n%s\n", initSlackForm.toString());
        System.out.printf("eval=%.3f isFeasible=%b\n", eval, isFeasible);
        assertTrue(Math.abs(eval - 0) < 1e-11);
        
        double[][] expectedA = new double[2][];
        expectedA[0] = new double[]{-1, 2, -1};
        expectedA[1] = new double[]{-1, 1, -5};
        double[] expectedB = new double[]{2, -4};
        double[] expectedC = new double[]{-1, 0, 0};
        double expectedV = 0;
        double[] expectedX = new double[]{0.000, 0.000, 0.000, 2.000, -4.000};
        int[] expectedNIndices = new int[]{0, 1, 2};
        int[] expectedBIndices = new int[]{3, 4};
        boolean expectedIsFeasible = false;
        assertEquals(expectedA.length, initSlackForm.a.length);
        assertEquals(expectedA[0].length, initSlackForm.a[0].length);
        assertEquals(expectedB.length, initSlackForm.b.length);
        assertEquals(expectedC.length, initSlackForm.c.length);
        assertEquals(expectedNIndices.length, initSlackForm.nIndices.length);
        assertEquals(expectedBIndices.length, initSlackForm.bIndices.length);
        assertTrue(Arrays.equals(expectedBIndices, initSlackForm.bIndices));
        assertTrue(Arrays.equals(expectedNIndices, initSlackForm.nIndices));
        double diff, tol = 1e-7;
        int i, j;
        for (i = 0; i < expectedA.length; ++i) {
            for (j = 0; j < expectedA[i].length; ++j) {
                diff = Math.abs(expectedA[i][j] - initSlackForm.a[i][j]);
                assertTrue(diff < tol);
            }
        }
        for (i = 0; i < expectedC.length; ++i) {
            diff = Math.abs(expectedC[i] - initSlackForm.c[i]);
            assertTrue(diff < tol);
        }
        for (i = 0; i < expectedB.length; ++i) {
            diff = Math.abs(expectedB[i] - initSlackForm.b[i]);
            assertTrue(diff < tol);
        }
        assertEquals(expectedIsFeasible, isFeasible);
        diff = Math.abs(expectedV - initSlackForm.v);
        assertTrue(diff < tol);
        for (i = 0; i < expectedX.length; ++i) {
            diff = Math.abs(expectedX[i] - initSlackForm.x[i]);
            assertTrue(diff < tol);
        }
        
        //=======
        SlackForm slackForm = lp.initializeSimplex(standForm);
    }
        
        /*
        Example that is not yet in Standard Form:
           minimize -2*x1 + 3*x2
           subject to x1 + x2 = 7
                      x1 - 2*x2 .leq. 4
                      x1 .geq. 0
        Converted to Standard Form:
           maximize  2*x1 - 3*x2 + 3*x3
               subject to x1 + x2 - x3 .leq. 7
                          x1 - 2*x2 + 2*x3 .leq. 4
                          -x1 - x2 + x3 .leq. -7
                          x1, x2, and x3 .geq. 0
        */
        
        double[][] a = new double[2][];
        a[0] = new double[]{1, 1};
        a[1] = new double[]{1, -2};
        double[] b = new double[]{7, 4};
        double[] c = new double[]{-2, 3};
        int[] constraintComparisons = new int[]{0, -1};
        boolean isMaximization = false;
        boolean[] nonnegativityConstraints = new boolean[]{true, false};
    
        StandardForm standForm = LinearProgramming.convertLinearProgramToStandardForm(
            isMaximization, a, b, c, 
            constraintComparisons, nonnegativityConstraints);
        
        //System.out.printf("standForm=\n%s\n", standForm.toString());
        /*
        maximize  2*x1 - 3*x2 + 3*x3
        subject to 
            x1 +   x2 -   x3 .leq. 7
            x1 - 2*x2 + 2*x3 .leq. 4
           -x1 -   x2 +   x3 .leq. -7
            x1, x2, and x3 .geq. 0
        */
        double[][] expectedA = new double[3][];
        expectedA[0] = new double[]{1, 1, -1};
        expectedA[1] = new double[]{1, -2, 2};
        expectedA[2] = new double[]{-1, -1, 1};
        double[] expectedB = new double[]{7, 4, -7};
        double[] expectedC = new double[]{2,-3, 3};
        double expectedV = 0;
        
        assertEquals(expectedA.length, standForm.a.length);
        assertEquals(expectedA[0].length, standForm.a[0].length);
        assertEquals(expectedB.length, standForm.b.length);
        assertEquals(expectedC.length, standForm.c.length);
        
        int i, j;
        double diff, tol=1e-7;
        for (i = 0; i < expectedB.length; ++i) {
            diff = Math.abs(expectedB[i] - standForm.b[i]);
            assertTrue(diff < tol);
        }
        for (i = 0; i < expectedC.length; ++i) {
            diff = Math.abs(expectedC[i] - standForm.c[i]);
            assertTrue(diff < tol);
        }
        diff = Math.abs(expectedV - standForm.v);
        assertTrue(diff < tol);
        for (i = 0; i < expectedA.length; ++i) {
            for (j = 0; j < expectedA[i].length; ++j) {
                diff = Math.abs(expectedA[i][j] - standForm.a[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        SlackForm slackForm = LinearProgramming.convertConstraints(standForm);
        //System.out.printf("standForm=\n%s\n", standForm.toString());
        //System.out.printf("slackForm=\n%s\n", slackForm.toString());
        
        // slackform b should be same as original linear problem
        assertEquals(b.length, slackForm.b.length);
        assertEquals(a.length, slackForm.a.length);
        for (i = 0; i < b.length; ++i) {
            diff = Math.abs(b[i] - slackForm.b[i]);
            assertTrue(diff < tol);
        }
        for (i = 0; i < a.length; ++i) {
            for (j = 0; j < a[i].length; ++j) {
                diff = Math.abs(a[i][j] - slackForm.a[i][j]);
                assertTrue(diff < tol);
            }
        }
    }
    
    public void testConvertConstraints() {
        /*
        sample from https://walkccc.me/CLRS/Chap29/29.3/
                Linear Program in Standard Form:
                maximize: 18x1 + 12.5x2
                subject to:
                            x1 + x2 .leq. 20
                            x1      .leq. 12
                                 x2 .leq. 16
                            x1,x2 .geq. 0
        
                Converted to Linear Program in Slack Form:
                maximize: 18x1 + 12.5x2
                subject to:
                            x3 = 20 - x1 - x2
                            x4 = 12 - x1
                            x5 = 16 - x2
                            x1,x2,x3, x4, x5 .geq. 0
        */
        int n = 2;
        int m = 3;
        double[] c = new double[]{18, 12.5};
        double[] b = new double[]{20, 12, 16};
        int[] nIndices = new int[]{0, 1};
        double v = 0;
        double[][] a = new double[m][];
        a[0] = new double[]{1, 1};
        a[1] = new double[]{1, 0};
        a[2] = new double[]{0, 1};
        
        double[] cExpected = Arrays.copyOf(c, c.length);
        double[] bExpected = Arrays.copyOf(b, b.length);
        int[] nIndicesExpected = Arrays.copyOf(nIndices, nIndices.length);
        int[] bIndicesExpected = new int[]{2, 3, 4};
        double vExpected = 0;
        double[][] aExpected = MatrixUtil.copy(a);
        
        StandardForm standForm = new StandardForm(a, b, c, v);
        SlackForm slackForm = LinearProgramming.convertConstraints(standForm);
        
        double diff, tol=1e-7;
        int i;
        assertEquals(bIndicesExpected.length,  slackForm.bIndices.length);
        assertTrue(Arrays.equals(bIndicesExpected, slackForm.bIndices));
        assertEquals(nIndicesExpected.length,  slackForm.nIndices.length);
        assertTrue(Arrays.equals(nIndicesExpected, slackForm.nIndices));
        
        assertEquals(cExpected.length, slackForm.c.length);
        assertEquals(bExpected.length, slackForm.b.length);
        assertEquals(aExpected.length, slackForm.a.length);
        assertEquals(aExpected[0].length, slackForm.a[0].length);
        
        for (i = 0; i < cExpected.length; ++i) {
            diff = Math.abs(cExpected[i] - slackForm.c[i]);
            assertTrue(diff < tol);
        }
        for (i = 0; i < bExpected.length; ++i) {
            diff = Math.abs(bExpected[i] - slackForm.b[i]);
            assertTrue(diff < tol);
        }
        int j;
        for (i = 0; i < aExpected.length; ++i) {
            for (j = 0; j < aExpected[i].length; ++j) {
                diff = Math.abs(aExpected[i][j] - slackForm.a[i][j]);
                assertTrue(diff < tol);
            }
        }
        diff = Math.abs(vExpected - slackForm.v);
        assertTrue(diff < tol);
    }
    
    public void testSolveUsingSimplexMethod2() {
        System.out.println("testSolveUsingSimplexMethod2");
        /*
        example from "Operations Research. Linear Programming",
           pg 79, sect 2.9.  OpenCourseWare, UPV/EHU
        
        Linear model
           max z = 6x1 +4x2 +5x3 +5x4 +0x5 +0x6 +0x7
           subject to
              x1 + x2 +  x3 + x4 + x5           = 3
             2x1 + x2 + 4x3 + x4      + x6      = 4
             x1 + 2x2 − 2x3 + 3x4          + x7 = 10
             x1, x2, x3, x4, x5, x6, x7, >= 0
        
        x = [1, 0, 0, 2, 0, 0, 3, 16]
        z = 16
        */
        
        double[] b = new double[]{3, 4, 10};
        double[] c = new double[]{6, 4, 5, 5, 0, 0, 0};
        double[][] a = new double[3][];
        a[0] = new double[]{1, 1, 1, 1, 1, 0, 0};
        a[1] = new double[]{2, 1, 4, 1, 0, 1, 0};
        a[2] = new double[]{1, 2, -2, 3, 0, 0, 1};
        double v = 0;
        int[] constraintComparisons = new int[]{0, 0, 0};
        boolean[] nonnegativityConstraints = new boolean[]{true, true, true,
            true, true, true, true};
        boolean isMaximization = true;
        
        LinearProgramming lp = new LinearProgramming();
        StandardForm standForm = LinearProgramming.convertLinearProgramToStandardForm(
            isMaximization, a, b, c, constraintComparisons, nonnegativityConstraints);
        
        //System.out.printf("standForm=\n%s\n", standForm.toString());
        
        SlackForm soln = lp.solveUsingSimplexMethod(standForm);
        
        //System.out.printf("soln=\n%s\n", soln.toString());
        //System.out.printf("z=%.3f\n", soln.evaluateObjective());
        
        /*
        v=16.000
        c=-1.000, -3.000, -4.000, -1.000, 0.000, -4.000, -1.000 
        b=1.000, 2.000, 3.000 
        a=
        0.000, 3.000, -1.000, 1.000, 0.000, -1.000, 1.000 
        1.000, -2.000, 2.000, -1.000, 0.000, 2.000, -1.000 
        -1.000, 1.000, -5.000, 2.000, 1.000, -5.000, 2.000 
        x=1.000, 0.000, 0.000, 2.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.000 
        nIndices=[1, 2, 4, 5, 6, 7, 8]
        bIndices=[0, 3, 9]
        STATE=OPTIMAL
        0, 1, 2, 3, --, 5, 6
        */
        
        // degenerate, ~3 eqns, 7 unknowns
        double tol = 1e-7, diff;
        double expectedZ = 16;
        double[] expectedX = new double[]{1, 0, 0, 2, 0, 0, 3, 16};
                                        //1, 0, 0, 2, 0, 0, 0, 0, 0, 3 
        int i;
        /*for (i = 0; i < expectedX.length; ++i) {
            diff = Math.abs(expectedX[i] - soln.x[i]);
            assertTrue(diff < tol);
        }*/
        diff = Math.abs(expectedZ - soln.evaluateObjective());
        assertTrue(diff < tol);
       
    }
    
    /*
    unique soln as have 3 eqns and 2 unknowns, though could need best fit
       maximize x1 + x2
       subject to 
            x1 +  2*x2  <= 4
           4*x1 +  2*x1  <=  12
            -x1 +  x2  <= 1
         x1, x2 >= 0
          
    x=[0.67, 2.67]
    z = 3.33
    */
    
    public void testSolveUsingSimplexMethod3() {
        System.out.println("testSolveUsingSimplexMethod3");
        /*
        example from "Operations Research. Linear Programming",
           pg 63, sect 2.7.  OpenCourseWare, UPV/EHU
        
        unique optimal solution exists.
        
        Linear model
           max z = 6*x1 + 4*x2 + 5*x3 + 5*x4 
           subject to
              x1 +  x2 +  x3 +  x4 ≤  3 
             2x1 +  x2 + 4x3 +  x4 ≤  4 
              x1 + 2x2 − 2x3 + 3x4 ≤ 10 
              x1,x2,x3,x4 ≥ 0
        
        Standard Form:
            max z = 6x1 + 4x2 + 5x3 + 5x4 + 0x5 + 0x6 + 0x7 
            subject to
              x1 +  x2 +  x3 +  x4 +  x5               =   3
             2x1 +  x2 + 4x3 +  x4        + x6         =   4
              x1 + 2x2 − 2x3 + 3x4              +   x7 =  10
              x1, x2, x3, x4, x5, x6, x7 ≥0
        x = [1, 0, 0, 2, 0, 0, 3, 16]
        z = 16
        */
        
        double[] b = new double[]{3, 4, 10};
        double[] c = new double[]{6, 4, 5, 5};
        double[][] a = new double[3][];
        a[0] = new double[]{1, 1, 1, 1};
        a[1] = new double[]{2, 1, 4, 1};
        a[2] = new double[]{1, 2, -2, 3};
        double v = 0;
        int[] constraintComparisons = new int[]{-1, -1, -1};
        boolean[] nonnegativityConstraints = new boolean[]{true, true, true,
            true};
        boolean isMaximization = true;
        
        LinearProgramming lp = new LinearProgramming();
        StandardForm standForm = LinearProgramming.convertLinearProgramToStandardForm(
            isMaximization, a, b, c, constraintComparisons, nonnegativityConstraints);
        
        //System.out.printf("standForm=\n%s\n", standForm.toString());
        
        SlackForm soln = lp.solveUsingSimplexMethod(standForm);
        
        //System.out.printf("soln=\n%s\n", soln.toString());
        //System.out.printf("z=%.3f\n", soln.evaluateObjective());
        
        /*
        [junit] soln=
        [junit] v=16.000
        [junit] c=-1.000, -3.000, -4.000, -1.000 
        [junit] b=1.000, 2.000, 3.000 
        [junit] a=
        [junit] 0.000, 3.000, -1.000, 1.000 
        [junit] 1.000, -2.000, 2.000, -1.000 
        [junit] -1.000, 1.000, -5.000, 2.000 
        [junit] x=1.000, 0.000, 0.000, 2.000, 0.000, 0.000, 3.000 
        [junit] nIndices=[1, 2, 4, 5]
        [junit] bIndices=[0, 3, 6]
        [junit] STATE=OPTIMAL
        [junit] 
        [junit] z=16.000
        */
        
        // degenerate, ~3 eqns, 7 unknowns
        double tol = 1e-7, diff;
        double expectedZ = 16;
        double[] expectedX = new double[]{1, 0, 0, 2, 0, 0, 3, 16};
                                        //1, 0, 0, 2, 0, 0, 0, 0, 0, 3 
        int i;
        /*for (i = 0; i < expectedX.length; ++i) {
            diff = Math.abs(expectedX[i] - soln.x[i]);
            assertTrue(diff < tol);
        }*/
        diff = Math.abs(expectedZ - soln.evaluateObjective());
        assertTrue(diff < tol);
       
    }

    private void assertExpected(double tol,
        SlackForm slackForm2, double expectedV, 
        double expectedZ,
        double[] expectedB, double[] expectedC, double[][] expectedA, 
        int[] expectedBIndices, int[] expectedNIndices, double[] expectedX) {
        
        int i;
        double diff;
        assertTrue(Math.abs(expectedV - slackForm2.v) < tol);
        assertEquals(expectedB.length, slackForm2.b.length);
        assertEquals(expectedC.length, slackForm2.c.length);
        assertEquals(expectedBIndices.length, slackForm2.bIndices.length);
        assertEquals(expectedNIndices.length, slackForm2.nIndices.length);
        assertEquals(expectedA.length, slackForm2.a.length);
        assertEquals(expectedA[0].length, slackForm2.a[0].length);
        assertEquals(expectedX.length, slackForm2.x.length);
        
        for (i = 0; i < expectedB.length; ++i) {
            diff = Math.abs(expectedB[i] - slackForm2.b[i]);
            assertTrue(diff < tol);
        }
        for (i = 0; i < expectedC.length; ++i) {
            diff = Math.abs(expectedC[i] - slackForm2.c[i]);
            assertTrue(diff < tol);
        }
        for (i = 0; i < expectedBIndices.length; ++i) {
            assertEquals(expectedBIndices[i], slackForm2.bIndices[i]);
        }
        for (i = 0; i < expectedNIndices.length; ++i) {
            assertEquals(expectedNIndices[i], slackForm2.nIndices[i]);
        }
        int j;
        for (i = 0; i < slackForm2.a.length; ++i) {
            for (j = 0; j < slackForm2.a[i].length; ++j) {
                diff = Math.abs(expectedA[i][j] - slackForm2.a[i][j]);
                assertTrue(diff < tol);
            }
        }
        for (i = 0; i < expectedX.length; ++i) {
            diff = Math.abs(expectedX[i] - slackForm2.x[i]);
            assertTrue(diff < tol);
        }
        diff = Math.abs(expectedZ - slackForm2.evaluateObjective());
        assertTrue(diff < tol);
    }
}
