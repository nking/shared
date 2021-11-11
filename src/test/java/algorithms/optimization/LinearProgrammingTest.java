package algorithms.optimization;

import algorithms.optimization.LinearProgramming.SlackForm;
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
        
        Slack Form example:
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
        int lIdx = 2;
        int eIdx = 0;
        LinearProgramming lp = new LinearProgramming();
        SlackForm slackForm2 = lp.pivot(slackForm, lIdx, eIdx);
        
        x = slackForm2.computeBasicSolution();
        double z = slackForm2.evaluateObjective();
        
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
        
        assertTrue(Math.abs(expectedV - slackForm2.v) < tol);
        assertEquals(expectedB.length, slackForm2.b.length);
        assertEquals(expectedC.length, slackForm2.c.length);
        assertEquals(expectedBIndices.length, slackForm2.bIndices.length);
        assertEquals(expectedNIndices.length, slackForm2.nIndices.length);
        assertEquals(expectedA.length, slackForm2.a.length);
        assertEquals(expectedA[0].length, slackForm2.a[0].length);
        assertEquals(expectedX.length, x.length);
        
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
        for (i = 0; i < a.length; ++i) {
            for (j = 0; j < a[i].length; ++j) {
                diff = Math.abs(expectedA[i][j] - slackForm2.a[i][j]);
                assertTrue(diff < tol);
            }
        }
        for (i = 0; i < expectedC.length; ++i) {
            diff = Math.abs(expectedX[i] - x[i]);
            assertTrue(diff < tol);
        }
        assertTrue(Math.abs(expectedZ - z) < tol);
    }
}
