package algorithms.misc;

import algorithms.util.FormatArray;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * @author nichole
 */
public class PolynomialRootSolverTest extends TestCase {
    
    public PolynomialRootSolverTest(String testName) {
        super(testName);
    }

    public void test00() throws IOException {
        //[4, 3, 2, 1] for 4*x^3 + 3*x^2 + 2*x + 1 = 0.
        //-3(p^2) + 4p - 1 => [-3, 4, -1]
        Complex[] a = PolynomialRootSolver.solveUsingMPSolve(new double[]{-3, 4, -1});

        //$6 * (p^4) - 14 * (p^3) + 10*p^2 - 2p = 0$
        Complex[] b = PolynomialRootSolver.solveUsingMPSolve(new double[]{6, -4, +10, -2});

        int w = 1;
    }
    public void test0() throws NotConvergedException, Exception {
        
        double tol = 1e-5;
        
        double[] coeffs = new double[]{1, 1, 0, -2};
        double[] expected = new double[]{1};
        
        Set<Complex> expectedC = new HashSet<Complex>();
        expectedC.add(new Complex(-1, 1));
        expectedC.add(new Complex(-1, -1));
        expectedC.add(new Complex(1, 0));
        
        /*
        >>> numpy.roots([1, 1, 0, -2])
            array([-1.+1.j, -1.-1.j,  1.+0.j])
        */
     
        assertTrue(Arrays.equals(new int[]{0, 1, 3}, PolynomialRootSolver.nonzero(coeffs, tol)));
        
        Complex[] roots = PolynomialRootSolver.solveUsingCompanionMatrix(coeffs);
        assertEquals(expectedC.size(), roots.length);
        
        double diff0, diff1;
        for (Complex rE : expectedC) {
            for (int i = 0; i < roots.length; ++i) {
                Complex r = roots[i];
                if (r == null) continue;
                diff0 = Math.abs(r.re() - rE.re());
                diff1 = Math.abs(r.im() - rE.im());
                if (diff0 < tol && diff1 < tol) {
                    roots[i] = null;
                }
            }
        }
        int u = 0;
        for (int i = 0; i < roots.length; ++i) {
            Complex r = roots[i];
            if (r != null) u++;
        }
        assertEquals(0, u);
        
        double[] realRoots = PolynomialRootSolver.solveForRealUsingCompanionMatrix(coeffs, tol);
        assertEquals(expected.length, realRoots.length);
        
        for (int i = 0; i < expected.length; ++i) {
            diff0 = Math.abs(expected[i] - realRoots[i]);
            assertTrue(diff0 < tol);
        }
       
        Complex[] roots2 = PolynomialRootSolver.solveUsingMPSolve(coeffs);
        //System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        for (Complex rE : expectedC) {
            for (int i = 0; i < roots2.length; ++i) {
                Complex r = roots2[i];
                if (r == null) continue;
                diff0 = Math.abs(r.re() - rE.re());
                diff1 = Math.abs(r.im() - rE.im());
                if (diff0 < tol && diff1 < tol) {
                    roots2[i] = null;
                }
            }
        }
        u = 0;
        for (int i = 0; i < roots2.length; ++i) {
            Complex r = roots2[i];
            if (r != null) u++;
        }
        assertEquals(0, u);
    }
}
