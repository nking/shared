package algorithms.statistics;

import algorithms.misc.Complex;
import algorithms.util.FormatArray;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * @author nichole
 */
public class CubicRootSolverTest extends TestCase {
    
    public static double eps = 1e-4;
    
    public CubicRootSolverTest() {
    }
    
    public void test0() throws NotConvergedException, Exception {
        double[] coeffs, dCoeffs;
        double[] expected;
        double[] roots;
        double a, b, c, p, q;
        double tol = 1e-3;
        double diff;
        int i;
        
        //https://www.mathemania.com/lesson/cardanos-formula-solving-cubic-equations/
        // x^3 - 6*x - 9 = 0
        //  discrimant is > 0.  
        // solution for depressed cubic is = 2+1 = 3
        coeffs = new double[]{1, 0, -6, -9};
        expected = new double[]{3};
        
        dCoeffs = algorithms.statistics.CubicRootSolver.calcDepressedCubicCoefficients(coeffs);
        assertEquals(5, dCoeffs.length);
        a = dCoeffs[0];
        p = dCoeffs[3];
        q = dCoeffs[4];
        assertTrue(Math.abs(p - -6) < tol);
        assertTrue(Math.abs(q - -9) < tol);
        
        roots = algorithms.statistics.CubicRootSolver.solveUsingDepressedCubic(p, q);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        Complex[] rootsC = algorithms.statistics.CubicRootSolver.solveUsingDepressedCubic0(p, q);
        assertTrue(Math.abs(expected[0] - rootsC[2].re()) < tol);
        assertTrue(Math.abs(0 - rootsC[2].im()) < tol);
        roots = algorithms.statistics.CubicRootSolver.realNonZeroOnly(rootsC);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        
        Complex[] roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        //System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        double[] dRoots2 = algorithms.statistics.PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[i] - dRoots2[i]);
            assertTrue(diff < tol);
        }
        
        int ph = 1;
        
        //https://www.mathemania.com/lesson/cardanos-formula-solving-cubic-equations/
        // x^3 - 15*x - 4 = 0.  
        // discriminant is < 0.   
        //solutions are approx 4, -2-sqrt(3), -2+sqrt(3)
        coeffs = new double[]{1, 0, -15, -4};
        expected = new double[]{4., -2 + Math.sqrt(3), -2 - Math.sqrt(3)};
        dCoeffs = algorithms.statistics.CubicRootSolver.calcDepressedCubicCoefficients(coeffs);
        assertEquals(5, dCoeffs.length);
        a = dCoeffs[0];
        p = dCoeffs[3];
        q = dCoeffs[4];
        
        roots = algorithms.statistics.CubicRootSolver.solveUsingDepressedCubic(p, q);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        dRoots2 = algorithms.statistics.PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        
        System.out.printf("mpsolve real=%s\n", FormatArray.toString(dRoots2, "%.3e"));
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[3-i-1] - dRoots2[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.math.ucdavis.edu/~kkreith/tutorials/sample.lesson/cardano.html
        //x3 + x2 - 2 = 0  
        // answer: x=1
        coeffs = new double[]{1, 1, 0, -2};
        expected = new double[]{1};
  
        //roots = CubicRootSolver.solve(coeffs);
        roots = algorithms.statistics.CubicRootSolver.solveUsingGeneral(coeffs);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        dRoots2 = algorithms.statistics.PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[i] - dRoots2[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        // x^3 - 6x^2 + 11x - 6 = 0
        // 3 real distinct roots (discriminant < 0 in context of CubicSolver)
        // solns x=1,2,3
        coeffs = new double[]{1, -6, 11, -6};
        expected = new double[]{3, 2, 1};
  
        roots = algorithms.statistics.CubicRootSolver.solve(coeffs);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        dRoots2 = algorithms.statistics.PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[3-i-1] - dRoots2[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        // x^3 - 5x^2 + 8x - 4 = 0
        // 3 real roots, but not distinct (2 are same)
        //     (discriminant == 0 in context of CubicSolver?)
        // solns are x=1,2
        coeffs = new double[]{1, -5, 8, -4};
        expected = new double[]{1, 2};
  
        roots = algorithms.statistics.CubicRootSolver.solve(coeffs);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        dRoots2 = algorithms.statistics.PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        System.out.printf("mpsolve real=%s\n", FormatArray.toString(dRoots2, "%.3e"));
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[i] - dRoots2[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 − 3x^2 + 3X - 1 =0
        // 3 real roots, but they're all the same
        // soln x=1
        coeffs = new double[]{1, -3, 3, -1};
        expected = new double[]{1};
  
        roots = algorithms.statistics.CubicRootSolver.solve(coeffs);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.7e"));
        dRoots2 = algorithms.statistics.PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        System.out.printf("mpsolve real=%s\n", FormatArray.toString(dRoots2, "%.7e"));
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[i] - dRoots2[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        // x^3 + x^2 + x - 3 = 0
        // has 1 real solution 
        // (discriminant > 0 in context of CubicSolver)
        // real solution is x=1
        coeffs = new double[]{1, 1, 1, -3};
        expected = new double[]{1};
  
        roots = algorithms.statistics.CubicRootSolver.solve(coeffs);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        dRoots2 = algorithms.statistics.PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[i] - dRoots2[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 - 5x^2 - 2x + 24 = 0
        // solutions x=-2,3,4
        coeffs = new double[]{1, -5, -2, 24};
        expected = new double[]{4,3,-2};
  
        roots = algorithms.statistics.CubicRootSolver.solve(coeffs);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        dRoots2 = algorithms.statistics.PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[3-i-1] - dRoots2[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 - 7x - 6 = 0
        // solutions x = -2, -1, 3
        coeffs = new double[]{1, 0, -7, -6};
        expected = new double[]{3, -1, -2};
        roots = algorithms.statistics.CubicRootSolver.solve(coeffs);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        dRoots2 = algorithms.statistics.PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[3-i-1] - dRoots2[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 - 4x^2 - 9x + 36 = 0
        // solutions x = -3, 3, 4
        coeffs = new double[]{1, -4, -9, 36};
        expected = new double[]{4,3,-3};
  
        roots = algorithms.statistics.CubicRootSolver.solve(coeffs);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        dRoots2 = algorithms.statistics.PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[3-i-1] - dRoots2[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 - 6x^2 - 6x - 7 = 0
        // one real solution 
        // solutions x = 7
        coeffs = new double[]{1, -6, -6, -7};
        expected = new double[]{7};
  
        roots = algorithms.statistics.CubicRootSolver.solve(coeffs);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        dRoots2 = algorithms.statistics.PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[i] - dRoots2[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 + 3x^2 + 3x + 1 = 0
        // one real solution 
        // solutions x = -1
        coeffs = new double[]{1, 3, 3, 1};
        expected = new double[]{-1};
  
        roots = CubicRootSolver.solve(coeffs);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        roots2 = algorithms.statistics.PolynomialRootSolver.solveUsingMPSolve(coeffs);
        System.out.printf("mpsolve roots=%s\n", FormatArray.toString(roots2, "%.3e"));
        dRoots2 = PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
        assertEquals(expected.length, dRoots2.length);
        for (i = 0; i < dRoots2.length; ++i) {
            diff = Math.abs(expected[i] - dRoots2[i]);
            assertTrue(diff < tol);
        }
    }
}
