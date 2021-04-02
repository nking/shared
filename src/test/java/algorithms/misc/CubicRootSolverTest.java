package algorithms.misc;

import algorithms.util.FormatArray;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * @author nichole
 */
public class CubicRootSolverTest extends TestCase {
    
    public CubicRootSolverTest() {
    }
    
    public void test0() throws NotConvergedException {
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
        
        dCoeffs = CubicRootSolver.calcDepressedCubicCoefficients(coeffs);
        assertEquals(5, dCoeffs.length);
        a = dCoeffs[0];
        p = dCoeffs[3];
        q = dCoeffs[4];
        assertTrue(Math.abs(p - -6) < tol);
        assertTrue(Math.abs(q - -9) < tol);
        
        roots = CubicRootSolver.solveUsingDepressedCubic(p, q);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.mathemania.com/lesson/cardanos-formula-solving-cubic-equations/
        // x^3 - 15*x - 4 = 0.  
        // discriminant is < 0.   
        //solutions are approx 4, -2-sqrt(3), -2+sqrt(3)
        coeffs = new double[]{1, 0, -15, -4};
        expected = new double[]{4., -2 + Math.sqrt(3), -2 - Math.sqrt(3)};
        dCoeffs = CubicRootSolver.calcDepressedCubicCoefficients(coeffs);
        assertEquals(5, dCoeffs.length);
        a = dCoeffs[0];
        p = dCoeffs[3];
        q = dCoeffs[4];
        
        roots = CubicRootSolver.solveUsingDepressedCubic(p, q);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        
        
        //https://www.math.ucdavis.edu/~kkreith/tutorials/sample.lesson/cardano.html
        //x3 + x2 - 2 = 0  
        // answer: x=1
        coeffs = new double[]{1, 1, 0, -2};
        expected = new double[]{1};
        /*
        >>> numpy.roots([1, 1, 0, -2])
            array([-1.+1.j, -1.-1.j,  1.+0.j])
        */
        double[][] m = new double[3][3];
        m[0] = new double[]{-1., -0.,  2.};
        m[1] = new double[]{1.,  0.,  0.};
        m[2] = new double[]{0.,  1.,  0.};
        EVD evd = no.uib.cipr.matrix.EVD.factorize(new DenseMatrix(m));
        System.out.printf("real eigs=\n%s\n", FormatArray.toString(evd.getRealEigenvalues(), "%.3e"));
        System.out.printf("imag eigs=\n%s\n", FormatArray.toString(evd.getImaginaryEigenvalues(), "%.3e"));
        System.out.flush();
        for (i = 0; i < evd.getRealEigenvalues().length; ++i) {
            Complex c1 = new Complex(evd.getRealEigenvalues()[i], evd.getImaginaryEigenvalues()[i]);
            System.out.printf("   c=%.3e\n", c1.abs());
        }
        
        dCoeffs = CubicRootSolver.calcDepressedCubicCoefficients(coeffs);
        p = dCoeffs[3];
        q = dCoeffs[4];
        double absq = Math.abs(q);
        double t0 = (-2.*absq/q) * Math.sqrt(-p/3.) *
            Math.cosh((1./3.)*MiscMath0.acosh(((-3*absq)/(2.*p)) * Math.sqrt(-3/p)));
        
        //roots = CubicRootSolver.solve(coeffs);
        roots = CubicRootSolver.solveUsingGeneral(coeffs);
        assertEquals(expected.length, roots.length);
        for (i = 0; i < roots.length; ++i) {
            diff = Math.abs(expected[i] - roots[i]);
            assertTrue(diff < tol);
        }
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        // x^3 - 6x^2 + 11x - 6 = 0
        // 3 real distinct roots (discriminant < 0 in context of CubicSolver)
        // solns x=1,2,3
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        // x^3 - 5x^2 + 8x - 4 = 0
        // 3 real roots, but not distinct (2 are same)
        //     (discriminant == 0 in context of CubicSolver?)
        // solns are x=1,2
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 − 3x^2 + 3X - 1 =0
        // 3 real roots, but they're all the same
        // soln x=1
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        // x^3 + x^2 + x - 3 = 0
        // has 1 real solution 
        // (discriminant > 0 in context of CubicSolver)
        // real solution is x=1
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 - 5x^2 - 2x + 24 = 0
        // solutions x=-2,3,4
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 - 7x - 7 = 0
        // solutions x = -2, -1, 3
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 - - 4x^2 - 9x + 3 = 0
        // solutions x = -3, 3, 4
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 - 6x^2 - 6x - 7 = 0
        // one real solution 
        // solutions x = 7
        
        //https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-cubicequations-2009-1.pdf
        //x^3 + 3x^2 + 3x + 1 = 0
        // one real solution 
        // solutions x = -1
        
        
    }
}
