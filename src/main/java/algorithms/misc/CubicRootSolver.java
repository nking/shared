
package algorithms.misc;

import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * cubic root methods
 * @author nichole
 */
public class CubicRootSolver {
   
    /**
     * solve for t given coefficients p and q using the depressed cubic root t^3 + p*t + q = 0.
     * following  
     * https://www.mathemania.com/lesson/cardanos-formula-solving-cubic-equations/
     * https://www.wikiwand.com/en/Cubic_equation#/google_vignette
     * https://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
     * NOTE: only solving for the real, rational roots.
     * @param p coefficient p for the first order term.
     * @param q coefficient q for the zero order term.
     * @return 
     */
    public static double[] solveUsingDepressedCubic(double p, double q) {
        
        double eps = 1.e-5;
        
        double discr = Math.pow(q/2., 2) + Math.pow(p/3., 3);
        
        /*
        TODO: handle characteristic 3 modifications
        
        characteristic 3

        When the coefficients belong to a field of characteristic 3, 
        the results be modified because of the involved divisions by 3.
        The main tool for modification is the fact that a multiple root is a 
        common root of the polynomial and its formal derivative. 
        In this characteristic, if the derivative is not a constant, 
        it has a single root, being linear in characteristic 3. 
        This allows computing the multiple root, and the third root can be 
        deduced from the sum of the roots, which is provided by Vietas formulas.
        In characteristic characteristic 3, the formula for a triple root involves a cube root.
        */           
                 
        if (Math.abs(discr) < eps) {
            // discriminant is 0.  
            //there are 3 real roots, not necessarily distinct (has a multiple root)
            if (Math.abs(p) < eps) {
                return new double[]{0};
            }
                       
            // Vieta's formula, Viete's formula
            double t1 = 3.*q/p;
            double t2 = -t1/2.;
            return new double[]{t1, t2};
        } else if (discr > 0) {
            // cardano's formula.  valid for discriminant == 0 and > 0
            // one real root and 2 non-real
      
            double pt2 = Math.sqrt((q*q/4.) + p*p*p/27.);
            double pt1 = -q/2.;
            double t = Math.pow(pt1 + pt2, 1./3.) + Math.pow(pt1 - pt2, 1./3.);
            return new double[]{t};
        } else {
            // discriminant < 0.  3 real distinct roots
            // trigonometric solution
            //https://www.wikiwand.com/en/Cubic_equation#/google_vignette
            assert(p < 0);
            
            double pt1 = 2.*Math.sqrt(-p/3);
            double pt2 = Math.acos(((3.*q)/(2.*p)) * Math.sqrt(-3./p));
            double tp3 = 2.*Math.PI/3.;
            double[] t = new double[3];
            for (int i = 0; i < 3; ++i) {
                t[i] = pt1 * Math.cos((1./3.)*pt2 - (tp3*i)); ;
            }
            return t;
            
            /*
            //for only 1 real root and p < 0, can use the hyperbolic solution:
            double t0;
            if (p < 0) {
                double absq = Math.abs(q);
                t0 = (-2.*absq/q) * Math.sqrt(-p/3.) *
                    Math.cosh((1./3.)*acosh(((-3*absq)/(2.*p)) * Math.sqrt(-3/p)));
            } else {
                t0 = -2.* Math.sqrt(p/3.) *
                    Math.sinh((1./3.)*asinh(((3*q)/(2.*p)) * Math.sqrt(3/p)));
            }
            */
        }
    }
    
    /**
     * calculates coefficients for the reduced form od the cubic equation
     * called the depressed cubic.
     * @param cubicCoeff
     * @return a, b, c, p, q.  returns null if cubicCoeff[0] = 0.
     */
    static double[] calcDepressedCubicCoefficients(double[] cubicCoeff) {
        if (cubicCoeff.length != 4) {
            throw new IllegalArgumentException("cubicCoeff must be length 4");
        }
        
        double a3 = cubicCoeff[0];
        double a2 = cubicCoeff[1];
        double a1 = cubicCoeff[2];
        double a0 = cubicCoeff[3];
        
        double tol = 1e-5;
        
        if (Math.abs(a3) < tol) {
            return null;
        }
        
        double a = a2/a3;
        double b = a1/a3;
        double c = a0/a3;
                
        // p = b - a^2/3
        // q = (2*a^3/27) - (a*b/3) + c
        
        double p = b - a*a/3;
        
        double q = (2.*a*a*a/27.) - (a*b/3.) + c;
        
        return new double[]{a, b, c, p, q};
    }
    
    /**
     * solve for the roots of the equation a*x^3 + b*x^2 + c*x + d = 0
     * using the depressed cubic, given coefficients a, b, c, d. returns x.
     * following 
     * https://www.mathemania.com/lesson/cardanos-formula-solving-cubic-equations/
     * https://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
     * NOTE: only solving for the real roots.
     * @param coeffs an array holding a, b, c , d for the equation 
     * a*x^3 + b*x^2 + c*x + d = 0.
     * @return x as the roots of the cubic equation a*x^3 + b*x^2 + c*x + d = 0
     */
    public static double[] solve(double[] coeffs) {
        if (coeffs.length != 4) {
            throw new IllegalArgumentException("coeffs must be length 4");
        }
        
        double[] dCoeffs = calcDepressedCubicCoefficients(coeffs);
                
        double a = dCoeffs[0];
        double b = dCoeffs[1];
        double c = dCoeffs[2];
        double p = dCoeffs[3];
        double q = dCoeffs[4];
        
        // if discriminant is 0 and characteristic is 2 or 3, need to use
        //    the general solution.
        
        double t[] = solveUsingDepressedCubic(p, q);
        
        double[] x = new double[t.length];
        for (int i = 0; i < x.length; ++i) {
            x[i] = t[i] - b/3.*a;
        }
        
        return x;
    }
    
    /**
     * solve for the roots of the equation a*x^3 + b*x^2 + c*x + d = 0
     * using the general formula, given coefficients a, b, c, d. returns x.
     * following
     * https://www.wikiwand.com/en/Cubic_equation#/General_cubic_formula
     * @param coeffs an array holding a, b, c , d for the equation 
     * a*x^3 + b*x^2 + c*x + d = 0.
     * @return x as the roots of the cubic equation a*x^3 + b*x^2 + c*x + d = 0
     */
    public static double[] solveUsingGeneral(double[] coeffs) throws NotConvergedException {
        
        //https://www.wikiwand.com/en/Cubic_equation#/General_cubic_formula
        
        if (coeffs.length != 4) {
            throw new IllegalArgumentException("coeffs must be length 4");
        }
                        
        double a = coeffs[0];
        double b = coeffs[1];
        double c = coeffs[2];
        double d = coeffs[3];
        
        double tol = 1e-5;
        
        /*
        TODO: handle characteristic 2 and 3 modifications
                
        characteristic 2 and 3
        
        When the coefficients belong to a field of characteristic 2 or 3, 
        the results be modified because of the involved divisions by 2 and 3.
        The main tool for modification is the fact that a multiple root is a 
        common root of the polynomial and its formal derivative. 
        In these characteristics, if the derivative is not a constant, 
        it has a single root, being linear in characteristic 3, 
        or the square of a linear polynomial in characteristic 2. 
        This allows computing the multiple root, and the third root can be 
        deduced from the sum of the roots, which is provided by Vietas formulas.
        In characteristic 2, the formula for a double root involves a square root, 
        and, in characteristic 3, the formula for a triple root involves a cube root.
        */              
        
        double discr = 18.*a*b*c*d - 4*b*b*b*d + b*b*c*c - 4*a*c*c*c - 27.*a*a*d*d;
        
        if (Math.abs(discr) < tol) {
            // discriminant is 0.  
            //there are 3 real roots, not necessarily distinct (has a multiple root)
            
            if (Math.abs(b*b - 3*a*c) < 0) {
                // either b^2 = 3ac, triple root but all the same
                
                return new double[]{-b/(3.*a)};
                
            } else {
                // or b^2 != 3ac, double root
                
                // a double root
                double x2 = (9.*a*d - b*c)/(2.*(b*b - 3.*a*c));
                
                // a simple root
                double x1 = (4.*a*b*c - 9.*a*a*d - b*b*b)/(a*(b*b - 3.*a*c));
                
                return new double[]{x2, x1};
            }
            
        } else {
            
            /*
            if (discr < 0) {
                // 3 distinct real roots
                //use the depressed cubic which is the trigonometric solution
                //then back substitution here
                
            } else {
                // discriminant > 0
                // 1 real root
            }  
            */
            double[] roots = PolynomialRootSolver.realRoots(coeffs);
            
            return roots;
        }
        
    }
}
