/*
The author of this code is Wilco Oelen and he offers it
freely without copyright, but asks that his pages are referenced
as the source if used.
He has a webpage with information on the polynomial software
he ported and more modern versions which require jini bindings:
https://woelen.homescience.net/science/math/exps/polynomials/
https://woelen.homescience.net/science/math/exps/polynomials/software.html
The code here is from the Java port of RPoly, CPoly and MPSolve 1996 algorithms:
https://woelen.homescience.net/science/math/exps/polynomials/software/polsolve.tgz

MPSolve is an implementation of the algorithms of Bini & Fiorentino 2000,
"Design, analysis, and implementation of a multiprecision polynomial rootfinder".
...Counting, isolating and approximating all roots in a given set S are the main 
goals that the algorithm provides. Automatic determination of multiplicities 
and the detection of real or imaginary roots can be selected as well. 
Polynomials having coef- ficients with a bounded precision may be processed too.
...
The algorithm is based on an adaptive strategy which automatically exploits 
any specific feature of the input polynomial, like its sparsity or the 
conditioning of its roots, in order to speed up the computation. 
... 

The resulting algorithm MPSolve, implemented in C, can deal with polynomials 
having real or complex coefficients with integer, rational or floating point 
real and imaginary parts with any number of digits. The algorithm can also 
process polynomials with approximately known coefficients.
*/
package thirdparty.net.oelen.polsolve.pzeros;

import thirdparty.net.oelen.polarith.Complex;
import thirdparty.net.oelen.polarith.DoubleComplex;
import thirdparty.net.oelen.polarith.DoubleDouble;


// This class is the class that has public methods and can be used by
// user programs. The other classes in this package (PZerosD, PZerosDD,
// and Convex) are not intended for use by user programs, they are 
// helper classes for PZeros.
// Below follows documentation of all public methods, which are provided
// by the PZeros class.

public strictfp class PZeros {
    private static int MAX_ITERATIONS = 10000;
    
    private int degree;
    private DoubleDouble[] coefsReal;
    private DoubleComplex[] coefsCplx;
    
    
    /**
     * Constructor that takes an array of real coefficients. Coefficients are
     * supplied in order of increasing power. Example: p(x) = A + B*x + C*x^2 +
     * D*x^3 + E*x^4 gives coefficients[] = {A, B, C, D, E}. The degree of the
     * polynomial is equal to the length of the supplied array minus 1.
     *
     * @param coef An array of real coefficients in order of increasing power.
     */
    public PZeros(double[] coef) {
        degree = coef.length - 1;
        coefsReal = new DoubleDouble[coef.length];
        for (int i=0; i<coef.length; i++) {
            coefsReal[i] = new DoubleDouble(coef[i]);
        }
        coefsCplx = null;
    }
    
    
    
    
    
    /**
     * Constructor that takes an array of real coefficients. Coefficients are
     * supplied in order of increasing power. Example: p(x) = A + B*x + C*x^2 +
     * D*x^3 + E*x^4 gives coefficients[] = {A, B, C, D, E}. The degree of the
     * polynomial is equal to the length of the supplied array minus 1.
     *
     * @param coef An array of real coefficients in order of increasing power.
     */
    public PZeros(String[] coef) {
        degree = coef.length - 1;
        coefsReal = new DoubleDouble[coef.length];
        for (int i=0; i<coef.length; i++) {
            coefsReal[i] = new DoubleDouble(coef[i]);
        }
        coefsCplx = null;
    }
    
    
    
    
    
    /**
     * Constructor that takes an array of real coefficients. Coefficients are
     * supplied in order of increasing power. Example: p(x) = A + B*x + C*x^2 +
     * D*x^3 + E*x^4 gives coefficients[] = {A, B, C, D, E}. The degree of the
     * polynomial is equal to the length of the supplied array minus 1.
     *
     * @param coef An array of real coefficients in order of increasing power.
     */
    public PZeros(DoubleDouble[] coef) {
        degree = coef.length - 1;
        coefsReal = new DoubleDouble[coef.length];
        for (int i=0; i<coef.length; i++) {
            coefsReal[i] = coef[i];
        }
        coefsCplx = null;
    }

    
    
    
    /**
     * Constructor that takes an array with the real part of coefficients and 
     * an array with the imaginary part of the coefficients. Coefficients are
     * supplied in order of increasing power. Example: p(x) = A + B*x + C*x^2 +
     * D*x^3 + E*x^4 gives coefficients[] = {A, B, C, D, E}.
     * The degree of the polynomial is determined by the length of the longest
     * supplied array. If the longest supplied array has N elements, then the
     * degree of the polynomial equals N-1. The shorter array is extended with
     * zero values for the higher powers. E.g. a polynomial with arrays
     * {1,2,3} and {11,22,33,44,55} has degree 4 and can be written as
     * (1+11i) + (2+22i)*x + (3+33i)*x^2 + 44i*x^3 + 55i*x^4
     *
     * @param coef_re An array containing the real part of the coefficients 
     * in order of increasing power. If the supplied array equals null, then
     * the degree is determined by the length of the other array and the real
     * part of all coefficients equals 0 in that case.
     * @param coef_im An array containing the imaginary part of the coefficients
     * in order of increasing power. If the supplied array equals null, then
     * the degree is determined by the length of the other array and the imaginary
     * part of all coefficients equals 0 in that case.
     */
    public PZeros(double[] coef_re, double[] coef_im) {
        if ((coef_re == null || coef_re.length==0) && 
            (coef_im == null || coef_im.length==0)) {
            throw new RuntimeException("Construction of PZeros with empty coefficient set.");
        }
        if (coef_im == null || coef_im.length==0) {
            degree = coef_re.length - 1;
            coefsReal = new DoubleDouble[coef_re.length];
            for (int i=0; i<coef_re.length; i++) {
                coefsReal[i] = new DoubleDouble(coef_re[i]);
            }
            coefsCplx = null;
            return;
        }
        
        if (coef_re == null) {
            coef_re = new double[0];
        }
        
        int degRe = coef_re.length - 1;
        int degIm = coef_im.length - 1;
        
        if (degRe < degIm) {
            degree = degIm;
            coefsCplx = new DoubleComplex[coef_im.length];
            int i;
            for (i=0; i<=degRe; i++) {
                coefsCplx[i] = new DoubleComplex(coef_re[i], coef_im[i]);
            }
            for (; i<coef_im.length; i++) {
                coefsCplx[i] = new DoubleComplex(DoubleDouble.ZERO, new DoubleDouble(coef_im[i]));
            }
        }
        else {
            degree = degRe;
            coefsCplx = new DoubleComplex[coef_re.length];
            int i;
            for (i=0; i<=degIm; i++) {
                coefsCplx[i] = new DoubleComplex(coef_re[i], coef_im[i]);
            }
            for (; i<coef_re.length; i++) {
                coefsCplx[i] = new DoubleComplex(coef_re[i]);
            }
        }
        
        coefsReal = null;
    }

    
    
    
    /**
     * Constructor that takes an array with the real part of coefficients and 
     * an array with the imaginary part of the coefficients. Coefficients are
     * supplied in order of increasing power. Example: p(x) = A + B*x + C*x^2 +
     * D*x^3 + E*x^4 can be created with coef[] = {A, B, C, D, E}.
     * The degree of the polynomial is determined by the length of the longest
     * supplied array. If the longest supplied array has N elements, then the
     * degree of the polynomial equals N-1. The shorter array is extended with
     * zero values for the higher powers. E.g. a polynomial with arrays
     * {1,2,3} and {11,22,33,44,55} has degree 4 and can be written as
     * (1+11i) + (2+22i)*x + (3+33i)*x^2 + 44i*x^3 + 55i*x^4
     *
     * @param coef_re An array containing the real part of the coefficients 
     * in order of increasing power. If the supplied array equals null, then
     * the degree is determined by the length of the other array and the real
     * part of all coefficients equals 0 in that case.
     * @param coef_im An array containing the imaginary part of the coefficients
     * in order of increasing power. If the supplied array equals null, then
     * the degree is determined by the length of the other array and the imaginary
     * part of all coefficients equals 0 in that case.
     */
    public PZeros(String[] coef_re, String[] coef_im) {
        if ((coef_re == null || coef_re.length==0) && 
            (coef_im == null || coef_im.length==0)) {
            throw new RuntimeException("Construction of PZeros with empty coefficient set.");
        }
        if (coef_im == null || coef_im.length==0) {
            degree = coef_re.length - 1;
            coefsReal = new DoubleDouble[coef_re.length];
            for (int i=0; i<coef_re.length; i++) {
                coefsReal[i] = new DoubleDouble(coef_re[i]);
            }
            coefsCplx = null;
            return;
        }
        
        if (coef_re == null) {
            coef_re = new String[0];
        }
        
        int degRe = coef_re.length - 1;
        int degIm = coef_im.length - 1;
        
        if (degRe < degIm) {
            degree = degIm;
            coefsCplx = new DoubleComplex[coef_im.length];
            int i;
            for (i=0; i<=degRe; i++) {
                coefsCplx[i] = new DoubleComplex(coef_re[i], coef_im[i]);
            }
            for (; i<coef_im.length; i++) {
                coefsCplx[i] = new DoubleComplex(DoubleDouble.ZERO, new DoubleDouble(coef_im[i]));
            }
        }
        else {
            degree = degRe;
            coefsCplx = new DoubleComplex[coef_re.length];
            int i;
            for (i=0; i<=degIm; i++) {
                coefsCplx[i] = new DoubleComplex(coef_re[i], coef_im[i]);
            }
            for (; i<coef_re.length; i++) {
                coefsCplx[i] = new DoubleComplex(coef_re[i]);
            }
        }
        
        coefsReal = null;
    }

    
    
    
    /**
     * Constructor that takes an array with the real part of coefficients and 
     * an array with the imaginary part of the coefficients. Coefficients are
     * supplied in order of increasing power. Example: p(x) = A + B*x + C*x^2 +
     * D*x^3 + E*x^4 can be created with coef[] = {A, B, C, D, E}.
     * The degree of the polynomial is determined by the length of the longest
     * supplied array. If the longest supplied array has N elements, then the
     * degree of the polynomial equals N-1. The shorter array is extended with
     * zero values for the higher powers. E.g. a polynomial with arrays
     * {1,2,3} and {11,22,33,44,55} has degree 4 and can be written as
     * (1+11i) + (2+22i)*x + (3+33i)*x^2 + 44i*x^3 + 55i*x^4
     *
     * @param coef_re An array containing the real part of the coefficients 
     * in order of increasing power. If the supplied array equals null, then
     * the degree is determined by the length of the other array and the real
     * part of all coefficients equals 0 in that case.
     * @param coef_im An array containing the imaginary part of the coefficients
     * in order of increasing power. If the supplied array equals null, then
     * the degree is determined by the length of the other array and the imaginary
     * part of all coefficients equals 0 in that case.
     */
    public PZeros(DoubleDouble[] coef_re, DoubleDouble[] coef_im) {
        if ((coef_re == null || coef_re.length==0) && 
            (coef_im == null || coef_im.length==0)) {
            throw new RuntimeException("Construction of PZeros with empty coefficient set.");
        }
        if (coef_im == null || coef_im.length==0) {
            degree = coef_re.length - 1;
            coefsReal = new DoubleDouble[coef_re.length];
            for (int i=0; i<coef_re.length; i++) {
                coefsReal[i] = coef_re[i];
            }
            coefsCplx = null;
            return;
        }
        
        if (coef_re == null) {
            coef_re = new DoubleDouble[0];
        }
        
        int degRe = coef_re.length - 1;
        int degIm = coef_im.length - 1;
        
        if (degRe < degIm) {
            degree = degIm;
            coefsCplx = new DoubleComplex[coef_im.length];
            int i;
            for (i=0; i<=degRe; i++) {
                coefsCplx[i] = new DoubleComplex(coef_re[i], coef_im[i]);
            }
            for (; i<coef_im.length; i++) {
                coefsCplx[i] = new DoubleComplex(DoubleDouble.ZERO, coef_im[i]);
            }
        }
        else {
            degree = degRe;
            coefsCplx = new DoubleComplex[coef_re.length];
            int i;
            for (i=0; i<=degIm; i++) {
                coefsCplx[i] = new DoubleComplex(coef_re[i], coef_im[i]);
            }
            for (; i<coef_re.length; i++) {
                coefsCplx[i] = new DoubleComplex(coef_re[i]);
            }
        }
        
        coefsReal = null;
    }
    
    
    
    
    
    /**
     * Constructor that takes an array of complex coefficients. Coefficients are
     * supplied in order of increasing power. Example: p(x) = A + B*x + C*x^2 +
     * D*x^3 + E*x^4 can be created with coef[] = {A, B, C, D, E}. The degree of the
     * polynomial is equal to the length of the supplied array minus 1.
     *
     * @param coef An array of complex coefficients in order of increasing power.
     */
    public PZeros(DoubleComplex[] coef) {
        degree = coef.length - 1;
        coefsCplx = new DoubleComplex[coef.length];
        for (int i=0; i<coef.length; i++) {
            coefsCplx[i] = coef[i];
        }
        coefsReal = null;
    }
    
    
    
    
    
    /**
     * Constructor that takes an array of complex coefficients. Coefficients are
     * supplied in order of increasing power. Example: p(x) = A + B*x + C*x^2 +
     * D*x^3 + E*x^4 can be created with coef[] = {A, B, C, D, E}. The degree of the
     * polynomial is equal to the length of the supplied array minus 1.
     *
     * @param coef An array of complex coefficients in order of increasing power.
     */
    public PZeros(Complex[] coef) {
        degree = coef.length - 1;
        coefsCplx = new DoubleComplex[coef.length];
        for (int i=0; i<coef.length; i++) {
            coefsCplx[i] = new DoubleComplex(coef[i]);
        }
        coefsReal = null;
    }

    
    
    
    /**
     * Returns the degree of the polynomial.
     *
     * @return The degree of the polynomial.
     */
    public int degree() {
        return degree;
    }
    
    
    
    
    
    
    /**
     * This method computes the roots of the polynomial and stores the roots
     * in preallocated arrays, which are passed as arguments. DoubleDouble
     * 105-bit precision is used for the calculations and the results also are
     * returned at DoubleDouble precision.
     * @param root Array, in which the zeros will be stored after computation of 
     * the zeros. This array must have a length of at least N elements, where N 
     * is the degree of the polynomial.
     * @param radius Array, which gives an indication of the accuracy of the found
     * roots. For each root, a radius is returned. The root is assured to be in the
     * disk with the corresponding radius, centered around the returned root.
     * @param err Array, which specifies whether the root and its corresponding
     * radius of accuracy could be determined correctly. If err[i] is true, then
     * the program did not converge for root[i].
     * @return Returns the degree of the polynomial if the computation succeeds,
     * and returns a value less than the degree of the polynomial if an error
     * occurs (e.g. convergence failure). When a value less than the degree of
     * the polynomial is returned, then only part (or none) of the roots could
     * be determined. If a negative value is returned, then the supplied input
     * is not correct:
     *   -1: Leading coefficient equals 0.
     *   -2: Coefficient for x^0 (constant coefficient) equals 0.
     *   -3: Ratio of smallest coefficient magnitude and largest coefficient
     *       magnitude is too large and will lead to underflow/overflow.
     */
    public int solve(DoubleComplex[] root, double[] radius, boolean[] err) {
        int status;
        int[] iter = new int[] {0};
        if (coefsReal != null) {
            status = pzeros(degree, coefsReal, MAX_ITERATIONS, root, radius, err, iter);
        }
        else {
            status = pzeros(degree, coefsCplx, MAX_ITERATIONS, root, radius, err, iter);
        }
        return status;
    }
    
    
    
    
    
    
    /**
     * This method computes the roots of the polynomial and stores the roots
     * in preallocated arrays, which are passed as arguments. Standard 53-bits
     * precision is used for the calculations and the results also are returned as
     * standard 53-bits precision numbers.
     * @param root Array, in which the zeros will be stored after computation of 
     * the zeros. This array must have a length of at least N elements, where N 
     * is the degree of the polynomial.
     * @param radius Array, which gives an indication of the accuracy of the found
     * roots. For each root, a radius is returned. The root is assured to be in the
     * disk with the corresponding radius, centered around the returned root.
     * @param err Array, which specifies whether the root and its corresponding
     * radius of accuracy could be determined correctly. If err[i] is true, then
     * the program did not converge for root[i].
     * @return Returns the degree of the polynomial if the computation succeeds,
     * and returns a value less than the degree of the polynomial if an error
     * occurs (e.g. convergence failure). When a value less than the degree of
     * the polynomial is returned, then only part (or none) of the roots could
     * be determined. If a negative value is returned, then the supplied input
     * is not correct:
     *   -1: Leading coefficient equals 0.
     *   -2: Coefficient for x^0 (constant coefficient) equals 0.
     *   -3: Ratio of smallest coefficient magnitude and largest coefficient
     *       magnitude is too large and will lead to underflow/overflow.
     */
    public int solve(Complex[] root, double[] radius, boolean[] err) {
        int status;
        int[] iter = new int[] {0};
        if (coefsReal != null) {
            status = pzeros(degree, DoubleDouble.toDouble(coefsReal), MAX_ITERATIONS, root, radius, err, iter);
        }
        else {
            status = pzeros(degree, DoubleComplex.toComplex(coefsCplx), MAX_ITERATIONS, root, radius, err, iter);
        }
        return status;
    }
    
    
    
    
    
    
    
    /***********************************************************************/
    /***********************************************************************/
    /********* Below follows the port of the original Fortran  *************/
    /********* program. It is all private to this module!     *************/
    /********* The code above is a wrapper for easy usage. *************/
    /***********************************************************************/
    /***********************************************************************/
    /***********************************************************************/
    
    
    
    // Driver code, which uses PolSolveD for normal double precision (53 bits),
    // and uses PolSolveDD for DoubleDouble precision (104 bits). This code
    // contains driver methods for solving real and complex polynomial equations.
    //
    // The driver methods for DoubleDouble precision use a special strategy. They
    // first use the normal precision code and if this cannot improve the solutions
    // anymore, or cannot find all solutions, then it uses the DoubleDouble code
    // to improve the situation. This strategy makes the use of DoubleDouble code
    // approximately 3 times as fast, compared with sole use of DoubleDouble. For
    // the great majority of polynomials, this strategy assures that the initial
    // search and isolation of roots is done at normal precision, while the slower
    // high precision code only is used at the final stages of refinement.
    
    private static final double BIG = Double.MAX_VALUE;     // near overflow of double
    private static final double EPS_D = 1.0 / (1l<<52);     // 52 bits of precision
    private static final double EPS_DD = EPS_D*EPS_D;       // 104 bits of precision
    

    
    // Driver method for polynomial with real coefficients at 105 bit precision.
    private static int pzeros(int n, DoubleDouble[] a, int itmax,
                             DoubleComplex[] root, double[] radius, boolean[] err, int[] iter) {
        double[] apoly = new double[n + 1];
        double[] apolyr = new double[n + 1];
        
        // If we have a linear equation, simply call the high precision code
        // and let it compute the zero at high precision.
        if (n == 1) {
            int status = PZerosDD.polzeros(n, a, n*EPS_DD, BIG, 1,
                    root, radius, err, iter,
                    true, apoly, apolyr, true);
            return status == 0 ? 1 : status;
        }
        
        // First call the lower precision code and prepare coefficient arrays
        // and root arrays for the lower precision numbers.
        double[] a_d = new double[a.length];
        for (int i=0; i<a.length; i++) {
            a_d[i] = a[i].doubleValue();
        }
        Complex[] root_d = new Complex[root.length];
        int status = PZerosD.polzeros(n, a_d, n*EPS_D, BIG, 3*itmax/5,
                root_d, radius, err, iter,
                true, apoly, apolyr, false);
        if (status != 0 && status != -4) {
            return status;
        }
        
        // Convert the found roots to the higher precision format and prepare
        // a boolean status vector for the high precision code. With these
        // high quality initial estimates continue the process of root finding
        // at high precision.
        for (int i=0; i<root.length; i++) {
            root[i] = new DoubleComplex(root_d[i].real(), root_d[i].imag());
            err[i] = (radius[i] >= 0.0);
        }
        status = PZerosDD.polzeros(n, a, n*EPS_DD, BIG, itmax,
                      root, radius, err, iter,
                      false, apoly, apolyr, true);
        return transformStatus(status, err);
    }
    
    
    
    
    // Driver method for polynomial with complex coefficients at 104 bit precision.
    private static int pzeros(int n, DoubleComplex[] a, int itmax,
                             DoubleComplex[] root, double[] radius, boolean[] err, int[] iter) {
        double[] apoly = new double[n + 1];
        double[] apolyr = new double[n + 1];
        
        // If we have a linear equation, simply call the high precision code
        // and let it compute the zero at high precision.
        if (n == 1) {
            int status = PZerosDD.polzeros(n, a, n*EPS_DD, BIG, 1,
                    root, radius, err, iter,
                    true, apoly, apolyr, true);
            return status == 0 ? 1 : status;
        }
        
        // First call the lower precision code and prepare coefficient arrays
        // and root arrays for the lower precision numbers.
        Complex[] a_d = new Complex[a.length];
        for (int i=0; i<a.length; i++) {
            a_d[i] = new Complex(a[i].real().doubleValue(), a[i].imag().doubleValue());
        }
        Complex[] root_d = new Complex[root.length];
        int status = PZerosD.polzeros(n, a_d, n*EPS_D, BIG, 3*itmax/5,
                root_d, radius, err, iter,
                true, apoly, apolyr, false);
        if (status != 0 && status != -4) {
            return status;
        }
        
        
        // Convert the found roots to the higher precision format and prepare
        // a boolean status vector for the high precision code. With these
        // high quality initial estimates continue the process of root finding
        // at high precision.
        for (int i=0; i<root.length; i++) {
            root[i] = new DoubleComplex(root_d[i].real(), root_d[i].imag());
            err[i] = (radius[i] >= 0.0);
        }
        status = PZerosDD.polzeros(n, a, n*EPS_DD, BIG, itmax,
                    root, radius, err, iter,
                    false, apoly, apolyr, true);
        return transformStatus(status, err);
    }
    
    
    

    // Driver method for polynomial with real coefficients at 53 bit precision.
    private static int pzeros(int n, double[] a, int itmax,
                             Complex[] root, double[] radius, boolean[] err, int[] iter) {
        double[] apoly = new double[n + 1];
        double[] apolyr = new double[n + 1];
        int status = PZerosD.polzeros(n, a, EPS_D, BIG, itmax,
                root, radius, err, iter,
                true, apoly, apolyr, true);
        return transformStatus(status, err);
    }
    
    
    
    
    // Driver method for polynomial with complex coefficients at 53 bit precision.
    private static int pzeros(int n, Complex[] a, int itmax,
                             Complex[] root, double[] radius, boolean[] err, int[] iter) {
        double[] apoly = new double[n + 1];
        double[] apolyr = new double[n + 1];
        int status = PZerosD.polzeros(n, a, EPS_D, BIG, itmax,
                root, radius, err, iter,
                true, apoly, apolyr, true);
        return transformStatus(status, err);
    }
    
    
    
    
    
    private static int transformStatus(int status, boolean[] err) {
        if (status == 0 || status == -4) {
            status = 0;
            for (boolean e : err) {
                if (!e) status++;
            }
        }
        return status;
    }
}
