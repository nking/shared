/*
 *   CPoly  -- A class that represents a polynomial equation.
 *
 *   Copyright (C) 2000-2004  by Joseph A. Huwaldt   
 *                 2014-2016  by Wilco Oelen
 *   All rights reserved.
 *   
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *   
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *   Or visit:  http://www.gnu.org/licenses/lgpl.html
 **/
package thirdparty.net.oelen.polsolve.jt;


/**
 * <p>
 * Represents a polynomial equation of the form 
 *           p(x) = A + B*x + C*x^2 + D*x^3 + ...
 * where A, B, C, etc are either real or complex coefficients.
 * </p>
 * <p>
 * Includes a method for finding the zeros of a complex polynomial by the three
 * stage complex algorithm of Jenkins and Traub. The method finds the zeros one
 * at a time in roughly increasing order of modulus and deflates the polynomial
 * to one of lower degree. This method is extremely fast and timing is quite
 * insensitive to distribution of zeros.
 * </p>
 *
 * <p>
 * Zero finder ported from FORTRAN version of algorithm 419 courtesy <a
 * href="http://www.netlib.org/">Netlib Repository</a>.
 * </p>
 *
 * <p>
 * This class is not thread-safe. </p>
 *
 * <p>
 * Modified by: Joseph A. Huwaldt 
 *              Wilco Oelen         </p>
 * 
 * @author Joseph A. Huwaldt    Date: July 15, 2000
 * @author Wilco Oelen
 * @version October 6, 2015
 *
 */
public strictfp class CPoly {
    // Array of complex coefficients in order of increasing power.
    private final double[] coef_re;
    private final double[] coef_im;

    //-----------------------------------------------------------------------------------
    
    
    // Default constructor made unavailable for calling environments.
    private CPoly() {
        coef_re = null;
        coef_im = null;
    }

    
    
    /**
     * Constructor that takes an array of real coefficients. Coefficients are
     * supplied in order of increasing power. Example: p(x) = A + B*x + C*x^2 +
     * D*x^3 + E*x^4 gives coefficients[] = {A, B, C, D, E}. The degree of the
     * polynomial is equal to the length of the supplied array minus 1.
     *
     * @param coef An array of real coefficients in order of increasing power.
     */
    public CPoly(double[] coef) {
        coef_re = coef.clone();
        coef_im = new double[coef_re.length];
        for (int i=0; i<coef_im.length; i++) {
            coef_im[i] = 0.0;
        }
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
    public CPoly(double[] coef_re, double[] coef_im) {
        if ((coef_re == null || coef_re.length==0) && 
            (coef_im == null || coef_im.length==0)) {
            throw new RuntimeException("Construction of CPoly with empty coefficient set.");
        }
        if (coef_re == null) {
            coef_re = new double[0];
        }
        if (coef_im == null) {
            coef_im = new double[0];
        }
        
        int degRe = coef_re.length - 1;
        int degIm = coef_im.length - 1;
        
        if (degRe == degIm) {
            this.coef_re = coef_re.clone();
            this.coef_im = coef_im.clone();
            return;
        }
        
        if (degRe < degIm) {
            this.coef_re = new double[coef_im.length];
            int i;
            for (i=0; i<=degRe; i++) {
                this.coef_re[i] = coef_re[i];
            }
            for (; i<coef_im.length; i++) {
                this.coef_re[i] = 0.0;
            }
            this.coef_im = coef_im.clone();
        }
        else {
            this.coef_re = coef_re.clone();
            this.coef_im = new double[coef_re.length];
            int i;
            for (i=0; i<=degIm; i++) {
                this.coef_im[i] = coef_im[i];
            }
            for (; i<coef_re.length; i++) {
                this.coef_im[i] = 0.0;
            }
        }
    }
    
    
    
    
    /**
     * Evaluates the polynomial at a complex point, specified by the
     * arguments.
     * @param xre The real part of the point at which the polynomial
     * is evaluated.
     * @param xim The imaginary part of the point at which the polynomial
     * is evaluated.
     * @return A double[] array, with the value at index 0 being the
     * real part of the evaluated value and the value at index 1 being
     * the imaginary part of the evaluated value.
     */
    public double[] eval(double xre, double xim) {
        int deg = coef_re.length - 1;
        double re = coef_re[deg];
        double im = coef_im[deg];
        for (int i=deg-1; i>=0; i--) {
            double re2 = re*xre - im*xim + coef_re[i];
            double im2 = re*xim + im*xre + coef_im[i];
            re = re2;
            im = im2;
        }
        return new double[] {re, im};
    }

    
    
    
    
    /**
     * Evaluates the polynomial at the given complex argument. This is
     * like the other eval() method, but in this one the result of the
     * evaluation is put in an array, which already is allocated in the
     * calling environment. This method allows many evaluations of the
     * polynomial without the need to allocate many small 2-element arrays
     * for storing the result.
     *
     * @param result A double-array with at least two elements, in which
     * the result of the evaluation is stored. The value at index 0 is
     * replaced by the real part of the polynomial value and the value at
     * index 1 is replaced by the imaginary part of the polynomial value.
     * If the supplied array contains more than 2 elements, then the elements
     * with index 2 or larger are not touched at all.
     * @param xre The real part of the argument at which the polynomial is
     * evaluated.
     * @param xim The imaginary part of the argument at which the polynomial is
     * evaluated.
     */
    public void eval(double[] result, double xre, double xim) {
        int deg = coef_re.length - 1;
        double re = coef_re[deg];
        double im = coef_im[deg];
        for (int i=deg-1; i>=0; i--) {
            double re2 = re*xre - im*xim + coef_re[i];
            double im2 = re*xim + im*xre + coef_im[i];
            re = re2;
            im = im2;
        }
        result[0] = re;
        result[1] = im;
    }

    
    
    
    
    /**
     * Evaluates the polynomial and the derivative of the polynomial
     * simultaneously for the given value of x.
     *
     * @param result A double-array with at least two elements, in which
     * the result of the evaluation is stored. The value at index 0 is
     * replaced by the real part of the polynomial value and the value at
     * index 1 is replaced by the imaginary part of the polynomial value.
     * If the supplied array contains more than 2 elements, then the elements
     * with index 2 or larger are not touched at all.
     * @param dresult A double-array with at least two elements, in which
     * the result of the evaluation is stored. The value at index 0 is
     * replaced by the real part of the derivative value and the value at
     * index 1 is replaced by the imaginary part of the derivative value.
     * If the supplied array contains more than 2 elements, then the elements
     * with index 2 or larger are not touched at all.
     * @param xre The real part of the argument at which the derivative is
     * evaluated.
     * @param xim The imaginary part of the argument at which the derivative is
     * evaluated.
     */
    public void eval_deriv(double[] result, double[] dresult, double xre, double xim) {
        int deg = coef_re.length - 1;
        
        // p = coef[deg]
        double pre = coef_re[deg];
        double pim = coef_im[deg];
        
        // dp = 0
        double dpre = 0.0;
        double dpim = 0.0;
        
        for (int i = deg - 1; i >= 0; i--) {
            // dp = dp*x + p
            double dpre2 = dpre * xre - dpim * xim + pre;
            double dpim2 = dpre * xim + dpim * xre + pim;
            dpre = dpre2;
            dpim = dpim2;
            
            // p = p*x + coef[i];
            double pre2 = pre * xre - pim * xim + coef_re[i];
            double pim2 = pre * xim + pim * xre + coef_im[i];
            pre = pre2;
            pim = pim2;
        }
        result[0] = pre;
        result[1] = pim;
        dresult[0] = dpre;
        dresult[1] = dpim;
    }
    
    
    
    
    /**
     * Get the degree of the polynomial.
     * @return The degree of the polynomial.
     */
    public int degree() {
        return coef_re.length - 1;
    }
    
    

    
    
    
    /**
     * Creates a String representation of this polynomial.
     *
     * @return The String representation of this object.
	*
     */
    @Override
    public String toString() {
        StringBuilder buffer = new StringBuilder();
        int NN = coef_re.length;
        for (int i=0; i<NN; i++) {
            if (Math.abs(coef_re[i]) + Math.abs(coef_im[i]) > 0.0) {
                String term = "       x^" + i + " * (" + toString(coef_re[i], coef_im[i]) + ")\n";
                buffer.append(term);
            }
        }
        return buffer.toString();
    }
    
    
    static private String toString(double re, double im) {
        if (im == 0.0) {
            return "" + re;
        }
        if (re == 0.0) {
            if (im < 0) {
                return "-i*" + (-im);
            }
            return "i*" + im;
        }
        
        double r = re;
        if (r < 0) {
            r = -r;
        }
        
        double i = im;
        boolean imIsNeg = false;
        if (i < 0) {
            i = -i;
            imIsNeg = true;
        }
        
        double r_i = r/i;
        if (r_i > 4e15) {
            return "" + re;
        }
        if (r_i < 2.5e-16) {
            return (imIsNeg ? "i*" : "-i*") + i;
        }
        
        return "" + re + (imIsNeg ? " + i*" : " - i*") + i;
    }

    
    
    
    
    /**
     * This method computes the roots of the polynomial and stores the roots
     * in preallocated arrays, which are passed as arguments.
     * @param zeror Array, in which the real parts of the zeros will be stored
     * after computation of the zeros. This array must have a length of at
     * least N elements, where N is the degree of the polynomial.
     * @param zeroi Array, in which the imaginary parts of the zeros will be stored
     * after computation of the zeros. This array must have a length of at
     * least N elements, where N is the degree of the polynomial.
     * @return Returns the degree of the polynomial if the computation succeeds,
     * and returns a value less than the degree of the polynomial if an error
     * occurs (e.g. convergence failure). When a value less than the degree of
     * the polynomial is returned, then only part (or none) of the roots could
     * be determined.
     */
    public int solve(double[] zeror, double[] zeroi) {
        int nzeros = cpoly(zeror, zeroi);
        return nzeros;
    }

    
    
    
    /***********************************************************************/
    /***********************************************************************/
    /********* Below follows the port of the original Fortran  *************/
    /********* program. It is all private to this module!     *************/
    /********* The code above is a wrapper for easy usage. *************/
    /***********************************************************************/
    /***********************************************************************/
    /***********************************************************************/
    
    


    // The base of the number system being used.
    private static final double BASE = 2;

    // The number of base digits in each floating-point number (double precision)
    private static final double kT = 53;

    // The maximum relative representation error.  Fortran code:  BASE**(1-T)
    private static final double ETA = Math.pow(BASE, 1 - kT);

    // Infinity.  FORTRAN code:  BASE*(1.0D0-BASE**(-T))*BASE**(M-1)
    private static final double INFIN = 1e40;   // Double.MAX_VALUE;  

    // The smallest number that can be represented.  Fortran code:  (BASE**(N+3))/BASE**3
    private static final double SMALNO = 1e-40; // Double.MIN_NORMAL; 

    // Error bounds on complex addition.
    private static final double ARE = ETA;

    // Error bounds on complex multiplication.
    private static final double MRE = 2.0 * Math.sqrt(2.0) * ETA;
    
    private static final double DEG_TO_RAD = 3.14159265358979323846 / 180;  // Degrees-to-radians conversion factor = PI/180
    private static final double COSR = Math.cos(94.0 * DEG_TO_RAD);         // = -0.069756474;
    private static final double SINR = Math.sin(94.0 * DEG_TO_RAD);         // = 0.99756405

    
    // **** The following are used only by the root finding routines below. ****
    // Temporary storage space for complex numbers with real and imaginary parts.
    private double PVr, PVi, Tr, Ti, Sr, Si, Zr, Zi;
    
    
    //-----------------------------------------------------------------------------------
    /*
     * <p>
     * Finds all the zeros of a complex polynomial. </p>
     *
     * <p>
     * This program finds all the zeros of a complex polynomial by the three
     * stage complex algorithm of Jenkins and Traub. The program finds the zeros
     * one at a time in roughly increasing order of modulus and deflates the
     * polynomial to one of lower degree. The program is extremely fast and
     * timing is quite insenstive to distribution of zeros.
     * </p>
     *
     * <p>
     * Ported from FORTRAN to Java by Joseph A. Huwaldt, July 20, 2000 and
     * further enhanced by Wilco Oelen, 2014, 2015. 
     * Changes by Wilco Oelen:
     *   - When the iteration can be stopped, then one more iteration step
     *     is performed. This increases the accuracy of the found roots in
     *     nearly all cases. This solves a known old problem of inaccurately
     *     determined roots, even for well-conditioned polynomials.
     *   - Introduced tighter error bounds for termination of iteration. 
     *     This strongly improves the accuracy of the roots in many cases.
     *   - Added some wrapper code to make the use of the software easier
     *     and really 'black box'.
     * </p>
     *
     * <p>
     * FORTRAN version from <a href="http://www.netlib.org/">Netlib
     * Repository</a>
     * where it is listed as "419.f". A PDF file describing the algorithm and
     * its history with references is also available from Netlib.
     * </p>
     *
     * <p>
     * FORTRAN version had the following note at the top. ALGORITHM 419
     * COLLECTED ALGORITHMS FROM ACM. ALGORITHM APPEARED IN COMM. ACM, VOL. 15,
     * NO. 02, P. 097. Original Algol 60 zpolyzerofinder version by Jenkins,
     * 1969. (as noted in PDF scan of original ACM document).
     * </p>
     *
     * @param zeros_re Array, which will be filled with the real part of the zeros.
     * @param zeros_im Array, which will be filled with the imaginary part of the zeros.
     * @return The degree of the polynomial if all is OK.
	*
     */
    private int cpoly(double[] zeros_re, double[] zeros_im) {
        
        int NN = coef_re.length;
        int degree = NN - 1;
        double[] Pr = new double[NN];
        double[] Pi = new double[NN];
        for (int i = 0; i < NN; ++i) {
            Pr[i] = coef_re[degree - i];
            Pi[i] = coef_im[degree - i];
        }

        //  Algorithm fails if the leading coefficient is zero.
        if (Pr[0] == 0.0 && Pi[0] == 0.0) {
            return -1;
        }

        //  Allocate memory for arrays used by this method.
        double[] Hr = new double[NN];
        double[] Hi = new double[NN];
        double[] QPr = new double[NN];
        double[] QPi = new double[NN];
        double[] QHr = new double[NN];
        double[] QHi = new double[NN];
        double[] SHr = new double[NN];
        double[] SHi = new double[NN];

        //  Initialization of variables.
        double XX = Math.sqrt(0.5);
        double YY = -XX;

        //  Remove zeros at the origin, if any.
        while (Pr[NN - 1] == 0.0 && Pi[NN - 1] == 0.0) {
            int idNN2 = degree - NN + 1;
            zeros_re[idNN2] = 0.0;
            zeros_im[idNN2] = 0.0;
            --NN;
        }

        //  Calculate the modulus of the coefficients.
        for (int i = 0; i < NN; ++i) {
            SHr[i] = cmod(Pr[i], Pi[i]);
        }

        //  Scale the polynomial if needed.
        int factorExponent = scale(NN, SHr);
        if (factorExponent != 0) {
            for (int i = 0; i < NN; ++i) {
                Pr[i] = Math.scalb(Pr[i], factorExponent);
                Pi[i] = Math.scalb(Pi[i], factorExponent);
            }
        }

        // Start the algorithm for one zero.
outer:
        while (true) {
            if (NN <= 2) {
                // Calculate the final zero and return.
                cdiv(-Pr[1], -Pi[1], Pr[0], Pi[0]);       // Outputs Tr, Ti.
                zeros_re[degree - 1] = Tr;                // Outputs Tr, Ti.
                zeros_im[degree - 1] = Ti;
                return degree;
            }

            //  Calculate a lower bound on the modulus of the zeros.
            for (int i = 0; i < NN; ++i) {
                SHr[i] = cmod(Pr[i], Pi[i]);
            }
            double bound = cauchy(NN, SHr, SHi);

            // Outer loop to control 2 major passes with different sequences of shifts.
            for (int cnt1 = 1; cnt1 <= 2; ++cnt1) {

                // First stage calculation, no shift.
                noShift(NN, 5, Pr, Pi, Hr, Hi);

                // Inner loop to select a shift.
                for (int cnt2 = 1; cnt2 <= 9; ++cnt2) {
                    //  Shift is chosen with a modulus bound and amplitude rotated
                    //  by 94 degrees from the previous shift.
                    double XXX = COSR * XX - SINR * YY;
                    YY = SINR * XX - COSR * YY;
                    XX = XXX;
                    Sr = bound * XX;
                    Si = bound * YY;

                    //  Second stage calculation, fixed shift.
                    boolean conv = fxShift(NN, 10 * cnt2, Pr, Pi, QPr, QPi, Hr, Hi,
                            QHr, QHi, SHr, SHi);		// Outputs Zr, Zi.
                    if (conv) {
                        //  If successful the zero is stored and the polynomial deflated.
                        int idNN2 = degree - NN + 1;
                        zeros_re[idNN2] = Zr;
                        zeros_im[idNN2] = Zi;
                        --NN;
                        for (int i = 0; i < NN; ++i) {
                            Pr[i] = QPr[i];
                            Pi[i] = QPi[i];
                        }

                        // The 2nd stage jumps directly back to 3rd stage iteration.
                        continue outer;
                    }

                    //  If iteration is unsuccessful, another shift is chosen.
                }

                //  If 9 shifts fail, the outer loop is repeated with another
                //  sequence of shifts.
            }

            //  The zero finder has failed on two major passes.  Return empty handed.
            return degree - (NN - 1);
        }

        // We can never get here.
    }
    
    
    

    /**
     * Evaluate this polynomial at a complex x and returns the generally complex
     * result as a pair of class variables. Class variables are used to avoid
     * the overhead of creating a "Complex" object during root finding. Sets the
     * class variables PVr, and PVi with the real and imaginary parts of the
     * result. Uses the method of Horner Recurrence.
     *
     * @param NN The number of coefficients to use in the evaluation.
     * @param Sr The real component of the point at which to evaluate this
     * polynomial.
     * @param Si The imaginary component of the point at which to evaluate this
     * polynomial.
     * @param Pr, Pi Real & Imaginary coefficients of the polynomial.
     * @param Qr, Qi Arrays to contain partial sums.
	*
     */
    private void polyEv(int NN, double Sr, double Si, double[] Pr, double[] Pi,
                        double[] Qr, double[] Qi) {
        // Begin evaluation.
        double pvr = Qr[0] = Pr[0];
        double pvi = Qi[0] = Pi[0];

        for (int i = 1; i < NN; ++i) {
            double temp = pvr;
            pvr = pvr * Sr - pvi * Si + Pr[i];
            pvi = temp * Si + pvi * Sr + Pi[i];

            Qr[i] = pvr;
            Qi[i] = pvi;
        }

        // Use a class variable to pass results back when doing root finding.
        PVr = pvr;
        PVi = pvi;
    }

    
    
    
    /**
     * Bounds the error in evaluating the polynomial by the method of Horner
     * Recurrance.
     *
     * @param NN The number of coefficients to use in the evaluation.
     * @param Qr Real part of partial sum from evaluate().
     * @param Qi Imaginary part of partial sum from evaluate().
     * @param MS Modulus of the point being evaluated.
     * @param MP Modulus of the polynomial value.
     * @param ARE Error bounds on complex addition.
     * @param MRE Error bounds on complex multiplication.
     *
     */
    private static double errEv(int NN, double[] Qr, double[] Qi, double MS, double MP,
                                double ARE, double MRE) {
        double E = cmod(Qr[0], Qi[0]) * MRE / (ARE + MRE);
        for (int i = 0; i < NN; ++i) {
            E = E * MS + cmod(Qr[i], Qi[i]);
        }

        return E*(ARE + MRE) - MP*MRE;
    }

    
    
    
    /**
     * Computes a lower bound on the moduli of the zeros of a polynomial.
     *
     * @param NN The number of coefficients to use in the evaluation.
     * @param PT The modulus of the coefficients of the polynomial.
     * @param Q Array filled in on output.
     *
     */
    private static double cauchy(int NN, double[] PT, double[] Q) {
        int NNm1 = NN - 1;

        PT[NNm1] = -PT[NNm1];

        //	Compute the upper estimate of bound.
        int N = NN - 1;
        int Nm1 = N - 1;
        double X = Math.exp((Math.log(-PT[NNm1]) - Math.log(PT[0])) / N);
        if (PT[Nm1] != 0.) {
            // If newton step at the origin is better, use it.
            double XM = -PT[NNm1] / PT[Nm1];
            if (XM < X) {
                X = XM;
            }
        }
        // Chop the interval (0,X) until F <= 0.
        while (true) {
            double XM = X * 0.1;
            double F = PT[0];
            for (int i = 1; i < NN; ++i) {
                F = F * XM + PT[i];
            }
            if (F <= 0.) {
                break;
            }
            X = XM;
        }
        double DX = X;

        // Do newton iteration until X converges to two decimal places.
        while (Math.abs(DX / X) > 0.005) {
            Q[0] = PT[0];
            for (int i = 1; i < NN; ++i) {
                Q[i] = Q[i - 1] * X + PT[i];
            }
            double F = Q[NNm1];
            double DF = Q[0];
            for (int i = 1; i < N; ++i) {
                DF = DF * X + Q[i];
            }
            DX = F / DF;
            X = X - DX;
        }

        return X;
    }

    
    
    
    /**
     * Returns a scale factor to multiply the coefficients of the polynomial.
     * The scaling is done to avoid overflow and to avoid undetected underflow
     * interfering with the convergence criterion. The factor is a power of the
     * BASE.
     *
     * @param NN The number of coefficients to use in the evaluation.
     * @param PT The modulus of the coefficients of the polynomial.
     *
     */
    private static int scale(int NN, double[] PT) {
        // Find the largest and the smallest moduli of coefficients.
        double hi = Math.sqrt(INFIN);
        double lo = SMALNO / ETA;
        double max = 0.;
        double min = INFIN;
        double X, sc;
        for (int i = 0; i < NN; ++i) {
            X = PT[i];
            if (X > max) {
                max = X;
            }
            if (X != 0.0 && X < min) {
                min = X;
            }
        }

        // Scale only if there are very large or very small components.
        if (min >= lo && max <= hi) {
            return 0;
        }

        X = lo / min;
        if (X > 1.0) {
            sc = X;
            if (INFIN / sc > max) {
                return 0;
            }
        } else {
            sc = 1.0 / Math.sqrt(max * min);
        }

        // Compute a scale factor, close to sc, but in such a way
        // that the factor is exact and a power of 2, such that
        // multiplication with this factor does not lead to loss
        // of any precision in the coefficients. Not the scale 
        // factor itself, but the exponent of 2 is returned.
        return Math.getExponent(sc) + 1;
    }

    
    
    
    /**
     * Complex division C = A/B, avoiding overflow. Results are stored in class
     * variables Tr and Ti to avoid overhead of creating a Complex object during
     * root finding. Results are stored in class variables Tr and Ti.
     *
     * @param Ar The real part of the complex numerator.
     * @param Ai The imaginary part of the complex numerator.
     * @param Br The real part of the complex denominator.
     * @param Bi The imaginary part of the complex denominator.
     *
     */
    private void cdiv(double Ar, double Ai, double Br, double Bi) {
        if (Br == 0. && Bi == 0.) {
            // Division by zero, result = infinity.
            Tr = INFIN;
            Ti = INFIN;
            return;
        }

        if (Math.abs(Br) >= Math.abs(Bi)) {
            double R = Bi / Br;
            double D = Br + R * Bi;
            Tr = (Ar + Ai * R) / D;
            Ti = (Ai - Ar * R) / D;

        } 
        else {
            double R = Br / Bi;
            double D = Bi + R * Br;
            Tr = (Ar * R + Ai) / D;
            Ti = (Ai * R - Ar) / D;
        }

    }

    
    
    
    /**
     * Calculates the modulus or magnitude of a complex number avoiding
     * overflow.
     *
     * Adapted from "Numerical Recipes in C: The Art of Scientific Computing"
     * 2nd Edition, pg 949, ISBN 0-521-43108-5. The NR algorithm is only
     * slightly different from the ACM algorithm, but the NR version appears to
     * be slightly more robust.
     *
     * @param re The real part of the complex number.
     * @param im The imaginary part of the complex number.
     *
     */
    private static double cmod(double re, double im) {
        double ans;
        re = Math.abs(re);
        im = Math.abs(im);

        if (re == 0.0) {
            ans = im;
        } 
        else if (im == 0.0) {
            ans = re;
        } 
        else if (re > im) {
            double temp = im / re;
            ans = re * Math.sqrt(1.0 + temp * temp);
        } 
        else {
            double temp = re / im;
            ans = im * Math.sqrt(1.0 + temp * temp);
        }

        return ans;
    }

    
    
    
    /**
     * Computes the derivative polynomial as the intial H polynomial and
     * computes L1 no-shift H polynomials.
     *
     * @param NN The number of coefficients to use in the evaluation.
     * @param L1 Number of Level 1 shifts to make.
     * @param Pr, Pi The coefficients of the polynomial.
     * @param Hr, Hi Arrays containing output ?
     *
     */
    private void noShift(int NN, int L1, double[] Pr, double[] Pi, double[] Hr, double[] Hi) {
        int N = NN - 1;
        int Nm1 = N - 1;
        int NNm1 = NN - 1;

        for (int i = 0; i < N; ++i) {
            double XNi = NNm1 - i;
            Hr[i] = XNi * Pr[i] / N;
            Hi[i] = XNi * Pi[i] / N;
        }

        for (int jj = 1; jj <= L1; ++jj) {
            if (cmod(Hr[Nm1], Hi[Nm1]) > ETA * 10.0 * cmod(Pr[Nm1], Pi[Nm1])) {
                // Divide the negative coefficient by the derivative.
                cdiv(-Pr[NNm1], -Pi[NNm1], Hr[Nm1], Hi[Nm1]);	// Outputs Tr, Ti.
                for (int i = 1; i <= Nm1; ++i) {
                    int j = NNm1 - i;
                    double T1 = Hr[j - 1];
                    double T2 = Hi[j - 1];
                    Hr[j] = Tr * T1 - Ti * T2 + Pr[j];
                    Hi[j] = Tr * T2 + Ti * T1 + Pi[j];
                }
                Hr[0] = Pr[0];
                Hi[0] = Pi[0];

            } 
            else {
                // If the constant term is essentially zero, shift H coefficients.
                for (int i = 1; i <= Nm1; ++i) {
                    int j = NNm1 - i;
                    Hr[j] = Hr[j - 1];
                    Hi[j] = Hi[j - 1];
                }
                Hr[0] = 0.;
                Hi[0] = 0.;
            }
        }

    }

    
    
    
    /**
     * Computes T = -P(S)/H(S). Sets class variables Tr, Ti.
     *
     * @param NN The number of coefficients to use in the evaluation.
     * @param Sr The real part of the point we are evaluating the polynomial at.
     * @param Si The imaginary part of the point we are evaluating the
     * polynomial at.
     * @param Hr, Hi Arrays containing ?
     * @param QHr, QHi Arrays containing partial sums of H(S) polynomial.
     * @return True if H(S) is essentially zero.
     *
     */
    private boolean calcT(int NN, double Sr, double Si, double[] Hr, double[] Hi,
            double[] QHr, double[] QHi) {
        int N = NN - 1;
        int Nm1 = N - 1;

        //	Evaluate H(S).
        double tempR = PVr;
        double tempI = PVi;
        polyEv(N, Sr, Si, Hr, Hi, QHr, QHi);
        double HVr = PVr;
        double HVi = PVi;
        PVr = tempR;
        PVi = tempI;

        // Is H(S) essentially zero?
        boolean nearZero = cmod(HVr, HVi) <= ARE * 10.0 * cmod(Hr[Nm1], Hi[Nm1]);
        if (nearZero) {
            Tr = 0.;
            Ti = 0.;
        } 
        else {
            cdiv(-PVr, -PVi, HVr, HVi);		// Outputs Tr, Ti.
        }
        return nearZero;
    }

    
    
    
    
    /**
     * Calculates the next shifted H polynomial.
     *
     * @param NN The number of coefficients to use in the evaluation.
     * @param bool Set to true if H(S) is essentially zero.
     * @param Hr, Hi Arrays containing ?
     * @param QPr, QPi
     * @param QHr, QHi
     *
     */
    private void nextH(int NN, boolean bool, double[] Hr, double[] Hi,
            double[] QPr, double[] QPi, double[] QHr, double[] QHi) {
        int N = NN - 1;

        if (!bool) {
            for (int j = 1; j < N; ++j) {
                double T1 = QHr[j - 1];
                double T2 = QHi[j - 1];
                Hr[j] = Tr * T1 - Ti * T2 + QPr[j];
                Hi[j] = Tr * T2 + Ti * T1 + QPi[j];
            }
            Hr[0] = QPr[0];
            Hi[0] = QPi[0];

        } 
        else {
            // If H(S) is zero, replace H with QH.
            for (int j = 1; j < N; ++j) {
                Hr[j] = QHr[j - 1];
                Hi[j] = QHi[j - 1];
            }
            Hr[0] = 0.;
            Hi[0] = 0.;
        }

    }

    
    
    
    
    /**
     * Carries out the third stage iteration. On entry Zr, Zi contains the
     * initial iteration. If the iteration converges it contains the final
     * iteration on exit. Also uses and sets class variables Sr, Si.
     *
     * @param NN The number of coefficients to use in the evaluation.
     * @param L3 Limit of steps in stage 3.
     * @param Pr, Pi The coefficients of the polynomial.
     * @param QPr, QPi
     * @param Hr, Hi Arrays containing ?
     * @param QHr, QHi
     * @return True if iteration converges.
     *
     */
    private boolean vrShift(int NN, int L3, double[] Pr, double[] Pi,
                            double QPr[], double QPi[],
                            double Hr[], double Hi[], double QHr[], double QHi[]) {
        boolean conv = false;
        boolean B = false;
        double OMP = 0., RelSTP = 0.;

        Sr = Zr;
        Si = Zi;

        // Main loop for stage three.
        // WILCO: Added nOK-logic to force one additional iteration
        // after it is decided that iteration can be terminated.
        int nOK = 0;
        for (int i = 1; i <= L3; ++i) {
            //	Evaluate P at S and test for convergence.
            polyEv(NN, Sr, Si, Pr, Pi, QPr, QPi);     // Outputs PVr, PVi.

            double MP = cmod(PVr, PVi);
            double MS = cmod(Sr, Si);
            if (MP <= 20.0 * errEv(NN, QPr, QPi, MS, MP, ARE, MRE)) {
                nOK++;
                if (nOK == 2) {
                    // Polynomial value is smaller in value than a bound on the error
                    // in evaluating P, terminate the iteration.
                    Zr = Sr;
                    Zi = Si;
                    return true;
                }
            }
            else {
                nOK = 0;
            }

            if (i == 1) {
                OMP = MP;
            }
            else {
                if (!B && MP >= OMP && RelSTP < 0.05) {
                    //	Iteration has stalled. Probably a cluster of zeros. Do 5 fixed
                    //	shift steps into the cluster to force one zero to dominate.
                    double TP = RelSTP;
                    B = true;
                    if (RelSTP < ETA) {
                        TP = ETA;
                    }
                    double R1 = Math.sqrt(TP);
                    double R2 = Sr * (1. + R1) - Si * R1;
                    Si = Sr * R1 + Si * (1. + R1);
                    Sr = R2;
                    polyEv(NN, Sr, Si, Pr, Pi, QPr, QPi);             // Outputs PVr, PVi.

                    for (int j = 0; j < 5; ++j) {
                        boolean bool = calcT(NN, Sr, Si, Hr, Hi, QHr, QHi);   // Outputs Tr, Ti.
                        nextH(NN, bool, Hr, Hi, QPr, QPi, QHr, QHi);
                    }

                    OMP = INFIN;
                } 
                else {
                    // Exit if polynomial value increases significantly.
                    if (MP * 0.1 > OMP) {
                        return conv;
                    }
                    OMP = MP;
                }
            }

            // Calculate next iteration.
            boolean bool = calcT(NN, Sr, Si, Hr, Hi, QHr, QHi);       // Outputs Tr, Ti.
            nextH(NN, bool, Hr, Hi, QPr, QPi, QHr, QHi);
            bool = calcT(NN, Sr, Si, Hr, Hi, QHr, QHi);       // Outputs Tr, Ti.
            if (!bool) {
                RelSTP = cmod(Tr, Ti) / cmod(Sr, Si);
                Sr = Sr + Tr;
                Si = Si + Ti;
            }
        }

        return conv;
    }

    
    
    
    /**
     * Computes L2 fixed-shift H polynomials and tests for convergence.
     * Initiates a variable-shift iteration and returns with the approximate
     * zero if successfull. Uses and sets the class variables Sr and Si. Sets
     * class variables Zr, Zi to approximate zero if convergence is true.
     *
     * @param NN The number of coefficients to use in the evaluation.
     * @param L2 Limit of fixed shift steps.
     * @param Pr, Pi The coefficients of the polynomial.
     * @param QPr, QPi
     * @param Hr, Hi Arrays containing ?
     * @param QHr, QHi
     * @return True if convergence of stage 3 iteration is successfull.
     *
     */
    private boolean fxShift(int NN, int L2, double[] Pr, double[] Pi,
            double QPr[], double QPi[],
            double Hr[], double Hi[], double QHr[], double QHi[],
            double SHr[], double SHi[]) {
        int N = NN - 1;

        // Evaluate Polynomial at S.
        polyEv(NN, Sr, Si, Pr, Pi, QPr, QPi);	// Outputs PVr, PVi.
        boolean test = true;
        boolean pasd = false;

        // Calculate 1st T = -P(S)/H(S).
        boolean bool = calcT(NN, Sr, Si, Hr, Hi, QHr, QHi);   // Outputs Tr, Ti.

        // Main loop for one 2nd stage step.
        for (int j = 1; j <= L2; ++j) {
            double OTr = Tr;
            double OTi = Ti;

            // Compute next H polynomial and new T.
            nextH(NN, bool, Hr, Hi, QPr, QPi, QHr, QHi);
            bool = calcT(NN, Sr, Si, Hr, Hi, QHr, QHi);    // Outputs Tr, Ti.
            Zr = Sr + Tr;
            Zi = Si + Ti;

            // Test for convergence unless stage 3 has failed once or
            // this is the last H polynomial.
            if (!bool && test && j != L2) {
                if (cmod(Tr - OTr, Ti - OTi) < 0.5 * cmod(Zr, Zi)) {
                    if (pasd) {
                        // The weak convergence test has been passed twice, start
                        // the third stage iteration, after saving the current H
                        // polynomial and shift.
                        for (int i = 0; i < N; ++i) {
                            SHr[i] = Hr[i];
                            SHi[i] = Hi[i];
                        }
                        double SVSr = Sr;
                        double SVSi = Si;
                        boolean conv = vrShift(NN, 10, Pr, Pi, QPr, QPi, Hr, Hi, QHr, QHi);    // Outputs Zr, Zi.
                        if (conv) {
                            return conv;
                        }

                        // The iteration failed to converge.  Turn off testing and
                        // restore H, S, PV and T.
                        test = false;
                        for (int i = 0; i < N; ++i) {
                            Hr[i] = SHr[i];
                            Hi[i] = SHi[i];
                        }
                        Sr = SVSr;
                        Si = SVSi;

                        polyEv(NN, Sr, Si, Pr, Pi, QPr, QPi);        // Outputs PVr, PVi.
                        bool = calcT(NN, Sr, Si, Hr, Hi, QHr, QHi);  // Outputs Tr, Ti.
                        continue;
                    }
                    pasd = true;
                } 
                else {
                    pasd = false;
                }
            }
        }

        // Attempt an iteration with final H polynomial from second stage.
        boolean conv = vrShift(NN, 10, Pr, Pi, QPr, QPi, Hr, Hi, QHr, QHi);    // Outputs Zr, Zi.
        return conv;
    }
}
