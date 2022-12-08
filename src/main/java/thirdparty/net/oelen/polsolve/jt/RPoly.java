/*
 *   RPoly  -- A class that represents a polynomial equation.
 *
 *   Copyright (C) 2014, 2015 by Wilco Oelen
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
 * Represents a polynomial equation of the form p(x) = A + B*x + C*x^2 + D*x^3 +
 * ... where A, B, C, etc are real coefficients.
 * </p>
 * <p>
 * Includes a method for finding the zeros of a real polynomial by a specialized
 * version of the algorithm of Jenkins and Traub, developed especially for poly-
 * nomials with real coefficients. This specialized method is approximately 3
 * times as fast as the general algorithm when applied to the same polynomial.
 * Like the general algorithm, the specialized method finds the zeros one at
 * a time in roughly increasing order of modulus and deflates the polynomial
 * to one of lower degree.
 * </p>
 *
 * <p>
 * Zero finder ported from FORTRAN version of algorithm 493 courtesy <a
 * href="http://www.netlib.org/">Netlib Repository</a>.
 * </p>
 *
 * <p>
 * This class is not thread-safe. </p>
 *
 * <p>
 * Modified by: Wilco Oelen         </p>
 * 
 * @author Wilco Oelen
 * @version October 6, 2015
 *
 */

public strictfp class RPoly {   // Use strictfp to assure behavior is same on all platforms

    private final double coef[];

    
    
    
    // Make this constructor unavailable for calling environments.
    private RPoly() {
        coef = null;
    }

    
        
    /**
     * Constructor of the polynomial, with the coefficients given as an array of
     * doubles.
     *
     @param coef The coefficients of the polynomial in the order of the
     * constant term to the term with the highest power, i.e. the polynomial is
     * coef[0] + coef[1]*x + coef[2]*x^2 ... coef[N]*x^N, where N+1 is the
     * length of the supplied array. Element coef[N] must be non-zero!
     */
    public RPoly(double[] coef) {
        this.coef = coef.clone();
    }

    
    
    
    
    /**
     * Evaluates the polynomial at the given double argument.
     *
     @param x The real value for which the polynomial must be evaluated. The
     * result is a real value.
     @return The value of the polynomial, evaluated at the given argument.
     */
    public double eval(double x) {
        int deg = coef.length - 1;
        double val = coef[deg];
        for (int i = deg - 1; i >= 0; i--) {
            val = val * x + coef[i];
        }
        return val;
    }

    
    
    
    
    /**
     * Evaluates the polynomial at the given complex argument.
     *
     @param xre The real part of the argument at which the polynomial is
     * evaluated.
     @param xim The imaginary part of the argument at which the polynomial is
     * evaluated.
     @return A double array, which contains the complex value of the
     * polynomial, evaluated at the given argument. The returned array contains
     * two real values. The value at index 0 is the real part of the return
     * value and the value at index 1 is the imaginary part of the return value.
     */
    public double[] eval(double xre, double xim) {
        int deg = coef.length - 1;
        double re = coef[deg];
        double im = 0.0;
        for (int i = deg - 1; i >= 0; i--) {
            double re2 = re * xre - im * xim + coef[i];
            double im2 = re * xim + im * xre;
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
     @param result A double-array with at least two elements, in which
     * the result of the evaluation is stored. The value at index 0 is
     * replaced by the real part of the polynomial value and the value at
     * index 1 is replaced by the imaginary part of the polynomial value.
     * If the supplied array contains more than 2 elements, then the elements
     * with index 2 or larger are not touched at all.
     @param xre The real part of the argument at which the polynomial is
     * evaluated.
     @param xim The imaginary part of the argument at which the polynomial is
     * evaluated.
     */
    public void eval(double[] result, double xre, double xim) {
        int deg = coef.length - 1;
        double pre = coef[deg];
        double pim = 0.0;
        for (int i = deg - 1; i >= 0; i--) {
            double pre2 = pre * xre - pim * xim + coef[i];
            double pim2 = pre * xim + pim * xre;
            pre = pre2;
            pim = pim2;
        }
        result[0] = pre;
        result[1] = pim;
    }

    
    
    
    
    /**
     * Evaluates the polynomial and the derivative of the polynomial
     * simultaneously for the given value of x.
     *
     @param result A double-array with at least two elements, in which
     * the result of the evaluation is stored. The value at index 0 is
     * replaced by the real part of the polynomial value and the value at
     * index 1 is replaced by the imaginary part of the polynomial value.
     * If the supplied array contains more than 2 elements, then the elements
     * with index 2 or larger are not touched at all.
     @param dresult A double-array with at least two elements, in which
     * the result of the evaluation is stored. The value at index 0 is
     * replaced by the real part of the derivative value and the value at
     * index 1 is replaced by the imaginary part of the derivative value.
     * If the supplied array contains more than 2 elements, then the elements
     * with index 2 or larger are not touched at all.
     @param xre The real part of the argument at which the derivative is
     * evaluated.
     @param xim The imaginary part of the argument at which the derivative is
     * evaluated.
     */
    public void eval_deriv(double[] result, double[] dresult, double xre, double xim) {
        int deg = coef.length - 1;
        
        // p = coef[deg]
        double pre = coef[deg];
        double pim = 0.0;
        
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
            double pre2 = pre * xre - pim * xim + coef[i];
            double pim2 = pre * xim + pim * xre;
            pre = pre2;
            pim = pim2;
        }
        result[0] = pre;
        result[1] = pim;
        dresult[0] = dpre;
        dresult[1] = dpim;
    }

    
    
    
    /**
     * Returns the degree of the polynomial.
     *
     @return The degree of the polynomial.
     */
    public int degree() {
        return coef.length - 1;
    }
    
    
    
    
    /**
     * Creates a String representation of this polynomial.
     *
     @return The String representation of this object.
	*
     */
    @Override
    public String toString() {
        StringBuilder buffer = new StringBuilder();

        int NN = coef.length;
        for (int i=0; i<NN; i++) {
            if (coef[i] < 0) {
                String term = "     - x^" + i + " * " + (-coef[i]) + "\n";
                buffer.append(term);
            }
            else if (coef[i] > 0) {
                String term = "       x^" + i + " * " + coef[i] + "\n";
                buffer.append(term);
            }
            else {
                // No row is generated for a zero term.
            }
        }

        return buffer.toString();
    }

    
    
    
    /**
     * This method computes the roots of the polynomial and stores the roots
     * in preallocated arrays, which are passed as arguments.
     @param zeror Array, in which the real parts of the zeros will be stored
     * after computation of the zeros. This array must have a length of at
     * least N elements, where N is the degree of the polynomial.
     @param zeroi Array, in which the imaginary parts of the zeros will be stored
     * after computation of the zeros. This array must have a length of at
     * least N elements, where N is the degree of the polynomial.
     @return Returns the degree of the polynomial if the computation succeeds,
     * and returns a value less than the degree of the polynomial if an error
     * occurs (e.g. convergence failure). When a value less than the degree of
     * the polynomial is returned, then only part (or none) of the roots could
     * be determined.
     */
    public int solve(double[] zeror, double[] zeroi) {
        p = new double[coef.length];
        int degree = p.length - 1;
        for (int i = 0; i < p.length; i++) {
            p[i] = coef[degree - i];
        }

        int nzeros = rpoly(degree, zeror, zeroi);
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
    
    
    
    /*      Jenkins-Traub real polynomial root finder.
     *
     *      Translation of TOMS493 from FORTRAN to Java. The code is
     *      reorganized strongly because Java does not have a goto
     *      statement. Besides the restructuring, there also are
     *      quite a few algorithmic changes. The original Fortran
     *      code is quite buggy. The most important changes are:
     *          - adopted scaling mechanism for polynomial coefficients
     *            from CPoly and removed the original malfunctioning
     *            scaling code of algorithm 493.
     *          - added logic to the termination of the iteration
     *            process, such that after the newly found root can be
     *            isolated, one more iteration step is performed, such
     *            that a more accurate result is obtained.
     *
     *      The calling conventions are slightly modified to return
     *      the number of roots found as the function value.
     *
     *      INPUT:
     *      op - double precision vector of coefficients in order of
     *              decreasing powers.
     *      degree - integer degree of polynomial
     *
     *      OUTPUT:
     *      zeror,zeroi - output double precision vectors of the
     *                    real and imaginary parts of the zeros.
     *
     *      RETURN:
     *      returnval: -1 if leading coefficient is zero, otherwise
     *                  number of roots found. 
     */

    private static final double DEG_TO_RAD = 3.14159265358979323846 / 180;  // Degrees-to-radians conversion factor = PI/180
    private static final double COSR = Math.cos(94.0 * DEG_TO_RAD);  // = -0.069756474
    private static final double SINR = Math.sin(94.0 * DEG_TO_RAD);  // = 0.99756405
    private static final double SQRT_0_5 = Math.sqrt(0.5);
    private static final double INFIN = 1e40;    // Double.MAX_VALUE;
    private static final double SMALNO = 1e-40;  // Double.MIN_NORMAL;
    private static final double ETA = 1.0 / (1l << 52);
    private static final double HI = Math.sqrt(INFIN);
    private static final double LO = SMALNO / ETA;
    private static final double ARE = ETA;
    private static final double MRE = ETA;

    
    private double[] p, qp, k, qk, svk;
    private double u, v;
    private double a, b;
    private double c, d;
    private double sr, a1;
    private double a3, a7, e, f, g, h;
    private double szr, szi, lzr, lzi;
    private int n;

    private int rpoly(int degree, double[] zeror, double[] zeroi) {
        double t, aa, bb, cc;
        double[] temp;
        double[] pt;
        double xx, yy, x, bnd;
        double xm, ff, df, dx;
        int cnt, nz, i, j, jj, nm1;
        boolean zerok;

        n = degree;
        
        /*  Initialization of constants for shift rotation. */
        xx = SQRT_0_5;
        yy = -xx;
        /*  Algorithm fails if the leading coefficient is zero. */
        if (p[0] == 0.0) {
            return -1;
        }
        /*  Remove the zeros at the origin, if any. */
        while (p[n] == 0.0) {
            j = degree - n;
            zeror[j] = 0.0;
            zeroi[j] = 0.0;
            n--;
        }
        if (n < 1) {
            return degree;
        }
        /*
         *  Allocate memory here
         */
        temp = new double[degree + 1];
        pt = new double[degree + 1];
        qp = new double[degree + 1];
        k = new double[degree + 1];
        qk = new double[degree + 1];
        svk = new double[degree + 1];

    outer:
        while (true) {
            /*  Start the algorithm for one zero. */
            if (n == 1) {
                zeror[degree - 1] = -p[1] / p[0];
                zeroi[degree - 1] = 0.0;
                n -= 1;
                return degree - n;
            }
            /*  Calculate the final zero or pair of zeros. */
            if (n == 2) {
                quad(p[0], p[1], p[2]);   // Returns its results through szr, szi, lzr, lzi
                zeror[degree - 2] = szr;
                zeroi[degree - 2] = szi;
                zeror[degree - 1] = lzr;
                zeroi[degree - 1] = lzi;
                n -= 2;
                return degree - n;
            }
            
            // WILCO: Replaced scaling code of the original Fortran code by
            // the scaling mechanism as is used in the original Fortran code
            // of CPoly. The code in the original RPoly scales the coefficients
            // to just above underflow level and this makes the code quite
            // sensitive to loss of accuracy. Probably this is one of the bugs
            // of the TOMS code, algorithm 493: http://www.netlib.org/toms/493
            // This code seems to have unnoticed errors, which may be lingering
            // for years already.
            int factorExponent = scale(p);
            if (factorExponent != 0) {
                for (i = 0; i <= n; i++) {
                    p[i] = Math.scalb(p[i], factorExponent);     // Scale polynomial with factor 2^factorExponent.
                }
            }

            /*  Compute lower bound on moduli of roots. */
            for (i = 0; i <= n; i++) {
                pt[i] = (Math.abs(p[i]));
            }
            pt[n] = -pt[n];
            /*  Compute upper estimate of bound. */
            x = Math.exp((Math.log(-pt[n]) - Math.log(pt[0])) / (double) n);
            /*  If Newton step at the origin is better, use it. */
            if (pt[n - 1] != 0.0) {
                xm = -pt[n] / pt[n - 1];
                if (xm < x) {
                    x = xm;
                }
            }
            /*  Chop the interval (0,x) until ff <= 0 */
            while (true) {
                xm = x * 0.1;
                ff = pt[0];
                for (i = 1; i <= n; i++) {
                    ff = ff * xm + pt[i];
                }
                if (ff <= 0.0) {
                    break;
                }
                x = xm;
            }
            dx = x;
            /*  Do Newton interation until x converges to two 
             *  decimal places. 
             */
            while (Math.abs(dx / x) > 0.005) {
                ff = pt[0];
                df = ff;
                for (i = 1; i < n; i++) {
                    ff = ff * x + pt[i];
                    df = df * x + ff;
                }
                ff = ff * x + pt[n];
                dx = ff / df;
                x -= dx;
            }
            bnd = x;
            /*  Compute the derivative as the initial k polynomial
             *  and do 5 steps with no shift.
             */
            nm1 = n - 1;
            for (i = 1; i < n; i++) {
                k[i] = (double) (n - i) * p[i] / (double) n;
            }
            k[0] = p[0];
            aa = p[n];
            bb = p[n - 1];
            zerok = (k[n - 1] == 0);
            for (jj = 0; jj < 5; jj++) {
                cc = k[n - 1];
                if (!zerok) {
                    /*  Use a scaled form of recurrence if value of k at 0 is nonzero. */
                    t = -aa / cc;
                    for (i = 0; i < nm1; i++) {
                        j = n - i - 1;
                        k[j] = t * k[j - 1] + p[j];
                    }
                    k[0] = p[0];
                    zerok = (Math.abs(k[n - 1]) <= Math.abs(bb) * ETA * 10.0);
                } else {
                    /*  Use unscaled form of recurrence. */
                    for (i = 0; i < nm1; i++) {
                        j = n - i - 1;
                        k[j] = k[j - 1];
                    }
                    k[0] = 0.0;
                    zerok = (k[n - 1] == 0.0);
                }
            }
            /*  Save k for restarts with new shifts. */
            for (i = 0; i < n; i++) {
                temp[i] = k[i];
            }
            /*  Loop to select the quadratic corresponding to each new shift. */
            for (cnt = 1; cnt <= 20; cnt++) {
                /*  Quadratic corresponds to a double shift to a            
                 *  non-real point and its complex conjugate. The point
                 *  has modulus bnd and amplitude rotated by 94 degrees
                 *  from the previous shift.
                 */
                double xxx = COSR * xx - SINR * yy;
                yy = SINR * xx + COSR * yy;
                xx = xxx;
                sr = bnd * xx;
                u = -2.0 * sr;
                v = bnd;
                nz = fxshfr(20 * cnt);
                if (nz != 0) {
                    /*  The second stage jumps directly to one of the third
                     *  stage iterations and returns here if successful.
                     *  Deflate the polynomial, store the zero or zeros and
                     *  return to the main algorithm.
                     */
                    j = degree - n;
                    zeror[j] = szr;
                    zeroi[j] = szi;
                    n -= nz;
                    for (i = 0; i <= n; i++) {
                        p[i] = qp[i];
                    }
                    if (nz != 1) {
                        zeror[j + 1] = lzr;
                        zeroi[j + 1] = lzi;
                    }
                    continue outer;
                }
                /*  Iteration is unsuccessful: another quadratic
                 *  is chosen after restoring k.
                 */
                for (i = 0; i < n; i++) {
                    k[i] = temp[i];
                }
            }

            return degree - n;
        }
    }

    
    
    
    // WILCO: Deviation from original algorithm 493. In the original
    // version of the algorithm the polynomial is scaled such that
    // all coefficients become very small, in the order of magnitude
    // of 10^-40. This new scale() method tries to keep coefficients 
    // around the order of magnitude 1.0.
    private static int scale(double[] PT) {
        // Find the largest and the smallest moduli of coefficients.
        int NN = PT.length;
        double max = 0.0;
        double min = INFIN;
        double X, sc;
        for (int i = 0; i < NN; ++i) {
            X = Math.abs(PT[i]);
            if (X > max) {
                max = X;
            }
            if (X != 0.0 && X < min) {
                min = X;
            }
        }

        // Scale only if there are very large or very small components.
        if (min >= LO && max <= HI) {
            return 0;
        }

        X = LO / min;
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

    
    
    private
    double ui, vi;  // Used from within fxshfr() and initialized in newest().
                    // ui and vi are used as input to quadit() and also are
                    // changed by quadit(), also by calling newest(). In a
                    // language like C, ui and vi could be local in fxshfr()
                    // and pointers could be passed to newest() and quadit(),
                    // but in Java they are made class variables.
    
    private
    double s;       // A similar issue as for ui and vi exists for the local
                    // variable s in fxshfr(). In the fortran code, this is
                    // passed as a VAR-parameter, which can be changed by the
                    // called function. In C this could be solved by passing
                    // a pointer to s. In Java, a class variable must be made
                    // of this, which only is used in the context of fxshfr()
                    // and methods called from that.
    
    /*  Computes up to L2 fixed shift k-polynomials,
     *  testing for convergence in the linear or quadratic
     *  case. Initiates one of the variable shift
     *  iterations and returns with the number of zeros
     *  found.
     */
    private int fxshfr(int l2) {
        // double ui, vi;    // Not used because of reasons explained above.
        // double s;         // Not used because of reasons explained above.
        double betas, betav, oss, ovv, ss, vv, ts, tv;
        double ots = 0, otv = 0, tvv, tss;
        int type, i;
        boolean vtry, stry;

        betav = 0.25;
        betas = 0.25;
        oss = sr;
        ovv = v;
        /*  Evaluate polynomial by synthetic division. */
        quadsd_ab(n, u, v, p, qp);
        type = calcsc();
        for (int j = 1; j <= l2; j++) {
            /*  Calculate next k polynomial and estimate v. */
            nextk(type);
            type = calcsc();
            newest(type);    // Fills ui and vi with values.
            vv = vi;
            /*  Estimate s. */
            ss = 0.0;
            if (k[n - 1] != 0.0) {
                ss = -p[n] / k[n - 1];
            }
            tv = 1.0;
            ts = 1.0;
            if (j != 1 && type != 3) {

                /*  Compute relative measures of convergence of s and v sequences. */
                if (vv != 0.0) {
                    tv = Math.abs((vv - ovv) / vv);
                }
                if (ss != 0.0) {
                    ts = Math.abs((ss - oss) / ss);
                }
                /*  If decreasing, multiply two most recent convergence measures. */
                tvv = 1.0;
                if (tv < otv) {
                    tvv = tv * otv;
                }
                tss = 1.0;
                if (ts < ots) {
                    tss = ts * ots;
                }
                /*  Compare with convergence criteria. */
                boolean vpass = (tvv < betav);
                boolean spass = (tss < betas);
                if (spass || vpass) {
                    /*  At least one sequence has passed the convergence test.
                     *  Store variables before iterating.
                     */
                    double svu = u;
                    double svv = v;
                    for (i = 0; i < n; i++) {
                        svk[i] = k[i];
                    }
                    s = ss;
                    /*  Choose iteration according to the fastest converging
                     *  sequence.
                     */
                    vtry = false;
                    stry = false;
                    int state = 20;
                    if (spass && (!vpass || tss < tvv)) {
                        state = 40;
                    }
                    while (state != 70) {
                        switch (state) {
                            case 20: {   // Iteration for quadratic factor.
                                int nz = quadit();    // Beware: Changes ui and vi!
                                if (nz > 0) {
                                    return nz;
                                }
                                //  Quadratic iteration has failed. Flag that it has
                                //  been tried and decrease the convergence criterion.
                                vtry = true;
                                betav *= 0.25;
                                //  Try linear iteration if it has not been tried and
                                //  the S sequence is converging.
                                if (!stry && spass) {
                                    for (i = 0; i < n; i++) {
                                        k[i] = svk[i];
                                    }
                                    state = 40;
                                }
                                else {
                                    state = 50;
                                }
                                break;
                            }
                            case 40: {   // Iteration for linear factor.   
                                int nz = realit();   // Beware: changes s!
                                if (nz > 0) {
                                    return nz;
                                }
                                //  Linear iteration has failed. Flag that it has been
                                //  tried and decrease the convergence criterion.
                                stry = true;
                                betas *= 0.25;
                                if (nz == -1) {
                                    // nz == -1, almost double real zero, 
                                    // we attempt quadratic iteration.
                                    ui = -(s + s);
                                    vi = s * s;
                                    state = 20;
                                }
                                else  {  // nz == 0
                                    state = 50;
                                } 
                                break;
                            }
                            case 50: {   // Restoration of variables.
                                u = svu;
                                v = svv;
                                for (i = 0; i < n; i++) {
                                    k[i] = svk[i];
                                }
                                // Try quadratic iteration if it has not been tried
                                // and the V sequence is converging.
                                if (vpass && !vtry) {
                                    state = 20;
                                } 
                                else {
                                    // Recompute QP and scalar values to continue the
                                    // second stage.
                                    quadsd_ab(n, u, v, p, qp);
                                    type = calcsc();
                                    state = 70;
                                }
                                break;
                            }
                            default:   // Should never be reached, if this is reached, then internal error!
                                throw new RuntimeException("RPoly internal error (fxshfr:state=" + state + "), contact developer of this code!");
                        }
                    }
                }
            }

            ovv = vv;
            oss = ss;
            otv = tv;
            ots = ts;
        }
        
        // No convergence.
        return 0;
    }
    
    
    
    
    
    /*  Variable-shift k-polynomial iteration for a
     *  quadratic factor converges only if the zeros are
     *  equimodular or nearly so.
     *  This method uses ui and vi as inputs and also uses
     *  these as outputs.
     *  nz - number of zeros found.
     */
    private int quadit() {
        double mp, omp = 0, ee, relstp = 0, t, zm;
        int type, i, j;
        boolean tried;

        tried = false;
        u = ui;
        v = vi;
        j = 0;

        /*  Main loop. */
        int nOK = 0;
        while (true) {
            quad(1.0, u, v);   // Returns roots through szr, szi, lzr, lzi
            
            //  Return if roots of the quadratic are real and not
            //  close to multiple or nearly equal and of opposite sign.
            if (Math.abs(Math.abs(szr) - Math.abs(lzr)) > 0.01 * Math.abs(lzr)) {
                return 0;
            }
            //  Evaluate polynomial by quadratic synthetic division.
            quadsd_ab(n, u, v, p, qp);
            mp = Math.abs(a - szr * b) + Math.abs(szi * b);
            //  Compute a rigorous bound on the rounding error in evaluating p.
            zm = Math.sqrt(Math.abs(v));
            ee = 2.0 * Math.abs(qp[0]);
            t = -szr * b;
            for (i = 1; i < n; i++) {
                ee = ee * zm + Math.abs(qp[i]);
            }
            ee = ee * zm + Math.abs(a + t);
            ee = (5.0 * MRE + 4.0 * ARE) * ee
                    - (5.0 * MRE + 2.0 * ARE) * (Math.abs(a + t) + Math.abs(b) * zm)
                    + 2.0 * ARE * Math.abs(t);
            //  Iteration has converged sufficiently if the
            //  polynomial value is less than 20 times this bound.
            
            // WILCO: Added nOK logic. Only after two iterations
            // within the error bound we accept the root. The
            // additional iteration makes the returned root more
            // accurate in nearly all cases. This strategy removes
            // another known issue of Algorithm 493. The original
            // Fortran code sometimes returns a set of zeros, which
            // is isolated, but which is quite inaccurate, sometimes
            // only to 2 or 3 decimal places. The extra iteration
            // makes this issue much less severe, although in some
            // rare cases, where roots are present in clusters, there
            // still may be such cases.
            if (mp <= 20.0 * ee) {
                nOK++;
                if (nOK == 2) {
                    return 2;   // Converged to two zeros.
                }
            } 
            else {
                nOK = 0;
            }
            j++;
            //  Stop iteration after 20 steps.
            if (j > 20) {
                return 0;  // No convergence.
            }
            if (j >= 2) {
                if (relstp > 0.01 || mp < omp || tried) {

                } 
                else {
                    //  A cluster appears to be stalling the convergence.
                    //  Five fixed shift steps are taken with a u,v close
                    //  to the cluster.
                    if (relstp < ETA) {
                        relstp = ETA;
                    }
                    relstp = Math.sqrt(relstp);
                    u = u - u * relstp;
                    v = v + v * relstp;
                    quadsd_ab(n, u, v, p, qp);
                    for (i = 0; i < 5; i++) {
                        type = calcsc();
                        nextk(type);
                    }
                    tried = true;
                    j = 0;   // reset loop counter
                }
            }

            omp = mp;

            //  Calculate next k polynomial and new u and v.
            type = calcsc();
            nextk(type);
            type = calcsc();
            newest(type);
            //  If vi is zero the iteration is not converging.
            if (vi == 0.0) {
                return 0;
            }
            relstp = Math.abs((vi - v) / vi);
            u = ui;
            v = vi;
        }
    }
    
    
    
    
    
    /*  Variable-shift H polynomial iteration for a real zero.
     *  return nz, number of zeros found
     *  
     *  sss - starting iterate
     *  In the Java version, the argument sss is removed, it is
     *  connected to the variable s from the calling environment.
     *  iflag - flag to indicate a pair of zeros near real axis, used
     *          as output variable in original code.
     *  In the Java version, iflag is removed and the return value nz
     *  is set to -1 to indicate that a double zero is encountered.
     */
    private int realit(/* DoublePar sss, IntPar iflag */) {
        double t = 0, sss;
        double omp = 0;
        int j;

        sss = s;
        int nz = 0;
        j = 0;
        /*  Main loop */
        int nOK = 0;
        while (true) {
            double pv = p[0];
            /*  Evaluate p at s. */
            qp[0] = pv;
            for (int i = 1; i <= n; i++) {
                pv = pv * sss + p[i];
                qp[i] = pv;
            }
            double mp = Math.abs(pv);
            //  Compute a rigorous bound on the error in evaluating p.
            double ms = Math.abs(sss);
            double ee = (MRE / (ARE + MRE)) * Math.abs(qp[0]);
            for (int i = 1; i <= n; i++) {
                ee = ee * ms + Math.abs(qp[i]);
            }
            //  Iteration has converged sufficiently if the polynomial
            //  value is less than 20 times this bound.

            // WILCO: Added nOK logic. Only after two iterations
            // within the error bound we accept the root. The
            // additional iteration makes the returned root more
            // accurate in nearly all cases.
            // See comment in quadit() above!
            if (mp <= 20.0 * ((ARE + MRE) * ee - MRE * mp)) {
                nOK++;
                if (nOK == 2) {
                    nz = 1;
                    szr = sss;
                    szi = 0.0;
                    return nz;
                }
            } 
            else {
                nOK = 0;
            }
            j++;
            //  Stop iteration after 10 steps.
            if (j > 10) {
                return nz;
            }
            if (j >= 2) {
                if (Math.abs(t) > 0.001 * Math.abs(sss - t) || mp < omp) {

                } 
                else {
                    //  A cluster of zeros near the real axis has been
                    //  encountered. Return with iflag set to initiate a
                    //  quadratic iteration.
                    nz = -1;     // Indicate that quadratic iteration is needed.
                    s = sss;     // In the original code, the var parameter sss is modified,
                                 // here we set the class variable s, because in the original
                                 // code this function is called with a reference to s.
                    return nz;
                }
            }
            //  Return if the polynomial value has increased significantly.
            omp = mp;
            //  Compute t, the next polynomial, and the new iterate.
            double kv = k[0];
            qk[0] = kv;
            for (int i = 1; i < n; i++) {
                kv = kv * sss + k[i];
                qk[i] = kv;
            }
            if (Math.abs(kv) <= Math.abs(k[n - 1]) * 10.0 * ETA) {         // HVE n -> n-1
                //  Use unscaled form.
                k[0] = 0.0;
                for (int i = 1; i < n; i++) {
                    k[i] = qk[i - 1];
                }
            } else {
                //  Use the scaled form of the recurrence if the value
                //  of k at s is nonzero.
                t = -pv / kv;
                k[0] = qp[0];
                for (int i = 1; i < n; i++) {
                    k[i] = t * qk[i - 1] + qp[i];
                }
            }
            kv = k[0];
            for (int i = 1; i < n; i++) {
                kv = kv * sss + k[i];
            }
            t = 0.0;
            if (Math.abs(kv) > Math.abs(k[n - 1] * 10.0 * ETA)) {
                t = -pv / kv;
            }
            sss += t;
        }
    }

    
    
    
    
    /*  This routine calculates scalar quantities used to
     *  compute the next k polynomial and new estimates of
     *  the quadratic coefficients.
     *  type - integer variable set here indicating how the
     *  calculations are normalized to avoid overflow.
     */
    private int calcsc() {
        int type;
        /*  Synthetic division of k by the quadratic 1,u,v */
        quadsd_cd(n - 1, u, v, k, qk);
        if (Math.abs(c) <= Math.abs(k[n - 1] * 100.0 * ETA)
                && Math.abs(d) <= Math.abs(k[n - 2] * 100.0 * ETA)) {
            type = 3;
            /*  Type=3 indicates the quadratic is almost a factor of k. */
            return type;
        }

        if (Math.abs(d) < Math.abs(c)) {
            type = 1;
            /*  Type=1 indicates that all formulas are divided by c. */
            e = a / c;
            f = d / c;
            g = u * e;
            h = v * b;
            a3 = a * e + (h / c + g) * b;
            a1 = b - a * (d / c);
            a7 = a + g * d + h * f;
            return type;
        }
        type = 2;
        /*  Type=2 indicates that all formulas are divided by d. */
        e = a / d;
        f = c / d;
        g = u * b;
        h = v * b;
        a3 = (a + g) * e + h * (b / d);
        a1 = b * f - a;
        a7 = (f + u) * a + h;
        return type;
    }

    
    
    
    
    /*  Computes the next k polynomials using scalars 
     *  computed in calcsc.
     */
    private void nextk(int type) {
        double temp;
        int i;

        if (type == 3) {
            /*  Use unscaled form of the recurrence if type is 3. */
            k[0] = 0.0;
            k[1] = 0.0;
            for (i = 2; i < n; i++) {
                k[i] = qk[i - 2];
            }
            return;
        }
        temp = a;
        if (type == 1) {
            temp = b;
        }
        if (Math.abs(a1) <= Math.abs(temp) * ETA * 10.0) {
            /*  If a1 is nearly zero then use a special form of the
             *  recurrence.
             */
            k[0] = 0.0;
            k[1] = -a7 * qp[0];
            for (i = 2; i < n; i++) {
                k[i] = a3 * qk[i - 2] - a7 * qp[i - 1];
            }
        }
        else {
            /*  Use scaled form of the recurrence. */
            a7 /= a1;
            a3 /= a1;
            k[0] = qp[0];
            k[1] = qp[1] - a7 * qp[0];
            for (i = 2; i < n; i++) {
                k[i] = a3 * qk[i - 2] - a7 * qp[i - 1] + qp[i];
            }
        }
    }

    
    
    
    
    
    /*  Compute new estimates of the quadratic coefficients
     *  using the scalars computed in calcsc.
     */
    private void newest(int type) {
        double a4, a5, b1, b2, c1, c2, c3, c4, temp;

        /* Use formulas appropriate to setting of type. */
        if (type == 3) {
            /*  If type=3 the quadratic is zeroed. */
            ui = 0.0;
            vi = 0.0;
            return;
        }
        if (type == 2) {
            a4 = (a + g) * f + h;
            a5 = (f + u) * c + v * d;
        } else {
            a4 = a + u * b + h * f;
            a5 = c + (u + v * f) * d;
        }
        /*  Evaluate new quadratic coefficients. */
        b1 = -k[n - 1] / p[n];
        b2 = -(k[n - 2] + b1 * p[n - 1]) / p[n];
        c1 = v * b2 * a1;
        c2 = b1 * a7;
        c3 = b1 * b1 * a3;
        c4 = c1 - c2 - c3;
        temp = a5 + b1 * a4 - c4;
        if (temp == 0.0) {
            ui = 0.0;
            vi = 0.0;
            return;
        }
        ui = u - (u * (c3 + c2) + v * (b1 * a1 + b2 * a7)) / temp;
        vi = v * (1.0 + c4 / temp);
    }

    
    
    
    
    /*  Divides p by the quadratic 1,u,v placing the quotient
     *  in q and the remainder in class variables a,b.
     */
    private void quadsd_ab(int nn, double u, double v, double[] p, double[] q) {
        double cc;
        int i;

        b = p[0];
        q[0] = b;
        a = p[1] - b * u;
        q[1] = a;
        for (i = 2; i <= nn; i++) {
            cc = p[i] - a * u - b * v;
            q[i] = cc;
            b = a;
            a = cc;
        }
    }
    /*  Divides p by the quadratic 1,u,v placing the quotient
     *  in q and the remainder in class variables c,d.
     */
    private void quadsd_cd(int nn, double u, double v, double[] p, double[] q) {
        double cc;
        int i;

        d = p[0];
        q[0] = d;
        c = p[1] - d * u;
        q[1] = c;
        for (i = 2; i <= nn; i++) {
            cc = p[i] - c * u - d * v;
            q[i] = cc;
            d = c;
            c = cc;
        }
    }

    
    
    
    
    
    /*  Calculate the zeros of the quadratic a*z^2 + b1*z + c.
     *  The quadratic formula, modified to avoid overflow, is used 
     *  to find the larger zero if the zeros are real and both
     *  are complex. The smaller real zero is found directly from 
     *  the product of the zeros c/a.
     */
    private void quad(double a, double b1, double c) {
        double b, d, e;

        if (a == 0.0) {         /* less than two roots */
            szr = (b1!=0.0) ? -c/b1 : 0.0;
            lzr = 0.0;
            szi = 0.0;
            lzi = 0.0;
            return;
        }
        
        if (c == 0.0) {         /* one real root, one zero root */
            szr = 0.0;
            lzr = -b1/a;
            szi = 0.0;
            lzi = 0.0;
            return;
        }
        
        /* Compute square root of abs(discriminant) avoiding overflow. */
        b = b1 / 2.0;
        double ab = Math.abs(b);
        double ac = Math.abs(c);
        if (ab < ac) {
            e = (c<0.0) ? -a : a;
            e = b * (b / ac) - e;
            d = Math.sqrt(Math.abs(e)) * Math.sqrt(ac);
        } 
        else {
            e = 1.0 - (a / b) * (c / b);
            d = Math.sqrt(Math.abs(e)) * ab;
        }
        if (e < 0.0) {   
            /* complex conjugate zeros */
            lzr = szr = -b / a;
            szi = Math.abs(d / a);
            lzi = -szi;
        } 
        else {
            /* real zeros, compute zero with largest magnitude first */
            if (b >= 0.0) {
                d = -d;
            }
            lzr = (-b + d) / a;
            szr = 0.0;
            if (lzr != 0.0) {
                szr = (c / lzr) / a;
            }
            szi = 0.0;
            lzi = 0.0;
        }
    }
}
