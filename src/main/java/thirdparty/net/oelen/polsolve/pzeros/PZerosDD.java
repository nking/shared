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
*/
package thirdparty.net.oelen.polsolve.pzeros;

// PolSolve code for DoubleDouble precision (105 bits)

import thirdparty.net.oelen.polarith.DoubleComplex;
import thirdparty.net.oelen.polarith.DoubleDouble;



/************************************************************************
*    NUMERICAL COMPUTATION OF THE ROOTS OF A POLYNOMIAL HAVING          *
*        COMPLEX COEFFICIENTS, BASED ON ABERTH'S METHOD.                *
*                      Version 1.4, June   1996                         *
*    (D. Bini, Dipartimento di Matematica, Universita' di Pisa)         *
*                         (bini@dm.unipi.it)                            *
*************************************************************************
* Work performed under the support of the ESPRIT BRA project 6846 POSSO *
*************************************************************************
***********         SUBROUTINES AND FUNCTIONS                 ***********
*************************************************************************
*  The following modules are listed:                                    *
*  POLZEROS  :  computes polynomial roots by means of Aberth's method   *
*    ABERTH  :  computes the Aberth correction                          *
*    NEWTON  :  computes p(x)/p'(x) by means of Ruffini-Horner's rule   *
*    START   :  Selects N starting points by means of Rouche's theorem  *
*    CNVEX   :  Computes the convex hull, used by START                 *
*    CMERGE  :  Used by CNVEX                                           *
*    LEFT    :  Used by CMERGE                                          *
*    RIGHT   :  Used by CMERGE                                          *
*    CTEST   :  Convexity test, Used by CMERGE                          *
*************************************************************************/

strictfp class PZerosDD {

    
    
    
    
    /************************************************************************
    *                             SUBROUTINE NEWTON                         *
    *************************************************************************
    * Compute  the Newton's correction, the inclusion radius (4) and checks *
    * the stop condition (3)                                                *
    *************************************************************************
    * Input variables:                                                      *
    *     N     : degree of the polynomial p(x)                             *
    *     POLY  : coefficients of the polynomial p(x)                       *
    *     APOLY : upper bounds on the backward perturbations on the         *
    *             coefficients of p(x) when applying Ruffini-Horner's rule  *
    *     APOLYR: upper bounds on the backward perturbations on the         *
    *             coefficients of p(x) when applying (2), (2')              *
    *     Z     : value at which the Newton correction is computed          *
    *     SMALL : the min positive real*8, SMALL=2**(-1022) for the IEEE.   *
    *************************************************************************
    * Output variables:                                                     *
    *     RADIUS: upper bound to the distance of Z from the closest root of *
    *             the polynomial computed according to (4).                 *
    *     CORR  : Newton's correction                                       *
    *     AGAIN : this variable is .true. if the computed value p(z) is     *
    *             reliable, i.e., (3) is not satisfied in Z. AGAIN is       *
    *             .false., otherwise.                                       *
     @param n
     @param again
     @param apoly
     @param apolyr
     @param poly
     @param z
     @param outIndex
     @param small
     @param eps
     @param radius
     @return 
    *************************************************************************/
    // Version for a polynomial with complex coefficients.
    private static DoubleComplex newton(int n, DoubleComplex[] poly, double[] apoly, double[] apolyr,
                                        DoubleComplex z, double small, double eps,
                                        /* out */ double[] radius, /* out */ boolean[] again,
                                        int outIndex) {
        int i;
        DoubleComplex p, p1;
        DoubleComplex corr;
        double ap, az;

        az = z.ComplexValue().abs();

        // If |z|<=1 then apply Ruffini-Horner's rule for p(z)/p'(z) 
        // and for the computation of the inclusion radius.
        // If |z|>1 then apply Ruffini-Horner's rule to the reversed polynomial
        // and use formula (2) for p(z)/p'(z). Analogously do for the inclusion
        // radius.
        if (az < 1.0) {
            p = poly[n];
            p1 = p;
            ap = apoly[n];
            for (i = n - 1; i >= 1; i--) {
                p = poly[i].add(p.mul(z));
                p1 = p.add(p1.mul(z));
                ap = apoly[i] + ap * az;
            }
            p = poly[0].add(p.mul(z));
            ap = apoly[0] + ap * az;
            corr = p.div(p1);

            double epsap = eps * ap;
            // double absp = p.abs().doubleValue();
            double absp = Math.abs(p.re.doubleValue()) + Math.abs(p.im.doubleValue());
            again[outIndex] = (absp > small + epsap);
            if (!again[outIndex]) {
                radius[outIndex] = n * (absp + epsap) / p1.abs().doubleValue();
            }
        } else // az >= 1.0
        {
            DoubleComplex zi, den, ppsp;
            double azi;

            zi = z.recip();
            azi = 1.0 / az;
            p = poly[0];
            p1 = p;
            ap = apolyr[n];
            for (i = n - 1; i >= 1; i--) {
                p = poly[n - i].add(p.mul(zi));
                p1 = p.add(p1.mul(zi));
                ap = apolyr[i] + ap * azi;
            }
            p = poly[n].add(p.mul(zi));
            ap = apolyr[0] + ap * azi;

            ppsp = p.mul(z).div(p1);
            den = ppsp.mul(n).sub1();
            corr = z.mul(ppsp.div(den));
            double epsap = eps * ap;
            // double absp = p.abs().doubleValue();
            double absp = Math.abs(p.real().doubleValue()) + Math.abs(p.imag().doubleValue());
            again[outIndex] = (absp > small + epsap);
            if (!again[outIndex]) {
                radius[outIndex] = ppsp.abs().doubleValue() + (epsap * az) / p1.abs().doubleValue();
                radius[outIndex] *= n / den.abs().doubleValue();
                radius[outIndex] *= az;
            }
        }
        
        return corr;
    }
    
    
    

    // Version for a polynomial with real coefficients.
    private static DoubleComplex newton(int n, DoubleDouble[] poly, double[] apoly, double[] apolyr,
                                        DoubleComplex z, double small, double eps,
                                        /* out */ double[] radius, /* out */ boolean[] again,
                                        int outIndex) {
        int i;
        DoubleComplex p, p1;
        double ap, az;

        az = z.abs().doubleValue();
        
        DoubleComplex corr;

        // Use a method like Horner's method, but developed especially for
        // evaluation of real polynomials for complex values.
        // This method requires only half amount of the number of real
        // multiplications as standard Horner's rule.
        // BIT 5 (1965), 142, see also G. Goertzel, AMM 65 (1958), 34-35
        // This method only works for polynomials of degree >= 2.

        // If |z|<=1 then apply modified Horner's rule for p(z)/p'(z) 
        // and for the computation of the inclusion radius.
        // If |z|>1 then apply modified Horner's rule to the reversed polynomial
        // and use formula (2) for p(z)/p'(z). Analogously do for the inclusion
        // radius.
        if (az < 1.0) {
            DoubleDouble r = z.real().mulPowerOf2(2);
            DoubleDouble s = z.sqrabs();
            DoubleDouble aa = poly[n];
            DoubleDouble bb = poly[n - 1];

            // Compute p(z) and |p|(|z|)
            ap = apoly[n - 1] + az * apoly[n];
            for (i = n - 2; i >= 0; i--) {
                // Next step for polynomial.
                DoubleDouble tmp = bb.add(r.mul(aa));
                bb = poly[i].sub(s.mul(aa));
                aa = tmp;

                // Next step for |p|(|z|)
                ap = apoly[i] + ap * az;
            }
            p = z.mul(aa).add(bb);

            // Compute p'(z)
            aa = poly[n].mul(n);
            bb = poly[n - 1].mul(n - 1);
            for (i = n - 2; i >= 1; i--) {
                DoubleDouble tmp = bb.add(r.mul(aa));
                bb = poly[i].mul(i).sub(s.mul(aa));
                aa = tmp;
            }
            p1 = z.mul(aa).add(bb);

            double epsap = eps * ap;
            corr = p.div(p1);
            //double absp = p.abs().doubleValue();
            double absp = Math.abs(p.real().doubleValue()) + Math.abs(p.imag().doubleValue());
            again[outIndex] = (absp > small + epsap);
            if (!again[outIndex]) {
                radius[outIndex] = n * (absp + epsap) / p1.abs().doubleValue();
            }
        } 
        else // az >= 1.0
        {
            DoubleComplex den, ppsp;

            DoubleComplex zi = z.recip();
            double azi = 1.0 / az;

            DoubleDouble r = zi.real().mulPowerOf2(2);
            DoubleDouble s = zi.sqrabs();
            DoubleDouble aa = poly[0];
            DoubleDouble bb = poly[1];

            // Compute q(zi) and |q|(|zi|), where
            // q equals p with reversed coefficients.
            ap = apolyr[n - 1] + azi * apolyr[n];
            for (i = n - 2; i >= 0; i--) {
                // Next step for polynomial.
                DoubleDouble tmp = bb.add(r.mul(aa));
                bb = poly[n - i].sub(s.mul(aa));
                aa = tmp;

                // Next step for |q|(|zi|)
                ap = apolyr[i] + ap * azi;
            }
            p = zi.mul(aa).add(bb);

            // Compute q'(zi)
            aa = poly[0].mul(n);
            bb = poly[1].mul(n - 1);
            for (i = n - 2; i >= 1; i--) {
                DoubleDouble tmp = bb.add(r.mul(aa));
                bb = poly[n - i].mul(i).sub(s.mul(aa));
                aa = tmp;
            }
            p1 = zi.mul(aa).add(bb);

            ppsp = p.mul(z).div(p1);
            den = ppsp.mul(n).sub1();
            corr = z.mul(ppsp.div(den));
            double epsap = eps * ap;
            // double absp = p.abs().doubleValue();
            double absp = Math.abs(p.real().doubleValue()) + Math.abs(p.imag().doubleValue());
            again[outIndex] = (absp > small + epsap);
            if (!again[outIndex]) {
                radius[outIndex] = (ppsp.abs().doubleValue() + (epsap * az) / p1.abs().doubleValue());
                radius[outIndex] *= n / den.abs().doubleValue();
                radius[outIndex] *= az;
            }
        }
        
        return corr;
    }

    
    
    
    
    /************************************************************************
    *                             SUBROUTINE ABERTH                         *
    *************************************************************************
    * Compute  the Aberth correction. To save time, the reciprocation of    *
    * ROOT(J)-ROOT(I) could be performed in single precision (complex*8)    *
    * In principle this might cause overflow if both ROOT(J) and ROOT(I)    *
    * have too small moduli.                                                *
    *************************************************************************
    * Input variables:                                                      *
    *     N     : degree of the polynomial                                  *
    *     ROOT  : vector containing the current approximations to the roots *
    *     J     : index of the component of ROOT with respect to which the  *
    *             Aberth correction is computed                             *
    *************************************************************************
    * Output variable:                                                      *
    *     ABCORR: Aberth's correction (compare (1))                         *
     @param n
     @param j
     @param root
     @return 
    *************************************************************************/
    
    /* The pseudocode below is left here for documentation purposes, the real
     * code is further below. By doing the aberth correction in normal double
     * precision, the algorithm hardly suffers but the speedup is considerable
     * (more than a factor 2 speedup for quad double arithmatic in a typical
     * case). Only the difference between the roots is computed, using full
     * precision (this may be necessary in case of (nearly) equal roots), but
     * the 1/(z-zj) correction can be computed in standard precision.

       static complex aberth(int n, int j, complex[] root) { 
           int i; 
           complex zj;
     
           complex abcorr = 0.0; 
           zj = root[j]; 
           for (i=0; i<n; i++) {
               if (i != j) {
                   complex z = zj - root[i];
                   abcorr += 1.0 / z;
               }
           }
           return abcorr;
       }

     */
    
    private static DoubleComplex aberth(int n, int j, DoubleComplex root[]) {
        DoubleDouble zj_re = root[j].re;
        DoubleDouble zj_im = root[j].im;

        double corr_re = 0.0;
        double corr_im = 0.0;
        for (int i = 0; i < n; i++) {
            if (i != j) {
                double z_re = zj_re.sub(root[i].re).doubleValue();
                double z_im = zj_im.sub(root[i].im).doubleValue();
                double absr = (z_re > 0) ? z_re : -z_re;
                double absi = (z_im > 0) ? z_im : -z_im;
                if (absr > absi) {
                    double r = z_im / z_re;
                    double den = z_re + r * z_im;
                    corr_re += 1.0/den;
                    corr_im -= r/den;
                } 
                else {
                    double r = z_re / z_im;
                    double den = z_im + r * z_re;
                    corr_re += r/den;
                    corr_im -= 1.0/den;
                }
            }
        }

        return new DoubleComplex(corr_re, corr_im);  // Aberth correction
    }
    
    
    

    /************************************************************************
    *                             SUBROUTINE START                          *
    *************************************************************************
    * Compute  the starting approximations of the roots                     *
    *************************************************************************
    * Input variables:                                                      *
    *     N     :  number of the coefficients of the polynomial             *
    *     A     :  moduli of the coefficients of the polynomial             *
    *     SMALL : the min positive real*8, SMALL=2**(-1074) for the IEEE.   *
    *     BIG   : the max real*8, BIG=2**1023 for the IEEE standard.        *
    * Output variables:                                                     *
    *     Y     :  starting approximations                                  *
    *     RADIUS:  if a component is -1 then the corresponding root has a   *
    *              too big or too small modulus in order to be represented  *
    *              as double float with no overflow/underflow               *
    *     NZ    :  number of roots which cannot be represented without      *
    *              overflow/underflow                                       *
    * Auxiliary variables:                                                  *
    *     H     :  needed for the computation of the convex hull            *
    *************************************************************************
    * This routine selects starting approximations along circles center at  *
    * 0 and having suitable radii. The computation of the number of circles *
    * and of the corresponding radii is performed by computing the upper    *
    * convex hull of the set (i,log(A(i))), i=1,...,n+1.                    *
     @param n
     @param a
     @param y
     @param radius
     @param big
     @return 
    *************************************************************************/
    private static int start(int n, double[] a, DoubleComplex[] y, double[] radius, double big) {
        int i, iold, nzeros, j, jj;
        double r, th, ang, temp;
        final double sigma = 0.7;
        final double pi2 = 6.2831853071796;
        double xbig = Math.log(big);
        double xsmall = -xbig;
        boolean[] h = new boolean[n + 1];
        int nz;   // number of zeros which cannot be represented without under/overflow.

        nz = 0;

        // Compute the logarithm A(I) of the moduli of the coefficients of
        // the polynomial and then the upper convex hull of the set (A(I),I)
        for (i = 0; i <= n; i++) {
            if (a[i] != 0.0) {
                a[i] = Math.log(a[i]);
            } else {
                a[i] = -1e30;   // -infinity
            }
        }
        Convex.cnvex(n, a, h);

        // Given the upper convex hull of the set (A(I),I) compute the moduli
        // of the starting approximations by means of Rouche's theorem
        iold = 0;
        th = pi2 / n;
        for (i = 1; i <= n; i++) {
            if (h[i]) {
                nzeros = i - iold;
                temp = (a[iold] - a[i]) / nzeros;
                boolean underflow = false, overflow = false;
                if (temp <= xsmall) {
                    // Radii are too small, underflow situation.
                    underflow = true;
                    r = 0.0;
                    nz += nzeros;
                } 
                else if (temp >= xbig) {
                    // Overflow situation.
                    overflow = true;
                    r = big;
                    nz += nzeros;
                } 
                else {
                    // In range of floating point arithmetic.
                    r = Math.exp(temp);
                }

                ang = pi2 / nzeros;
                for (j = iold; j < i; j++) {
                    jj = j - iold;
                    if (underflow || overflow) {
                        radius[j] = -1.0;
                        y[j] = new DoubleComplex(r);
                    } 
                    else {
                        double re = r * Math.cos(ang * jj + th * i + sigma);
                        double im = r * Math.sin(ang * jj + th * i + sigma);
                        y[j] = new DoubleComplex(re, im);
                    }
                }

                // Make the current vertex the basis for the
                // new set of roots.
                iold = i;
            }
        }
        
        return nz;
    }
    
    

    /************************************************************************
    *                             SUBROUTINE SORT                           *
    *************************************************************************
    *   SORT  the vector X, according to nonincreasing real parts,          *
    *   the same permutation is applied to vectors Y and E.                 *
     @param n
     @param x
     @param y
     @param e
    *************************************************************************/
    private static void sort(int n, DoubleComplex[] x, double[] y, boolean[] e) {
        for (int k = 0; k < n - 1; k++) {
            DoubleDouble amax = x[k].real();
            int imax = k;
            for (int i = k + 1; i < n; i++) {
                DoubleDouble rxi = x[i].real();
                if (amax.compareTo(rxi) < 0) {
                    amax = rxi;
                    imax = i;
                }
            }

            DoubleComplex temp = x[k];
            x[k] = x[imax];
            x[imax] = temp;

            double yt = y[k];
            y[k] = y[imax];
            y[imax] = yt;

            boolean et = e[k];
            e[k] = e[imax];
            e[imax] = et;
        }
    }

    
    
    
    
    /************************************************************************
    *********************** SUBROUTINE POLZEROS *****************************
    *************************************************************************
    *                        GENERAL COMMENTS                               *
    *************************************************************************
    *  This routine approximates the roots of   the  polynomial             *
    *  p(x)=a(n)x^n+a(n-1)x^(n-1)+...+a(0), a(j)=cr(j)+I ci(j), I**2=-1,    *
    *  where a(0) and a(n) are nonzero.                                     *
    *  The coefficients are complex*16 numbers. The routine is fast, robust *
    *  against overflow, and allows to deal with polynomials of any degree. *
    *  Overflow situations are very unlikely and may occurr if there exist  *
    *  simultaneously coefficients of moduli close to BIG and close to      *
    *  SMALL, i.e., the greatest and the smallest positive real*8 numbers,  *
    *  respectively. In this limit situation the program outputs a warning  *
    *  message. The computation can be speeded up by performing some side   *
    *  computations in single precision, thus slightly reducing the         *
    *  robustness of the program (see the comments in the routine ABERTH).  *
    *  Besides a set of approximations to the roots, the program delivers a *
    *  set of a-posteriori error bounds which are guaranteed in the most    *
    *  part of cases. In the situation where underflow does not allow to    *
    *  compute a guaranteed bound, the program outputs a warning message    *
    *  and sets the bound to 0. In the situation where the root cannot be   *
    *  represented as a complex*16 number the error bound is set to -1.     *
    *************************************************************************
    *  The computation is performed by means of Aberth's method             *
    *  according to the formula                                             *
    *           x(i)=x(i)-newt/(1-newt*abcorr), i=1,...,n             (1)   *
    *  where newt=p(x(i))/p'(x(i)) is the Newton correction and abcorr=     *
    *  =1/(x(i)-x(1))+...+1/(x(i)-x(i-1))+1/(x(i)-x(i+1))+...+1/(x(i)-x(n)) *
    *  is the Aberth correction to the Newton method.                       *
    *************************************************************************
    *  The value of the Newton correction is computed by means of the       *
    *  synthetic division algorithm (Ruffini-Horner's rule) if |x|<=1,      *
    *  otherwise the following more robust (with respect to overflow)       *
    *  formula is applied:                                                  *
    *                    newt=1/(n*y-y**2 R'(y)/R(y))                 (2)   *
    *  where                                                                *
    *                    y=1/x                                              *
    *                    R(y)=a(0)*y**n+...+a(n-1)*y+a(n)            (2')   *
    *  This computation is performed by the routine NEWTON.                 *
    *************************************************************************
    *  The starting approximations are complex numbers that are             *
    *  equispaced on circles of suitable radii. The radius of each          *
    *  circle, as well as the number of roots on each circle and the        *
    *  number of circles, is determined by applying Rouche's theorem        *
    *  to the functions a(k)*x**k and p(x)-a(k)*x**k, k=0,...,n.            *
    *  This computation is performed by the routine START.                  *
    *************************************************************************
    *                              STOP CONDITION                           *
    *************************************************************************
    * If the condition                                                      *
    *                     |p(x(j))|<EPS s(|x(j)|)                      (3)  *
    * is satisfied,    where      s(x)=s(0)+x*s(1)+...+x**n * s(n),         *
    * s(i)=|a(i)|*(1+3.8*i),  EPS is the machine precision (EPS=2**-53      *
    * for the IEEE arithmetic), then the approximation x(j) is not updated  *
    * and the subsequent iterations (1)  for i=j are skipped.               *
    * The program stops if the condition (3) is satisfied for j=1,...,n,    *
    * or if the maximum number NITMAX of  iterations   has   been reached.  *
    * The condition (3) is motivated by a backward rounding error analysis  *
    * of the Ruffini-Horner rule, moreover the condition (3) guarantees     *
    * that the computed approximation x(j) is an exact root of a slightly   *
    * perturbed polynomial.                                                 *
    *************************************************************************
    *             INCLUSION DISKS, A-POSTERIORI ERROR BOUNDS                *
    *************************************************************************
    * For each approximation x of a root, an a-posteriori absolute error    *
    * bound r is computed according to the formula                          *
    *                   r=n(|p(x)|+EPS s(|x|))/|p'(x)|                 (4)  *
    * This provides an inclusion disk of center x and radius r containing a *
    * root.                                                                 *
    *************************************************************************
    *************************************************************************
    **************       MEANING OF THE INPUT VARIABLES         *************
    *************************************************************************
    *************************************************************************
    *                                                                       *
    *  -- N     : degree of the polynomial.                                 *
    *  -- POLY  : complex vector of N+1 components, POLY(i) is the          *
    *           coefficient of x**(i), i=0,1,...,N of the polynomial p(x)   *
    *  -- EPS   : machine precision of the floating point arithmetic used   *
    *            by the computer, EPS=2**(-53)  for the IEEE standard.      *
    *  -- BIG   : the max real*8, BIG=2**1024 for the IEEE standard.        *
    *  -- SMALL : the min positive real*8, SMALL=2**(-1022) for the IEEE.   *
    *  -- NITMAX: the max number of allowed iterations.                     *
    *************************************************************************
    *************************************************************************
    **************      MEANING OF THE OUTPUT VARIABLES         *************
    *************************************************************************
    *************************************************************************
    *  ROOT   : complex vector of N components, containing the              *
    *           approximations to the roots of p(x).                        *
    *  RADIUS : real vector of N components, containing the error bounds to *
    *           the approximations of the roots, i.e. the disk of center    *
    *           ROOT(i) and radius RADIUS(i) contains a root of p(x), for   *
    *           i=1,...,N. RADIUS(i) is set to -1 if the corresponding root *
    *           cannot be represented as floating point due to overflow or  *
    *           underflow.                                                  *
    *  ERR    : vector of N components detecting an error condition;        *
    *           ERR(j)=.TRUE. if after NITMAX iterations the stop condition *
    *                         (3) is not satisfied for x(j)=ROOT(j);        *
    *           ERR(j)=.FALSE.  otherwise, i.e., the root is reliable,      *
    *                         i.e., it can be viewed as an exact root of a  *
    *                         slightly perturbed polynomial.                *
    *           The vector ERR is used also in the routine convex hull for  *
    *           storing the abscissae of the vertices of the convex hull.   *
    *  ITER   : number of iterations performed.                             *
    *************************************************************************
    *************************************************************************
    *************    MEANING OF THE AUXILIARY VARIABLES         *************
    *************************************************************************
    *************************************************************************
    *  APOLY  : real vector of N+1 components used to store the moduli of   *
    *           the coefficients of p(x) and the coefficients of s(x) used  *
    *           to test the stop condition (3).                             *
    *  APOLYR : real vector of N+1 components used to test the stop         *
    *           condition                                                   *
    *************************************************************************
    ******         WARNING:   2 is the output unit                    *******
     @param n
     @param init
     @param eps
     @param poly
     @param radius
     @param root
     @param nitmax
     @param big
     @param err
     @param iter
     @param needSort
     @param apolyr
     @param apoly
     @return 
    *************************************************************************/
    // Version for polynomial with real coefficients.
    static int polzeros(int n, DoubleDouble[] poly, double eps, double big,
                        int nitmax, DoubleComplex[] root, double[] radius, boolean[] err,
                        /* out */ int[] iter, boolean init, double[] apoly, double[] apolyr, boolean needSort) {
        int i;
        int nzeros;
        double amax;
        double small = 1.0 / big;

        if (init) {
            // Check consistency of data
            if (poly[n].isZero()) {
                return -1;
            }
            if (poly[0].isZero()) {
                return -2;
            }

            // Compute the moduli of the coefficients
            // and the largest of the absolute values.
            amax = 0.0;
            for (i = 0; i <= n; i++) {
                apoly[i] = poly[i].abs().doubleValue();
                if (amax < apoly[i]) {
                    amax = apoly[i];
                }
                apolyr[i] = apoly[i];
            }

            if (amax > big / (4 * n + 1)) {
                return -3;
            }

            if (n == 1) {
                // Linear equation, subsequent code is not called.
                root[0] = new DoubleComplex(poly[0].div(poly[1]).neg());
                err[0] = false;
                iter[0] = 1;
                radius[0] = eps * root[0].abs().doubleValue();
                return 0;
            }

            // Initialize
            for (i = 0; i < n; i++) {
                radius[i] = 0.0;
            }

            // Select the starting points. We use the array apolyr
            // as the array of absolute values of the polynomial
            // coefficients. The storage for apolyr is overwritten
            // by start(), but that is of no cencern. After calling
            // start(), the apoly and apolyr arrays are derived
            // from the apoly array values.
            nzeros = start(n, apolyr, root, radius, big);

            // Compute the coefficients of the backward-error polynomial.
            for (i = 0; i <= n; i++) {
                apolyr[n - i] = apoly[i] * (3.8 * (n - i) + 1.0);
                apoly[i] = apoly[i] * (3.8 * i + 1.0);
            }

            for (i = 0; i < n; i++) {
                err[i] = (radius[i] >= 0.0);
            }

            // Set number of iterations equal to 0.
            iter[0] = 0;
        } 
        else {
            nzeros = 0;
            for (i = 0; i < n; i++) {
                if (!err[i]) {
                    nzeros++;
                }
            }
        }

        // Starts Aberth's iterations.
        // Beware: The vector err[] is now used as a flags vector
        // again[] in order to determine whether we need to iterate
        // again or not.
        for (; iter[0] < nitmax; iter[0]++) {
            for (i = 0; i < n; i++) {
                if (err[i]) {
                    DoubleComplex corr = newton(n, poly, apoly, apolyr, root[i], small, eps, radius, err, i);
                    DoubleComplex abcorr = aberth(n, i, root);
                    root[i] = root[i].sub(corr.div(corr.mul(abcorr).sub1().neg()));
                    if (!err[i]) {
                        // We have found a new root, increment the total
                        // number of roots and advance to the next root.
                        nzeros++;
                        if (nzeros == n) {
                            if (needSort) {
                                sort(n, root, radius, err);
                            }
                            return 0;
                        }
                    }
                }
            }
        }

        return -4;
    }

    
    // Version for polynomial with complex coefficients.
    static int polzeros(int n, DoubleComplex[] poly, double eps, double big,
                        int nitmax, DoubleComplex[] root, double[] radius, boolean[] err,
                        /* out */ int[] iter, boolean init, double[] apoly, double[] apolyr, boolean needSort) {
        int i;
        int nzeros;
        double small = 1.0 / big;

        if (init) {
            // Check consistency of data
            if (poly[n].isZero()) {
                return -1;
            }
            if (poly[0].isZero()) {
                return -2;
            }

            // Compute the moduli of the coefficients
            // and the largest of the absolute values.
            double amax = 0.0;
            for (i = 0; i <= n; i++) {
                apoly[i] = poly[i].abs().doubleValue();
                if (amax < apoly[i]) {
                    amax = apoly[i];
                }
                apolyr[i] = apoly[i];
            }

            if (amax > big / (4 * n + 1)) {
                return -3;
            }

            if (n == 1) {
                // Linear equation, subsequent code is not called.
                root[0] = poly[0].div(poly[1]).neg();
                err[0] = false;
                iter[0] = 1;
                radius[0] = eps * root[0].abs().doubleValue();
                return 0;
            }

            // Initialize
            for (i = 0; i < n; i++) {
                radius[i] = 0.0;
            }

            // Select the starting points. We use the array apolyr
            // as the array of absolute values of the polynomial
            // coefficients. The storage for apolyr is overwritten
            // by start(), but that is of no cencern. After calling
            // start(), the apoly and apolyr arrays are derived
            // from the apoly array values.
            nzeros = start(n, apolyr, root, radius, big);

            // Compute the coefficients of the backward-error polynomial.
            for (i = 0; i <= n; i++) {
                apolyr[n - i] = apoly[i] * (3.8 * (n - i) + 1.0);
                apoly[i] = apoly[i] * (3.8 * i + 1.0);
            }

            for (i = 0; i < n; i++) {
                err[i] = (radius[i] >= 0.0);
            }

            // Set number of iterations equal to 0.
            iter[0] = 0;
        }
        else {
            nzeros = 0;
            for (i = 0; i < n; i++) {
                if (!err[i]) {
                    nzeros++;
                }
            }
        }

        // Starts Aberth's iterations.
        // Beware: The vector err[] is now used as a flags vector
        // again[] in order to determine whether we need to iterate
        // again or not.
        for (; iter[0] < nitmax; iter[0]++) {
            for (i = 0; i < n; i++) {
                if (err[i]) {
                    DoubleComplex corr = newton(n, poly, apoly, apolyr, root[i], small, eps, radius, err, i);
                    DoubleComplex abcorr = aberth(n, i, root);
                    root[i] = root[i].sub(corr.div(corr.mul(abcorr).sub1().neg()));
                    if (!err[i]) {
                        // We have found a new root, increment the total
                        // number of roots and advance to the next root.
                        nzeros++;
                        if (nzeros == n) {
                            if (needSort) {
                                sort(n, root, radius, err);
                            }
                            return 0;
                        }
                    }
                }
            }
        }

        return -4;
    }
    
}
