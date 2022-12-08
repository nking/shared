package thirdparty.net.oelen.polsolve.pzeros;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import junit.framework.TestCase;
import thirdparty.net.oelen.polarith.Complex;
import thirdparty.net.oelen.polarith.DoubleComplex;
import thirdparty.net.oelen.polarith.DoubleDouble;

public class PZerosTest extends TestCase {

    private DoubleDouble coefReal[];
    private DoubleComplex coefComplex[];

    // Constructor, used privately.
    private void run0(DoubleDouble[] coefs) {
        coefReal = coefs;
        coefComplex = null;
    }

    /**
     * Constructor, used for testing purposes.
     *
     @param testPolynomial A String array, as specified in the TestPolynomials
     * class. This is converted to a polynomial, useful for solving/testing.
     */
    public void run1(String[] testPolynomial) {
        List<DoubleDouble> roots_re = new ArrayList<>();
        List<DoubleDouble> roots_im = new ArrayList<>();
        if (testPolynomial[0].equals("RR")) {
            // Roots Real
            for (int i = 1; i < testPolynomial.length; i++) {
                String[] num = testPolynomial[i].split(";");
                roots_re.add(new DoubleDouble(num[0]));
                if (num.length >= 2) {
                    roots_re.add(new DoubleDouble(num[0]));
                    roots_im.add(new DoubleDouble(num[1]));
                    roots_im.add(new DoubleDouble(num[1]).neg());
                } else {
                    roots_im.add(DoubleDouble.ZERO);
                }
            }
            int degree = roots_re.size();
            coefReal = new DoubleDouble[degree + 1];
            coefComplex = null;
        } else if (testPolynomial[0].equals("RC")) {
            // Roots Complex
            for (int i = 1; i < testPolynomial.length; i++) {
                String[] num = testPolynomial[i].split(";");
                roots_re.add(new DoubleDouble(num[0]));
                roots_im.add(num.length >= 2 ? new DoubleDouble(num[1]) : DoubleDouble.ZERO);
            }
            int degree = roots_re.size();
            coefComplex = new DoubleComplex[degree + 1];
            coefReal = null;
        } else if (testPolynomial[0].equals("CR")) {
            // Coefs Real
            int ncoefs = testPolynomial.length - 1;
            coefReal = new DoubleDouble[ncoefs];
            for (int i = 1; i < testPolynomial.length; i++) {
                coefReal[i - 1] = new DoubleDouble(testPolynomial[i]);
            }
            coefComplex = null;
        } else if (testPolynomial[0].equals("CC")) {
            // Coefs Complex
            int ncoefs = testPolynomial.length - 1;
            coefReal = null;
            coefComplex = new DoubleComplex[ncoefs];
            for (int i = 1; i < testPolynomial.length; i++) {
                String[] num = testPolynomial[i].split(";");
                DoubleDouble re = new DoubleDouble(num[0]);
                DoubleDouble im = (num.length >= 2) ? new DoubleDouble(num[1]) : DoubleDouble.ZERO;
                coefComplex[i - 1] = new DoubleComplex(re, im);
            }
        } else {
            coefReal = null;
            coefComplex = null;
        }

        if (!roots_re.isEmpty()) {
            if (testPolynomial[0].equals("RR")) {
                coefReal[0] = DoubleDouble.ONE;

                for (int k = 0; k < roots_re.size(); k++) {
                    if (roots_im.get(k).isZero()) {
                        DoubleDouble alpha = roots_re.get(k).neg();
                        coefReal[k + 1] = coefReal[k];
                        for (int j = k; j >= 1; j--) {
                            coefReal[j] = coefReal[j - 1].add(alpha.mul(coefReal[j]));
                        }
                        coefReal[0] = coefReal[0].mul(alpha);
                    } else {
                        DoubleDouble alpha = roots_re.get(k).sqr().add(roots_im.get(k).sqr());
                        DoubleDouble beta = roots_re.get(k).mulPowerOf2(2).neg();

                        coefReal[k + 2] = coefReal[k];
                        if (k > 0) {
                            coefReal[k + 1] = coefReal[k - 1].add(coefReal[k].mul(beta));
                            for (int j = k; j >= 2; j--) {
                                coefReal[j] = coefReal[j - 2].add(coefReal[j - 1].mul(beta)).add(coefReal[j].mul(alpha));
                            }
                            coefReal[1] = coefReal[0].mul(beta).add(coefReal[1].mul(alpha));
                        } else {
                            coefReal[1] = coefReal[0].mul(beta);
                        }
                        coefReal[0] = coefReal[0].mul(alpha);

                        // Skip conjugate root
                        k++;
                    }
                }
            } else {
                coefComplex[0] = DoubleComplex.ONE;

                for (int k = 0; k < roots_re.size(); k++) {
                    DoubleDouble alpha_re = roots_re.get(k).neg();
                    DoubleDouble alpha_im = roots_im.get(k).neg();
                    DoubleComplex alpha = new DoubleComplex(alpha_re, alpha_im);
                    coefComplex[k + 1] = coefComplex[k];
                    for (int j = k; j >= 1; j--) {
                        coefComplex[j] = coefComplex[j - 1].add(alpha.mul(coefComplex[j]));
                    }
                    coefComplex[0] = alpha.mul(coefComplex[0]);
                }
            }
        }
    }

    /**
     * Constructor, mainly used for testing purposes. This constructor can be
     * used to create a real polynomial with a set of real roots. Create
     * coefficients of polynomial, based on supplied roots. This constructor
     * determines the real coefficients of
     *
     * factor * (z-roots[0]) * (z-roots[1]) * . . . * (z-roots[n-1]).
     *
     @param factor The constant scaling factor, mentioned above.
     @param roots The set of complex roots. The degree of the polynomial is
     * equal to the length of this array.
     */
    public void run00(DoubleDouble factor, DoubleDouble[] roots) {
        coefReal = new DoubleDouble[roots.length + 1];
        coefComplex = null;

        coefReal[0] = factor;
        for (int r = 0; r < roots.length; r++) {
            DoubleDouble p = roots[r].neg();
            coefReal[r + 1] = coefReal[r];
            for (int s = r; s > 0; s--) {
                coefReal[s] = coefReal[s - 1].add(p.mul(coefReal[s]));    // coefReal[s-1] + p*coefReal[s];
            }
            coefReal[0] = coefReal[0].mul(p);
        }
    }

    /**
     * Constructor, mainly used for testing purposes. This constructor can be
     * used to create a real polynomial with a set of real roots. Create
     * coefficients of polynomial, based on supplied roots. This constructor
     * determines the real coefficients of
     *
     * factor * (z-roots[0]) * (z-roots[1]) * . . . * (z-roots[n-1]).
     *
     @param factor The constant scaling factor, mentioned above.
     @param roots The set of complex roots. The degree of the polynomial is
     * equal to the length of this array.
     */
    public void run2(double factor, double[] roots) {
        coefReal = new DoubleDouble[roots.length + 1];
        coefComplex = null;

        coefReal[0] = new DoubleDouble(factor);
        for (int r = 0; r < roots.length; r++) {
            DoubleDouble p = new DoubleDouble(roots[r]).neg();
            coefReal[r + 1] = coefReal[r];
            for (int s = r; s > 0; s--) {
                coefReal[s] = coefReal[s - 1].add(p.mul(coefReal[s]));    // coefReal[s-1] + p*coefReal[s];
            }
            coefReal[0] = coefReal[0].mul(p);
        }
    }

    /**
     * Constructor, mainly used for testing purposes. This constructor can be
     * used to create a real polynomial with a set of real roots. Create
     * coefficients of polynomial, based on supplied roots. This constructor
     * determines the real coefficients of
     *
     * factor * (z-roots[0]) * (z-roots[1]) * . . . * (z-roots[n-1]).
     *
     @param factor The constant scaling factor, mentioned above.
     @param roots The set of complex roots. The degree of the polynomial is
     * equal to the length of this array.
     */
    public void run01(String factor, String[] roots) {
        coefReal = new DoubleDouble[roots.length + 1];
        coefComplex = null;

        coefReal[0] = new DoubleDouble(factor);
        for (int r = 0; r < roots.length; r++) {
            DoubleDouble p = new DoubleDouble(roots[r]).neg();
            coefReal[r + 1] = coefReal[r];
            for (int s = r; s > 0; s--) {
                coefReal[s] = coefReal[s - 1].add(p.mul(coefReal[s]));    // coefReal[s-1] + p*coefReal[s];
            }
            coefReal[0] = coefReal[0].mul(p);
        }
    }

    /**
     * Constructor, mainly used for testing purposes. This constructor can be
     * used to create a polynomial with a set of non-conjugate complex roots.
     * Create coefficients of polynomial, based on supplied roots. This
     * constructor determines the complex coefficients of
     *
     * factor * (z-roots[0]) * (z-roots[1]) * . . . * (z-roots[n-1]).
     *
     @param factor The constant scaling factor, mentioned above.
     @param roots The set of complex roots. The degree of the polynomial is
     * equal to the length of this array.
     */
    public void run02(DoubleComplex factor, DoubleComplex[] roots) {
        coefComplex = new DoubleComplex[roots.length + 1];
        coefReal = null;

        coefComplex[0] = factor;
        for (int r = 0; r < roots.length; r++) {
            DoubleComplex p = roots[r].neg();
            coefComplex[r + 1] = coefComplex[r];
            for (int s = r; s > 0; s--) {
                coefComplex[s] = coefComplex[s - 1].add(p.mul(coefComplex[s]));    // coefComplex[s-1] + p*coefComplex[s];
            }
            coefComplex[0] = coefComplex[0].mul(p);
        }
    }

    /**
     * Constructor, mainly used for testing purposes. This constructor can be
     * used to create a polynomial with a set of non-conjugate complex roots.
     * Create coefficients of polynomial, based on supplied roots. This
     * constructor determines the complex coefficients of
     *
     * factor * (z-roots[0]) * (z-roots[1]) * . . . * (z-roots[n-1]).
     *
     @param factor The constant scaling factor, mentioned above.
     @param roots The set of complex roots. The degree of the polynomial is
     * equal to the length of this array.
     */
    public void run03(Complex factor, Complex[] roots) {
        coefComplex = new DoubleComplex[roots.length + 1];
        coefReal = null;

        coefComplex[0] = new DoubleComplex(factor);
        for (int r = 0; r < roots.length; r++) {
            DoubleComplex p = new DoubleComplex(roots[r].neg());
            coefComplex[r + 1] = coefComplex[r];
            for (int s = r; s > 0; s--) {
                coefComplex[s] = coefComplex[s - 1].add(p.mul(coefComplex[s]));    // coefComplex[s-1] + p*coefComplex[s];
            }
            coefComplex[0] = coefComplex[0].mul(p);
        }
    }

    /**
     * Constructor, mainly used for testing purposes. This constructor can be
     * used to create a polynomial with a set of real roots and a set of
     * conjugate complex roots.
     *
     @param factor_re The real part of the constant scaling factor, which is
     * used to multiply the monic polynomial, derived from the given roots (see
     * below).
     @param factor_im The imaginary part of the constant scaling factor, which
     * is used to multiply the monic polynomial, derived from the given roots.
     @param roots_re The real part of the set of complex roots.
     @param roots_im The imaginary part of the set of complex roots. The
     * parameter roots_re and roots_im must have the same length. The degree of
     * the derived polynomial is equal to the length of these arrays.
     */
    public void run04(String factor_re, String factor_im, String[] roots_re, String[] roots_im) {
        coefComplex = new DoubleComplex[roots_re.length + 1];
        coefReal = null;

        coefComplex[0] = new DoubleComplex(factor_re, factor_im);
        for (int r = 0; r < roots_re.length; r++) {
            DoubleComplex p = new DoubleComplex(roots_re[r], roots_im[r]).neg();
            coefComplex[r + 1] = coefComplex[r];
            for (int s = r; s > 0; s--) {
                coefComplex[s] = coefComplex[s - 1].add(p.mul(coefComplex[s]));    // coefComplex[s-1] + p*coefComplex[s];
            }
            coefComplex[0] = coefComplex[0].mul(p);
        }
    }

    /**
     * Constructor, mainly used for testing purposes. This constructor can be
     * used to create a polynomial with a set of real roots and a set of
     * conjugate complex roots.
     *
     @param factor_re The real part of the constant scaling factor, which is
     * used to multiply the monic polynomial, derived from the given roots (see
     * below).
     @param factor_im The imaginary part of the constant scaling factor, which
     * is used to multiply the monic polynomial, derived from the given roots.
     @param roots_re The real part of the set of complex roots.
     @param roots_im The imaginary part of the set of complex roots. The
     * parameter roots_re and roots_im must have the same length. The degree of
     * the derived polynomial is equal to the length of these arrays.
     */
    public void run3(double factor_re, double factor_im, double[] roots_re, double[] roots_im) {
        coefComplex = new DoubleComplex[roots_re.length + 1];
        coefReal = null;

        coefComplex[0] = new DoubleComplex(factor_re, factor_im);
        for (int r = 0; r < roots_re.length; r++) {
            DoubleComplex p = new DoubleComplex(roots_re[r], roots_im[r]).neg();
            coefComplex[r + 1] = coefComplex[r];
            for (int s = r; s > 0; s--) {
                coefComplex[s] = coefComplex[s - 1].add(p.mul(coefComplex[s]));    // coefComplex[s-1] + p*coefComplex[s];
            }
            coefComplex[0] = coefComplex[0].mul(p);
        }
    }

    /**
     * A wrapper method which calls the PZeros-solver to obtain the roots of the
     * polynomial. This tester derives the coefficients of the polynomial to be
     * solved from the double double precision coefficients, stored in the
     * tester. It uses double double precision arithmatic as well for computing
     * roots.
     *
     @param roots A pre-allocated array, in which the roots are stored. This
     * array must have a length of at least the degree of the polynomial.
     @param rad Array, which gives an indication of the accuracy of the found
     * roots. For each root, a radius is returned. The root is assured to be in
     * the disk with the corresponding radius, centered around the returned
     * root.
     @param err Array, which specifies whether the root and its corresponding
     * radius of accuracy could be determined correctly. If err[i] is true, then
     * the program did not converge for root[i].
     @return Returns the degree of the polynomial if the computation succeeds,
     * and returns a value less than the degree of the polynomial if an error
     * occurs (e.g. convergence failure). When a value less than the degree of
     * the polynomial is returned, then only part (or none) of the roots could
     * be determined. If a negative value is returned, then the supplied input
     * is not correct: -1: Leading coefficient equals 0. -2: Coefficient for x^0
     * (constant coefficient) equals 0. -3: Ratio of smallest coefficient
     * magnitude and largest coefficient magnitude is too large and will lead to
     * underflow/overflow.
     */
    public int solve(DoubleComplex[] roots, double[] rad, boolean[] err) {
        PZeros pol = (coefReal != null) ? new PZeros(coefReal) : new PZeros(coefComplex);
        return pol.solve(roots, rad, err);
    }

    /**
     * A wrapper method which calls the PZeros-solver to obtain the roots of the
     * polynomial. This tester derives the coefficients of the polynomial to be
     * solved from the double double precision coefficients, stored in the
     * tester. It uses double double precision arithmatic as well for computing
     * roots.
     *
     @param roots A pre-allocated array, in which the roots are stored. This
     * array must have a length of at least the degree of the polynomial.
     @return Returns the degree of the polynomial if the computation succeeds,
     * and returns a value less than the degree of the polynomial if an error
     * occurs (e.g. convergence failure). When a value less than the degree of
     * the polynomial is returned, then only part (or none) of the roots could
     * be determined. If a negative value is returned, then the supplied input
     * is not correct: -1: Leading coefficient equals 0. -2: Coefficient for x^0
     * (constant coefficient) equals 0. -3: Ratio of smallest coefficient
     * magnitude and largest coefficient magnitude is too large and will lead to
     * underflow/overflow.
     */
    public int solve(DoubleComplex[] roots) {
        PZeros pol = (coefReal != null) ? new PZeros(coefReal) : new PZeros(coefComplex);
        boolean[] err = new boolean[degree()];
        double[] rad = new double[degree()];
        return pol.solve(roots, rad, err);
    }

    /**
     * A wrapper method which calls the PZeros-solver to obtain the roots of the
     * polynomial. This tester derives the coefficients of the polynomial to be
     * solved from the double double precision coefficients, stored in the
     * tester. It uses normal double precision arithmatic for computing roots.
     *
     @param roots A pre-allocated array, in which the roots are stored. This
     * array must have a length of at least the degree of the polynomial.
     @param rad Array, which gives an indication of the accuracy of the found
     * roots. For each root, a radius is returned. The root is assured to be in
     * the disk with the corresponding radius, centered around the returned
     * root.
     @param err Array, which specifies whether the root and its corresponding
     * radius of accuracy could be determined correctly. If err[i] is true, then
     * the program did not converge for root[i].
     @return Returns the degree of the polynomial if the computation succeeds,
     * and returns a value less than the degree of the polynomial if an error
     * occurs (e.g. convergence failure). When a value less than the degree of
     * the polynomial is returned, then only part (or none) of the roots could
     * be determined. If a negative value is returned, then the supplied input
     * is not correct: -1: Leading coefficient equals 0. -2: Coefficient for x^0
     * (constant coefficient) equals 0. -3: Ratio of smallest coefficient
     * magnitude and largest coefficient magnitude is too large and will lead to
     * underflow/overflow.
     */
    public int solve(Complex[] roots, double[] rad, boolean[] err) {
        PZeros pol = (coefReal != null) ? new PZeros(DoubleDouble.toDouble(coefReal)) : new PZeros(DoubleComplex.toComplex(coefComplex));
        return pol.solve(roots, rad, err);
    }

    /**
     * A wrapper method which calls the PZeros-solver to obtain the roots of the
     * polynomial. This tester derives the coefficients of the polynomial to be
     * solved from the double double precision coefficients, stored in the
     * tester. It uses normal double precision arithmatic for computing roots.
     *
     @param roots A pre-allocated array, in which the roots are stored. This
     * array must have a length of at least the degree of the polynomial.
     @return Returns the degree of the polynomial if the computation succeeds,
     * and returns a value less than the degree of the polynomial if an error
     * occurs (e.g. convergence failure). When a value less than the degree of
     * the polynomial is returned, then only part (or none) of the roots could
     * be determined. If a negative value is returned, then the supplied input
     * is not correct: -1: Leading coefficient equals 0. -2: Coefficient for x^0
     * (constant coefficient) equals 0. -3: Ratio of smallest coefficient
     * magnitude and largest coefficient magnitude is too large and will lead to
     * underflow/overflow.
     */
    public int solve(Complex[] roots) {
        PZeros pol = (coefReal != null) ? new PZeros(DoubleDouble.toDouble(coefReal)) : new PZeros(DoubleComplex.toComplex(coefComplex));
        boolean[] err = new boolean[degree()];
        double[] rad = new double[degree()];
        return pol.solve(roots, rad, err);
    }

    /**
     * Returns the degree of the polynomial.
     *
     @return The degree of the polynomial.
     */
    int degree() {
        return -1 + ((coefReal != null) ? coefReal.length : coefComplex.length);
    }

    /**
     * Creates a new PZerosTester, with the real part of the coefficients
     * retained only. The imaginary part is set to 0.
     *
     @return PZerosTester with the imaginary part of the coefficients set to
     * 0.
     */
    public void keepRealPart() {
        if (coefComplex == null) {
            return;
        }

        DoubleDouble[] coef = new DoubleDouble[coefComplex.length];
        for (int i = 0; i < coef.length; i++) {
            coef[i] = coefComplex[i].re;
        }
        run0(coef);
    }

    public static void main(String[] args) {
        //mainT();   // Polynomials from TestPolynomials.java
        //mainWilco();
        //mainD();
        //mainDCircle();  // 53-bits precision random polynomials
        //    mainDD(); // 105-bit precision random polynomials
        //mainDDCircle();
    }

    public void testT() {

        List<String> names = new ArrayList<>();
        List<String[]> pols = new ArrayList<>();
        PolynomialsUtil.getPolynomials(names, pols);
        for (int i = 0; i < names.size(); i++) {
            run1(pols.get(i));
            DoubleComplex[] roots = new DoubleComplex[degree()];
            int degree = solve(roots);
            if (degree < degree()) {
                System.err.println("CONVERGENCE ERROR FOR : " + names.get(i));
                continue;
            }

            System.out.println(names.get(i));
            System.out.println("---------------");
            for (int j = 0; j < degree; j++) {
                System.out.println(roots[j]);
            }
            System.out.println();
        }

        {
            // Double double 105-bit precision
            run2(1, new double[]{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25});

            DoubleComplex[] root = new DoubleComplex[degree()];
            int result = solve(root);
            Arrays.asList(root).forEach(r -> System.out.println("" + r.real() + "   " + r.imag()));
        }
        System.out.println();
        {
            // Normal 53-bit precision
            run3(1, 1, new double[]{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22},
                    new double[]{0, 2, 1, 4, 4, 6, 7, 1, 9, -3, -5, 12, 13, 14, -5, 16, -7, 18, -9, 20, -1, 22});
            Complex[] root = new Complex[degree()];
            int result = solve(root);
            Arrays.asList(root).forEach(r -> System.out.println("" + r.real() + "   " + r.imag()));
        }
    }

    public void testWilco() {
        double[] rre = {-3.02, -3.86, -1.94, -0.41, -2.27, 4.97, -2.39, -2.37, -1.68, -2.17, -2.31, -3.13, -3.34, -3.03, 0.3};
        double[] rim = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.28, 1.12, -3.93, -4.87};

        run3(9.643264256234523, 0.856475436745643, rre, rim);

        Complex[] roots = new Complex[degree()];
        double[] radius = new double[degree()];
        boolean[] err = new boolean[degree()];
        int result = solve(roots, radius, err);
        Arrays.asList(roots).forEach(System.out::println);
    }

    public void testD() {
        int nErrors = 0;
        int nInaccs = 0;
        Random rnd = new Random(System.currentTimeMillis());

        int nMaxIter = 5000;//500000;
        for (int i = 0; i < nMaxIter; i++) {

            int degr = (int) (4 + 8 * rnd.nextDouble());
            int degc = (int) (4 + 8 * rnd.nextDouble());
            if (degc % 2 == 1) {
                degc++;
            }
            double[] r_re = new double[degr + degc];
            double[] r_im = new double[degr + degc];
            for (int k = 0; k < degr; k++) {
                r_re[k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_re[k] == 0.0) {
                    r_re[k] = 0.005;
                }
                r_im[k] = 0.0;
            }
            for (int k = 0; k < degc / 2; k++) {
                r_re[degr + 2 * k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_re[degr + 2 * k] == 0.0) {
                    r_re[degr + 2 * k] = 0.005;
                }
                r_im[degr + 2 * k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_im[degr + 2 * k] == 0.0) {
                    r_im[degr + 2 * k] = 0.005;
                }
                r_re[degr + 2 * k + 1] = r_re[degr + 2 * k];
                r_im[degr + 2 * k + 1] = -r_im[degr + 2 * k];
            }
            run3(0.1 + 10 * rnd.nextDouble(), 0.1 + 10 * rnd.nextDouble(), r_re, r_im);
            keepRealPart();   // Can be used for making a real polynomial, because the roots
            // are either all real, or are conjugates.

            Complex[] roots = new Complex[degree()];
            double[] radius = new double[degree()];
            boolean[] err = new boolean[degree()];
            int result = solve(roots, radius, err);

            if (result != roots.length) {
                System.err.println("CONVERGENCE ERROR DEGREE = " + degree());
                for (int k = 0; k < degree(); k++) {
                    System.err.println("   " + r_re[k] + " + i*" + r_im[k]);
                }
                nErrors++;
                continue;
            }

            boolean hasError = false;
            boolean hasInacc = false;
            for (int k = 0; k < degree(); k++) {
                Info info = getClosestDistanceFromRoot(r_re[k], r_im[k], roots);
                if (err[info.idx]) {
                    System.err.println("ERROR FOR ROOT FOR REAL POLY " + i + ": ROOT = " + r_re[k] + " + i*" + r_im[k] + "  DIFF = " + info.diff + "   RAD = " + radius[info.idx]);
                    hasError = true;
                }
                if (info.diff > 1e-3) {
                    System.err.println("INACCURATE ROOT FOR REAL POLY " + i + ": ROOT = " + r_re[k] + " + i*" + r_im[k] + "  DIFF = " + info.diff + "   RAD = " + radius[info.idx]);
                    hasInacc = true;
                }
            }
            if (hasError || hasInacc) {
                for (int kk = 0; kk < degree(); kk++) {
                    System.err.println("    root = " + r_re[kk] + " + i*" + r_im[kk]);
                }
                if (hasError) {
                    nErrors++;
                }
                if (hasInacc) {
                    nInaccs++;
                }
            }
        }

        // Complex polynomials
        for (int i = 0; i < nMaxIter; i++) {

            int degr = (int) (4 + 8 * rnd.nextDouble());
            int degc = (int) (4 + 8 * rnd.nextDouble());
            double[] r_re = new double[degr + degc];
            double[] r_im = new double[degr + degc];
            for (int k = 0; k < degr; k++) {
                r_re[k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_re[k] == 0.0) {
                    r_re[k] = 0.005;
                }
                r_im[k] = 0.0;
            }
            for (int k = 0; k < degc; k++) {
                r_re[degr + k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_re[degr + k] == 0.0) {
                    r_re[degr + k] = 0.005;
                }
                r_im[degr + k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_im[degr + k] == 0.0) {
                    r_im[degr + k] = 0.005;
                }
            }
            run3(0.1 + 10 * rnd.nextDouble(), 0.1 + 10 * rnd.nextDouble(), r_re, r_im);

            Complex[] roots = new Complex[degree()];
            double[] radius = new double[degree()];
            boolean[] err = new boolean[degree()];
            int result = solve(roots, radius, err);

            if (result != roots.length) {
                System.err.println("CONVERGENCE ERROR DEGREE = " + degree());
                for (int k = 0; k < degree(); k++) {
                    System.err.println("   " + r_re[k] + " + i*" + r_im[k]);
                }
                nErrors++;
                continue;
            }

            boolean hasInacc = false;
            for (int k = 0; k < degr + degc; k++) {
                Info info = getClosestDistanceFromRoot(r_re[k], r_im[k], roots);
                if (info.diff > 1e-3) {
                    System.err.println("INACCURATE ROOT FOR CPLX POLY " + i + ": ROOT = " + r_re[k] + " + i*" + r_im[k] + "  DIFF = " + info.diff + "   RAD = " + radius[info.idx]);
                    hasInacc = true;
                }
            }
            if (hasInacc) {
                for (int kk = 0; kk < degree(); kk++) {
                    System.err.println("    root = " + r_re[kk] + " + i*" + r_im[kk]);
                }
                nInaccs++;
            }
        }
        System.out.println("ERRORS = " + nErrors);
        System.out.println("INACCURACIES = " + nInaccs);
    }

    public void testDCircle() {
        int nErrors = 0;
        int nInaccs = 0;
        Random rnd = new Random(System.currentTimeMillis());
        int nMaxIter = 5000;//500000
        for (int i = 0; i < nMaxIter; i++) {

            int degc = (int) (8 + 15 * rnd.nextDouble());
            if (degc % 2 == 1) {
                degc++;
            }
            double[] r_re = new double[degc];
            double[] r_im = new double[degc];
            for (int k = 0; k < degc / 2; k++) {
                double phi = Math.PI * rnd.nextDouble();
                double r = 0.9 + 0.2 * rnd.nextDouble();
                r_re[2 * k] = r * Math.cos(phi);
                r_im[2 * k] = r * Math.sin(phi);
                r_re[2 * k + 1] = r_re[2 * k];
                r_im[2 * k + 1] = -r_im[2 * k];
            }
            run3(0.1 + 10 * rnd.nextDouble(), 0.1 + 10 * rnd.nextDouble(), r_re, r_im);
            keepRealPart();   // Can be used for making a real polynomial, because the roots
            // are either all real, or are conjugates.

            Complex[] roots = new Complex[degree()];
            double[] radius = new double[degree()];
            boolean[] err = new boolean[degree()];
            int result = solve(roots, radius, err);

            if (result != roots.length) {
                System.err.println("CONVERGENCE ERROR DEGREE = " + degree());
                for (int k = 0; k < degree(); k++) {
                    System.err.println("   " + r_re[k] + " + i*" + r_im[k]);
                }
                nErrors++;
                continue;
            }

            boolean hasError = false;
            boolean hasInacc = false;
            for (int k = 0; k < degree(); k++) {
                Info info = getClosestDistanceFromRoot(r_re[k], r_im[k], roots);
                if (err[info.idx]) {
                    System.err.println("ERROR FOR ROOT FOR REAL POLY " + i + ": ROOT = " + r_re[k] + " + i*" + r_im[k] + "  DIFF = " + info.diff + "   RAD = " + radius[info.idx]);
                    hasError = true;
                }
                if (info.diff > 1e-3) {
                    System.err.println("INACCURATE ROOT FOR REAL POLY " + i + ": ROOT = " + r_re[k] + " + i*" + r_im[k] + "  DIFF = " + info.diff + "   RAD = " + radius[info.idx]);
                    hasInacc = true;
                }
            }
            if (hasError || hasInacc) {
                for (int kk = 0; kk < degree(); kk++) {
                    System.err.println("    root = " + r_re[kk] + " + i*" + r_im[kk]);
                }
                if (hasError) {
                    nErrors++;
                }
                if (hasInacc) {
                    nInaccs++;
                }
            }
        }
        System.out.println("ERRORS = " + nErrors);
        System.out.println("INACCURACIES = " + nInaccs);
    }

    public void testDD() {
        int nErrors = 0;
        int nInaccs = 0;
        Random rnd = new Random(System.currentTimeMillis());
        int nMaxIter = 5000;
        for (int i = 0; i < nMaxIter; i++) {
            int degr = (int) (6 + 20 * rnd.nextDouble());
            int degc = (int) (6 + 20 * rnd.nextDouble());
            if (degc % 2 == 1) {
                degc++;
            }
            double[] r_re = new double[degr + degc];
            double[] r_im = new double[degr + degc];
            for (int k = 0; k < degr; k++) {
                r_re[k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_re[k] == 0.0) {
                    r_re[k] = 0.005;
                }
                r_im[k] = 0.0;
            }
            for (int k = 0; k < degc / 2; k++) {
                r_re[degr + 2 * k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_re[degr + 2 * k] == 0.0) {
                    r_re[degr + 2 * k] = 0.005;
                }
                r_im[degr + 2 * k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_im[degr + 2 * k] == 0.0) {
                    r_im[degr + 2 * k] = 0.005;
                }
                r_re[degr + 2 * k + 1] = r_re[degr + 2 * k];
                r_im[degr + 2 * k + 1] = -r_im[degr + 2 * k];
            }
            run3(0.1 + 10 * rnd.nextDouble(), 0.1 + 10 * rnd.nextDouble(), r_re, r_im);
            keepRealPart();   // Can be used for making a real polynomial, because the roots
            // are either all real, or are conjugates.

            DoubleComplex[] roots = new DoubleComplex[degree()];
            double[] radius = new double[degree()];
            boolean[] err = new boolean[degree()];
            int result = solve(roots, radius, err);

            if (result != roots.length) {
                System.err.println("CONVERGENCE ERROR DEGREE = " + degree());
                for (int k = 0; k < degr + degc; k++) {
                    System.err.println("   " + r_re[k] + " i*" + r_im[k]);
                }
                nErrors++;
                continue;
            }

            boolean hasError = false;
            boolean hasInacc = false;
            for (int k = 0; k < degree(); k++) {
                Info info = PZerosTest.getClosestDistanceFromRoot(r_re[k], r_im[k], roots);
                if (err[info.idx]) {
                    System.err.println("ERROR FOR ROOT FOR REAL POLY " + i + ": ROOT = " + roots[info.idx].toString(10, 8, false) + "  DIFF = " + info.diff + "   RAD = " + radius[info.idx]);
                    hasError = true;
                }
                if (info.diff > 1e-6) {
                    System.err.println("INACCURATE ROOT FOR REAL POLY " + i + ": ROOT = " + roots[info.idx].toString(10, 8, false) + "  DIFF = " + info.diff + "   RAD = " + radius[info.idx]);
                    hasInacc = true;
                }
            }
            if (hasError || hasInacc) {
                for (int kk = 0; kk < degree(); kk++) {
                    System.err.println("    root = " + new DoubleComplex(r_re[kk], r_im[kk]).toString(10, 4, false));
                }
                if (hasError) {
                    nErrors++;
                }
                if (hasInacc) {
                    nInaccs++;
                }
            }
        }

        for (int i = 0; i < nMaxIter; i++) {

            int degr = (int) (6 + 20 * rnd.nextDouble());
            int degc = (int) (6 + 20 * rnd.nextDouble());
            double[] r_re = new double[degr + degc];
            double[] r_im = new double[degr + degc];
            for (int k = 0; k < degr; k++) {
                r_re[k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_re[k] == 0.0) {
                    r_re[k] = 0.005;
                }
                r_im[k] = 0.0;
            }
            for (int k = 0; k < degc; k++) {
                r_re[degr + k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_re[degr + k] == 0.0) {
                    r_re[degr + k] = 0.005;
                }
                r_im[degr + k] = ((int) (-500 + 1000 * rnd.nextDouble())) * 0.01;
                if (r_im[degr + k] == 0.0) {
                    r_im[degr + k] = 0.005;
                }
            }
            run3(0.1 + 10 * rnd.nextDouble(), 0.1 + 10 * rnd.nextDouble(), r_re, r_im);

            DoubleComplex[] roots = new DoubleComplex[degree()];
            double[] radius = new double[degree()];
            boolean[] err = new boolean[degree()];
            int result = solve(roots, radius, err);

            if (result != roots.length) {
                System.err.println("CONVERGENCE ERROR DEGREE = " + degree());
                for (int k = 0; k <= degree(); k++) {
                    System.err.println("   " + r_re[k] + " i*" + r_im[k]);
                }
                nErrors++;
                continue;
            }

            boolean hasError = false;
            boolean hasInacc = false;
            for (int k = 0; k < degr + degc; k++) {
                Info info = PZerosTest.getClosestDistanceFromRoot(r_re[k], r_im[k], roots);
                if (err[info.idx]) {
                    System.err.println("ERROR FOR ROOT FOR CPLX POLY " + i + ": ROOT = " + roots[info.idx].toString(10, 8, false) + "  DIFF = " + info.diff + "   RAD = " + radius[info.idx]);
                    hasError = true;
                }
                if (info.diff > 1e-6) {
                    System.err.println("INACCURATE ROOT FOR CPLX POLY " + i + ": ROOT = " + roots[info.idx].toString(10, 8, false) + "  DIFF = " + info.diff + "   RAD = " + radius[info.idx]);
                    hasInacc = true;
                }
            }
            if (hasError || hasInacc) {
                for (int kk = 0; kk < degree(); kk++) {
                    System.err.println("    root = " + new DoubleComplex(r_re[kk], r_im[kk]).toString(10, 4, false));
                }
                if (hasError) {
                    nErrors++;
                }
                if (hasInacc) {
                    nInaccs++;
                }
            }
        }
        System.out.println("ERRORS = " + nErrors);
        System.out.println("INACCURACIES = " + nInaccs);
    }

    public void testDDCircle() {
        int nErrors = 0;
        int nInaccs = 0;
        Random rnd = new Random(System.currentTimeMillis());
        int nMaxIter = 5000;//50000
        for (int i = 0; i < nMaxIter; i++) {
            int degc = (int) (12 + 40 * rnd.nextDouble());
            if (degc % 2 == 1) {
                degc++;
            }
            double[] r_re = new double[degc];
            double[] r_im = new double[degc];
            for (int k = 0; k < degc / 2; k++) {
                double phi = Math.PI * rnd.nextDouble();
                double r = 0.9 + 0.2 * rnd.nextDouble();
                r_re[2 * k] = r * Math.cos(phi);
                r_im[2 * k] = r * Math.sin(phi);
                r_re[1 + 2 * k] = r * Math.cos(phi);
                r_im[1 + 2 * k] = -r * Math.sin(phi);
            }
            run3(0.1 + 10 * rnd.nextDouble(), 0.1 + 10 * rnd.nextDouble(), r_re, r_im);
            keepRealPart();   // Can be used for making a real polynomial, because the roots
            // are either all real, or are conjugates.

            DoubleComplex[] roots = new DoubleComplex[degree()];
            double[] radius = new double[degree()];
            boolean[] err = new boolean[degree()];
            int result = solve(roots, radius, err);

            if (result != roots.length) {
                System.err.println("CONVERGENCE ERROR DEGREE = " + degree());
                for (int k = 0; k < degc; k++) {
                    System.err.println("   " + r_re[k] + " i*" + r_im[k]);
                }
                nErrors++;
                continue;
            }

            boolean hasError = false;
            boolean hasInacc = false;
            for (int k = 0; k < degree(); k++) {
                Info info = PZerosTest.getClosestDistanceFromRoot(r_re[k], r_im[k], roots);
                if (err[info.idx]) {
                    System.err.println("ERROR FOR ROOT FOR CIRCLE POLY " + i + ": ROOT = " + roots[info.idx].toString(10, 8, false) + "  DIFF = " + info.diff + "   RAD = " + radius[info.idx]);
                    hasError = true;
                }
                if (info.diff > 1e-6) {
                    System.err.println("INACCURATE ROOT FOR CIRCLE POLY " + i + ": ROOT = " + roots[info.idx].toString(10, 8, false) + "  DIFF = " + info.diff + "   RAD = " + radius[info.idx]);
                    hasInacc = true;
                }
            }
            if (hasError || hasInacc) {
                for (int kk = 0; kk < degree(); kk++) {
                    System.err.println("    root = " + new DoubleComplex(r_re[kk], r_im[kk]).toString(10, 4, false));
                }
                if (hasError) {
                    nErrors++;
                }
                if (hasInacc) {
                    nInaccs++;
                }
            }
        }
        System.out.println("ERRORS = " + nErrors);
        System.out.println("INACCURACIES = " + nInaccs);
    }

    static Info getClosestDistanceFromRoot(double re, double im, DoubleComplex[] r) {
        double abs = Math.abs(re) + Math.abs(im);
        double minDiff = 1e10;
        int idx = -1;
        for (int i = 0; i < r.length; i++) {
            double diff = (Math.abs(re - r[i].real().doubleValue()) + Math.abs(im - r[i].imag().doubleValue())) / abs;
            if (diff < minDiff) {
                minDiff = diff;
                idx = i;
            }
        }
        return new Info(minDiff, idx);
    }

    static Info getClosestDistanceFromRoot(double re, double im, Complex[] r) {
        double abs = Math.abs(re) + Math.abs(im);
        double minDiff = 1e10;
        int idx = -1;
        for (int i = 0; i < r.length; i++) {
            double diff = (Math.abs(re - r[i].real()) + Math.abs(im - r[i].imag())) / abs;
            if (diff < minDiff) {
                minDiff = diff;
                idx = i;
            }
        }
        return new Info(minDiff, idx);
    }

    public static class Info {

        double diff;
        int idx;

        Info(double diff, int idx) {
            this.diff = diff;
            this.idx = idx;
        }
    }
}
