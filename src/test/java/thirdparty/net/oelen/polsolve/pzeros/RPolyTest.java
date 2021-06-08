package thirdparty.net.oelen.polsolve.pzeros;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import junit.framework.TestCase;
import thirdparty.net.oelen.polarith.DoubleDouble;
import thirdparty.net.oelen.polsolve.jt.RPoly;


// This is a class for testing purposes of RPoly. The polynomial root
// finder can be tested by supplying a polynomial with known roots.
// This tester class allows creating polynomials, based on a given set
// of known roots. The tester multiplies the factors of the form
// (x - ri), with ri the i-th root to obtain a monic polynomial of
// degree equal to the number of roots. A separate scaling factor can
// be supplied to multiply all coefficients, so that non-monic 
// polynomials also can be tested.
// The computation of the coefficients is done in double double
// precision to assure that inaccuracies of the tester class do
// not add false alarms about the polynomial root finder.

public class RPolyTest extends TestCase {

    private DoubleDouble coef[];
    
    /**
     * Constructor, used for testing purposes.
     * 
     * @param testPolynomial A String array, as specified in the TestPolynomials
     * class. This is converted to a polynomial, useful for solving/testing.
     */
    public void init0(String[] testPolynomial) {
        List<DoubleDouble> roots_re = new ArrayList<>();
        List<DoubleDouble> roots_im = new ArrayList<>();
        if (testPolynomial[0].equals("RR")) {
            // Roots Real
            for (int i=1; i<testPolynomial.length; i++) {
                String[] num = testPolynomial[i].split(";");
                roots_re.add(new DoubleDouble(num[0]));
                if (num.length >= 2) {
                    roots_re.add(new DoubleDouble(num[0]));
                    roots_im.add(new DoubleDouble(num[1]));
                    roots_im.add(new DoubleDouble(num[1]).neg());
                }
                else {
                    roots_im.add(DoubleDouble.ZERO);
                }
            }
            int degree = roots_re.size();
            coef = new DoubleDouble[degree+1];
        }
        else if (testPolynomial[0].equals("RC")) {
            // Roots Complex
            coef = null;
            throw new RuntimeException("RPolyTester does not support non-conjugate complex roots.");
        }
        else if (testPolynomial[0].equals("CR")) {
            // Coefs Real
            int ncoefs = testPolynomial.length - 1;
            coef = new DoubleDouble[ncoefs];
            for (int i=1; i<testPolynomial.length; i++) {
                coef[i-1] = new DoubleDouble(testPolynomial[i]);
            }
        }
        else if (testPolynomial[0].equals("CC")) {
            // Coefs Complex
            coef = null;
            throw new RuntimeException("RPolyTester does not support complex coefficients.");
        }
        else {
            coef = null;
        }
        
        if (!roots_re.isEmpty()) {
            coef[0] = DoubleDouble.ONE;

            for (int k=0; k<roots_re.size(); k++) {
                if (roots_im.get(k).isZero()) {
                    DoubleDouble alpha = roots_re.get(k).neg();
                    coef[k + 1] = coef[k];
                    for (int j = k; j >= 1; j--) {
                        coef[j] = coef[j - 1].add(alpha.mul(coef[j]));
                    }
                    coef[0] = coef[0].mul(alpha);
                }
                else {
                    DoubleDouble alpha = roots_re.get(k).sqr().add(roots_im.get(k).sqr());
                    DoubleDouble beta = roots_re.get(k).mulPowerOf2(2).neg();

                    coef[k + 2] = coef[k];
                    if (k > 0) {
                        coef[k+1] = coef[k-1].add(coef[k].mul(beta));
                        for (int j = k; j >= 2; j--) {
                            coef[j] = coef[j - 2].add(coef[j - 1].mul(beta)).add(coef[j].mul(alpha));
                        }
                        coef[1] = coef[0].mul(beta).add(coef[1].mul(alpha));
                    } else {
                        coef[1] = coef[0].mul(beta);
                    }
                    coef[0] = coef[0].mul(alpha);
                    
                    // Skip conjugate root
                    k++;
                }
            }
        }
    }
    
    
    
    /**
     * Constructor, mainly used for testing purposes. The coefficients
     * of the polynomial are provided as a single array of real values.
     *
     * @param coefs The array of coefficients of the polynomial. The degree N
     * of the polynomial is equal to coefs.length-1, coefs[0] is the constant
     * term and coefs[N] is the coefficients of x^N.
     */
    public void init1(double[] coefs) {
        coef = new DoubleDouble[coefs.length];
        for (int i=0; i<coefs.length; i++) {
            coef[i] = new DoubleDouble(coefs[i]);
        }
    }
    
    
    
    /**
     * Constructor, mainly used for testing purposes. With this, a polynomial
     * can be created, based on a set of given real roots and a constant scaling
     * factor.
     *
     * @param coefN The constant scaling factor, which is used to multiply the
     * monic polynomial, derived from the given roots (see below).
     * @param roots A set of roots, which define a monic polynomial, whose
     * degree is equal to the length of the array.
     */
    public void init2(double coefN, double[] roots) {
        coef = new DoubleDouble[roots.length + 1];
        coef[0] = new DoubleDouble(coefN);
        if (roots.length == 0) {
            return;
        }

        for (int k = 0; k < roots.length; k++) {
            DoubleDouble alpha = new DoubleDouble(-roots[k]);
            coef[k + 1] = coef[k];
            for (int j = k; j >= 1; j--) {
                coef[j] = coef[j - 1].add(alpha.mul(coef[j]));
            }
            coef[0] = coef[0].mul(alpha);
        }
    }

    
    
    
    
    /**
     * Constructor, mainly used for testing purposes. This constructor can be
     * used to create a polynomial with a set of real roots and a set of
     * conjugate complex roots.
     *
     * @param coefN The constant scaling factor, which is used to multiuply the
     * monic polynomial, derived from the given roots (see below).
     * @param roots A set of real roots, which define part of the polynomial.
     * @param roots_re The real part of a set of conjugate complex roots.
     * @param roots_im The imaginary part of a set of conjugate complex roots.
     * For each pair of numbers roots_re[k], roots_im[k], two zeros are added!
     * One of them is roots_re[k] + i*roots_im[k] and the other of them is
     * roots_re[k] - i*roots_im[k]. This assures that the polynomial ends up
     * with real coefficients. The length of the arrays roots_re and roots_im
     * must be equal. The total degree of the polynomial is roots.length +
     * 2*roots_re.length.
     */
    public void init3(double coefN, double[] roots, double[] roots_re, double[] roots_im) {
        coef = new DoubleDouble[roots.length + 2 * roots_re.length + 1];
        coef[0] = new DoubleDouble(coefN);
        if (coef.length == 1) {
            return;
        }

        for (int k = 0; k < roots.length; k++) {
            DoubleDouble alpha = new DoubleDouble(-roots[k]);
            coef[k + 1] = coef[k];
            for (int j = k; j >= 1; j--) {
                coef[j] = coef[j - 1].add(alpha.mul(coef[j]));
            }
            coef[0] = coef[0].mul(alpha);
        }

        int nCoefs = 1 + roots.length;
        for (int k = 0; k < roots_re.length; k++) {
            DoubleDouble alpha = new DoubleDouble(roots_re[k]).sqr().add(new DoubleDouble(roots_im[k]).sqr());
            DoubleDouble beta = new DoubleDouble(-2 * roots_re[k]);

            coef[nCoefs + 1] = coef[nCoefs - 1];
            if (nCoefs > 1) {
                coef[nCoefs] = coef[nCoefs - 2].add(coef[nCoefs - 1].mul(beta));
                for (int j = nCoefs - 1; j >= 2; j--) {
                    coef[j] = coef[j - 2].add(coef[j - 1].mul(beta)).add(coef[j].mul(alpha));
                }
                coef[1] = coef[0].mul(beta).add(coef[1].mul(alpha));
            } else {
                coef[1] = coef[0].mul(beta);
            }
            coef[0] = coef[0].mul(alpha);

            nCoefs += 2;
        }
    }
    
    
    
    /**
     * A wrapper method which calls the RPoly-solver to obtain
     * the roots of the polynomial. This tester derives the 
     * coefficients of the polynomial to be solved from the double
     * double precision coefficients, stored in the tester.
     * @param roots_re A pre-allocated array, in which the real part
     * of the roots is stored. This array must have a length of at
     * least the degree of the polynomial.
     * @param roots_im A pre-allocated array, in which the imaginary
     * part of the roots is stored. This array must have a length of at
     * least the degree of the polynomial.
     * @return If all zeros could be isolated, then the degree of the
     * polynomial is returned. If a value less than the degree of the
     * polynomial is returned, then only part of the roots could be
     * isolated, or none at all.
     */
    public int solve(double[] roots_re, double[] roots_im) {
        double[] c = new double[coef.length];
        for (int i=0; i<coef.length; i++) {
            c[i] = coef[i].doubleValue();
        }
        return new RPoly(c).solve(roots_re, roots_im);
    }

    
    
    
    /**
     * Returns the degree of the polynomial.
     *
     * @return The degree of the polynomial.
     */
    int degree() {
        return coef.length - 1;
    }
    
    
    public void testT() {
        
        List<String> names = new ArrayList<>();
        List<String[]> pols = new ArrayList<>();
        PolynomialsUtil.getPolynomials(names, pols);
        for (int i=0; i<names.size(); i++) {
            String[] polynom = pols.get(i);
            if (PolynomialsUtil.hasRealCoefs(polynom)) {
                init0(polynom);
                double[] roots_re = new double[degree()];
                double[] roots_im = new double[degree()];
                int degree = solve(roots_re, roots_im);
                if (degree < degree()) {
                    System.err.println("CONVERGENCE ERROR FOR : " + names.get(i));
                    continue;
                }

                System.out.println(names.get(i));
                System.out.println("---------------");
                for (int j=0; j<degree; j++) {
                    System.out.println(roots_re[j] + "    " + roots_im[j]);
                }
                System.out.println();
            }
        }
        
        {
            // A badly conditioned polynomial with a double root at -2.17 and the
            // other roots also quite close to each other.
            System.out.println("TEST 1");
            System.out.println("------------------");
            init3(6.989607453371678, new double[] {-2.17, -2.17, -1.61, 1.97, 2.59, 2.91},
                                                               new double[] {2.68, 3.01}, new double[] {1.24, 0.42});
            double[] re = new double[degree()];
            double[] im = new double[degree()];
            int n = solve(re, im);
            for (int i = 0; i < n; i++) {
                System.out.println(re[i] + "    " + im[i]);
            }
            System.out.println();
            System.out.println();
        }
        

        {
            // A 20th-degree polynomial with cluster of 5 different but very close
            // roots around -0.01.
            System.out.println("TEST 2");
            System.out.println("------------------");
            //double[] coefs = new double[] {329.3261054405172, -405.11718605351564, 184.92496613009922, -37.74096495644895, 3.005185848918438};
            //double[] coefs = new double[]{-4.060751701902939, 195.100237984308, 395.0668903110914, 92.66756862268753, -131.97108873941974, -12.407124023977762, 9.693065643732627};
            //double[] coefs = new double[] {184.92496613009922*1e20, -184.92496613009922*1e20, -184.92496613009922, 184.92496613009922};
            double[] coefs = new double[] {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
            init1(coefs);
            double[] re = new double[degree()];
            double[] im = new double[degree()];

            int n = solve(re, im);
            for (int i = 0; i < n; i++) {
                System.out.println("" + re[i] + "    " + im[i]);
            }
            System.out.println();
            System.out.println();
        }

        {
            // A badly conditioned polynomial with a double root at -2.17 and the
            // other roots also quite close to each other.
            System.out.println("TEST 3");
            System.out.println("------------------");
            init3(3.8511589483222877, new double[] {-2.08, -3.03, -2.88, -2.85, -4.13, -2.85, -2.75},
                                             new double[] {-0.55, 2.70}, new double[] {-2.93, 0.05});
            double[] re = new double[degree()];
            double[] im = new double[degree()];
            int n = solve(re, im);
            for (int i = 0; i < n; i++) {
                System.out.println(re[i] + "    " + im[i]);
            }
            System.out.println();
            System.out.println();
        }
        
        {
            System.out.println("TEST 4");
            System.out.println("------------------");
            init3(6.989607453371678, new double[] {1, 1e-10, 1e-20, 1e-30, 1e-40, 1e-50},
                                             new double[] {2.68, 3.01}, new double[] {1.24, 0.42});
            double[] re = new double[degree()];
            double[] im = new double[degree()];
            int n = solve(re, im);
            for (int i = 0; i < n; i++) {
                System.out.println(re[i] + "    " + im[i]);
            }
            System.out.println();
            System.out.println();
        }
    }
     
    
    
    
    public void testRnd(String[] args) {
        {
            int NPOLY = 1000000;
            // This is a test with random polynomials whose degree
            // varies between 8 and 23. A million polynomials are
            // chosen at random. Random roots are selected between
            // -5 and 5 for the real and imaginary parts, in steps
            // of 0.01. If computed roots differ more than 0.001 
            // relatively from supplied roots then they are reported
            // as inaccurate. A total count is maintained of the
            // number of polynomials which has at least one root
            // whose approximation is inaccurate.
            System.out.println("TEST WITH " + NPOLY + " RANDOM POLYNOMIALS");
            System.out.println("-------------------------------------------");
            int nErrors = 0;
            int nInaccs = 0;

            for (int i = 0; i < NPOLY; i++) {
                int degr = (int) (4 + 8 * Math.random());
                int degc = (int) (4 + 8 * Math.random());
                if (degc % 2 == 1) {
                    degc++;
                }
                double[] r = new double[degr];
                double[] r_re = new double[degc / 2];
                double[] r_im = new double[degc / 2];
                for (int k = 0; k < degr; k++) {
                    r[k] = ((int) (-500 + 1000 * Math.random())) * 0.01;
                    if (r[k] == 0.0) {
                        r[k] = 0.01;
                    }
                }
                for (int k = 0; k < degc / 2; k++) {
                    r_re[k] = ((int) (-500 + 1000 * Math.random())) * 0.01;
                    if (r_re[k] == 0.0) {
                        r_re[k] = 0.01;
                    }
                    r_im[k] = ((int) (1 + 500 * Math.random())) * 0.01;
                }
                init3(0.1 + 10 * Math.random(), r, r_re, r_im);

                double[] roots_re = new double[degr + degc];
                double[] roots_im = new double[degr + degc];

                int degree = solve(roots_re, roots_im);
                if (degree < degr + degc) {
                    System.err.println("CONVERGENCE ERROR DEGREE = " + degree);
                    for (int k = 0; k <= degr + degc; k++) {
                        System.err.println("   " + coef[k].doubleValue());
                    }
                    nErrors++;
                    continue;
                }

                boolean hasInacc = false;
                for (int k = 0; k < degr; k++) {
                    double diff = distanceFromClosestRoots(r[k], 0.0, roots_re, roots_im);
                    if (diff > 0.001) {
                        System.err.println("INACCURATE REAL ROOT FOR POLY " + i + ": ROOT = " + r[k] + "  DIFF = " + diff);
                        hasInacc = true;
                    }
                }

                for (int k = 0; k < degc / 2; k++) {
                    double diff = distanceFromClosestRoots(r_re[k], r_im[k], roots_re, roots_im);
                    if (diff > 0.001) {
                        System.err.println("INACCURATE CPLX ROOT FOR POLY " + i + ": ROOT = " + r_re[k] + " + i*" + r_im[k] + "  DIFF = " + diff);
                        hasInacc = true;
                    }
                }
                if (hasInacc) {
                    for (int kk = 0; kk < degr; kk++) {
                        System.err.println("    root = " + r[kk]);
                    }
                    for (int kk = 0; kk < degc / 2; kk++) {
                        System.err.println("    root = " + r_re[kk] + " +/- i*" + r_im[kk]);
                    }
                    nInaccs++;
                }
            }

            System.err.println("ERRORS = " + nErrors);
            System.err.println("INACCURACIES = " + nInaccs);
        }
    }
        
    public void testRndCircle(String[] args) {
        {
            Random rnd = new SecureRandom();
            int NPOLY = 1000000;
            // This is a test with random polynomials whose degree
            // varies between 8 and 23. A million polynomials are
            // chosen at random. Random roots are selected on the
            // unit circle in the complex plane. This test is 
            // harder on RPoly than the mainRnd() test.
            System.out.println("TEST WITH " + NPOLY + " RANDOM CIRCLE POLYNOMIALS");
            System.out.println("-------------------------------------------");
            int nErrors = 0;
            int nInaccs = 0;

            for (int i = 0; i < NPOLY; i++) {
                int degc = (int) (8 + 15 * rnd.nextDouble());
                if (degc % 2 == 1) {
                    degc++;
                }
                double[] r_re = new double[degc / 2];
                double[] r_im = new double[degc / 2];
                for (int k = 0; k < degc / 2; k++) {
                    double phi = Math.PI * rnd.nextDouble();
                    double r = 0.9 + 0.2 * rnd.nextDouble();
                    r_re[k] = r*Math.cos(phi);
                    r_im[k] = r*Math.sin(phi);
                }
                init3(0.1 + 10 * rnd.nextDouble(), new double[] {}, r_re, r_im);

                double[] roots_re = new double[degc];
                double[] roots_im = new double[degc];

                int degree = solve(roots_re, roots_im);
                if (degree < degc) {
                    System.err.println("CONVERGENCE ERROR DEGREE = " + degree);
                    for (int k = 0; k <= degc; k++) {
                        System.err.println("   " + coef[k].doubleValue());
                    }
                    nErrors++;
                    continue;
                }

                boolean hasInacc = false;
                for (int k = 0; k < degc / 2; k++) {
                    double diff = distanceFromClosestRoots(r_re[k], r_im[k], roots_re, roots_im);
                    if (diff > 0.001) {
                        System.err.println("INACCURATE CPLX ROOT FOR CIRCLE POLY " + i + ": ROOT = " + r_re[k] + " + i*" + r_im[k] + "  DIFF = " + diff);
                        hasInacc = true;
                    }
                }
                if (hasInacc) {
                    for (int kk = 0; kk < degc / 2; kk++) {
                        System.err.println("    root = " + r_re[kk] + " +/- i*" + r_im[kk]);
                    }
                    nInaccs++;
                }
            }

            System.err.println("ERRORS = " + nErrors);
            System.err.println("INACCURACIES = " + nInaccs);
        }
    }
    
    
    

    public static double distanceFromClosestRoots(double re, double im, double[] rre, double[] rim) {
        double abs = Math.abs(re) + Math.abs(im);
        double minDiff = 1e10;
        for (int i = 0; i < rre.length; i++) {
            double diff = (Math.abs(re - rre[i]) + Math.abs(im - rim[i])) / abs;
            if (diff < minDiff) {
                minDiff = diff;
            }
        }
        return minDiff;
    }

}


/*
This is the result of the original RPoly code as given in the fortran listing.
It fails quite badly for many of the roots. After the modifications of the
RPoly code, as described in the source code, the result is as demonstrated by
TEST1 in the listing above. 

INACCURATE REAL ROOT FOR POLY 414588: ROOT = 2.91  DIFF = 0.18912043531310468
INACCURATE REAL ROOT FOR POLY 414588: ROOT = 1.97  DIFF = 0.015395527012868934
INACCURATE REAL ROOT FOR POLY 414588: ROOT = 2.59  DIFF = 0.15926462830019822
INACCURATE CPLX ROOT FOR POLY 414588: ROOT = 3.01 + i*0.42  DIFF = 0.043392041534168604
INACCURATE CPLX ROOT FOR POLY 414588: ROOT = 2.68 + i*1.24  DIFF = 0.001758006473823295
    coef = 2534.7843646885017
    coef = -65805.59048345718
    coef = 59925.51361433123
    coef = 26138.27852232832
    coef = -45011.37206808903
    coef = 7581.344649940044
    coef = 8909.326647984455
    coef = -4079.4132717432954
    coef = -15.151022596301273
    coef = 375.0553463404709
    coef = -90.44552044662953
    coef = 6.989607453371678
    root = 2.91
    root = 1.97
    root = -2.17
    root = 2.59
    root = 0.04
    root = -1.61
    root = -2.17
    root = 3.01 +/- i*0.42
    root = 2.68 +/- i*1.24
*/
