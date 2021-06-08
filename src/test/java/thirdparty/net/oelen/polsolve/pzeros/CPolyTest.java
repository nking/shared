package thirdparty.net.oelen.polsolve.pzeros;

import java.util.ArrayList;
import java.util.List;
import thirdparty.net.oelen.polarith.DoubleDouble;
import thirdparty.net.oelen.polsolve.jt.CPoly;
import junit.framework.TestCase;

// This class is used for testing CPoly. Although CPoly works at
// standard 53 bit precision, computation of coefficients from given
// sets of roots is done at 105 bit precision and only when calling
// the solver, the coefficients are truncated at 53 bits.
// For the tests, many polynomials are specified by means of its
// zeros and coefficients are constructed and the excess precision
// in this process assures that the solver does get coefficients
// at maximum attainable precision and that certain issues with the
// found roots are not due to excessive introduction of numerical
// noise into the coefficients.

public class CPolyTest extends TestCase {
    
    
    private  DoubleDouble[] coef_re;
    private  DoubleDouble[] coef_im;

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
            coef_re = new DoubleDouble[degree+1];
            coef_im = new DoubleDouble[degree+1];
        }
        else if (testPolynomial[0].equals("RC")) {
            // Roots Complex
            for (int i=1; i<testPolynomial.length; i++) {
                String[] num = testPolynomial[i].split(";");
                roots_re.add(new DoubleDouble(num[0]));
                roots_im.add((num.length >= 2) ? new DoubleDouble(num[1]) : DoubleDouble.ZERO);
            }
            int degree = roots_re.size();
            coef_re = new DoubleDouble[degree+1];
            coef_im = new DoubleDouble[degree+1];
        }
        else if (testPolynomial[0].equals("CR")) {
            // Coefs Real
            int ncoefs = testPolynomial.length - 1;
            coef_re = new DoubleDouble[ncoefs];
            coef_im = new DoubleDouble[ncoefs];
            for (int i=1; i<testPolynomial.length; i++) {
                coef_re[i-1] = new DoubleDouble(testPolynomial[i]);
                coef_im[i-1] = DoubleDouble.ZERO;
            }
        }
        else if (testPolynomial[0].equals("CC")) {
            // Coefs Complex
            int ncoefs = testPolynomial.length - 1;
            coef_re = new DoubleDouble[ncoefs];
            coef_im = new DoubleDouble[ncoefs];
            for (int i=1; i<testPolynomial.length; i++) {
                String[] num = testPolynomial[i].split(";");
                coef_re[i-1] = new DoubleDouble(num[0]);
                coef_im[i-1] = (num.length >= 2) ? new DoubleDouble(num[1]) : DoubleDouble.ZERO;
            }
        }
        else {
            coef_re = coef_im = null;
        }
        
        if (!roots_re.isEmpty()) {
            coef_re[0] = DoubleDouble.ONE;
            coef_im[0] = DoubleDouble.ZERO;

            for (int k=0; k<roots_re.size(); k++) {
                DoubleDouble alpha_re = roots_re.get(k).neg();
                DoubleDouble alpha_im = roots_im.get(k).neg();
                coef_re[k+1] = coef_re[k];
                coef_im[k+1] = coef_im[k];
                for (int j=k; j>=1; j--) {
                    DoubleDouble re_j = coef_re[j-1].add(alpha_re.mul(coef_re[j])).sub(alpha_im.mul(coef_im[j]));
                    DoubleDouble im_j = coef_im[j-1].add(alpha_re.mul(coef_im[j])).add(alpha_im.mul(coef_re[j]));
                    coef_re[j] = re_j;
                    coef_im[j] = im_j;
                }
                DoubleDouble re_0 = alpha_re.mul(coef_re[0]).sub(alpha_im.mul(coef_im[0]));
                DoubleDouble im_0 = alpha_re.mul(coef_im[0]).add(alpha_im.mul(coef_re[0]));
                coef_re[0] = re_0;
                coef_im[0] = im_0;
            }
        }
    }
    
    
    
    /**
     * Constructor, mainly used for testing purposes. This constructor can be
     * used to create a polynomial with a set of real roots and a set of
     * conjugate complex roots.
     *
     * @param coefN_re The real part of the constant scaling factor, which is 
     * used to multiuply the monic polynomial, derived from the given roots (see below).
     * @param coefN_im The imaginary part of the constant scaling factor, which
     * is used to multiuply the monic polynomial, derived from the given roots.
     * @param roots_re The real part of the set of complex roots.
     * @param roots_im The imaginary part of the set of complex roots.
     * The parameter roots_re and roots_im must have the same length. The degree
     * of the derived polynomial is equal to the length of these arrays.
     */
    public void init1(double coefN_re, double coefN_im, double[] roots_re, double[] roots_im) {
        coef_re = new DoubleDouble[roots_re.length + 1];
        coef_im = new DoubleDouble[roots_re.length + 1];
        coef_re[0] = new DoubleDouble(coefN_re);
        coef_im[0] = new DoubleDouble(coefN_im);
        if (roots_re.length == 0) {
            return;
        }
        
        for (int k=0; k<roots_re.length; k++) {
            DoubleDouble alpha_re = new DoubleDouble(-roots_re[k]);
            DoubleDouble alpha_im = new DoubleDouble(-roots_im[k]);
            coef_re[k+1] = coef_re[k];
            coef_im[k+1] = coef_im[k];
            for (int j=k; j>=1; j--) {
                DoubleDouble re_j = coef_re[j-1].add(alpha_re.mul(coef_re[j])).sub(alpha_im.mul(coef_im[j]));
                DoubleDouble im_j = coef_im[j-1].add(alpha_re.mul(coef_im[j])).add(alpha_im.mul(coef_re[j]));
                coef_re[j] = re_j;
                coef_im[j] = im_j;
            }
            DoubleDouble re_0 = alpha_re.mul(coef_re[0]).sub(alpha_im.mul(coef_im[0]));
            DoubleDouble im_0 = alpha_re.mul(coef_im[0]).add(alpha_im.mul(coef_re[0]));
            coef_re[0] = re_0;
            coef_im[0] = im_0;
        }
    }
    
    
    
    
    
    
    
    /**
     * A wrapper method which calls the CPoly-solver to obtain
     * the roots of the polynomial. This tester derives the 
     * coefficients of the polynomial to be solved from the double
     * double precision coefficients, stored in the tester.
     * @param r_re A pre-allocated array, in which the real part
     * of the roots is stored. This array must have a length of at
     * least the degree of the polynomial.
     * @param r_im A pre-allocated array, in which the imaginary
     * part of the roots is stored. This array must have a length of at
     * least the degree of the polynomial.
     * @return If all zeros could be isolated, then the degree of the
     * polynomial is returned. If a value less than the degree of the
     * polynomial is returned, then only part of the roots could be
     * isolated, or none at all.
     */
    public int solve(double[] r_re, double[] r_im) {
        double[] cr = new double[coef_re.length];
        double[] ci = new double[coef_im.length];
        for (int i=0; i<coef_re.length; i++) {
            cr[i] = coef_re[i].doubleValue();
        }
        for (int i=0; i<coef_im.length; i++) {
            ci[i] = coef_im[i].doubleValue();
        }
        CPoly pol = new CPoly(cr, ci);
        return pol.solve(r_re, r_im);
    }

    
    
    
    /**
     * Returns the degree of the polynomial.
     *
     * @return The degree of the polynomial.
     */
    int degree() {
        return coef_re.length - 1;
    }
    
    public void testWilco() {
        mainWilco();
    }
    public void testT() {
        mainT();
    }
    public void testRnd() {
        mainRnd();
    }
    public void test2() {
        main2();
    }
    
    public static void main(String[] args) {
        //mainWilco();
        //mainT();    // Polynomials from TestPolynomials.java
      //mainRnd();  // Random polynomials.
        //main2();
    }    
       
    // Tests polynomials from the TestPolynomials.java suite.
    public void mainT() {
        List<String> names = new ArrayList<>();
        List<String[]> pols = new ArrayList<>();
        PolynomialsUtil.getPolynomials(names, pols);
        for (int i=0; i<names.size(); i++) {
            init0(pols.get(i));
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
    
    
    
    
    
    // A tester main method, which solves 1000000 equations with
    // randomly chosen roots and which compares the found roots
    // with the original roots from which the polynomial is derived.
    // Of the random roots, 500000 are chosen from a square in 
    // the complex plane [-5, 5] + i*[-5, 5]. The other 500000 are
    // chosen randomly from a ring with radius between 0.9 and 1.1.
    //
    // Polynomials are constructed from a set of known roots. The
    // polynomials are solved and then the computed roots are 
    // compared to the known roots. If the relative error of a
    // computed root is more than 0.001 times the absolute value
    // of the known root, then a warning is issued. At the end of
    // the test, the software tells how many polynomials there were
    // with at least one inaccurate root.
    public void mainRnd() {
        int nErrors = 0;
        int nInaccs = 0;
        int nIterMax = 50000;//500000;
        for (int i=0; i<nIterMax; i++) {
            int degr = (int)(4 + 8*Math.random());
            int degc = (int)(4 + 8*Math.random());
            double[] r_re = new double[degr + degc];
            double[] r_im = new double[degr + degc];
            for (int k=0; k<degr; k++) {
                r_re[k] = ((int)(-500 + 1000*Math.random()))*0.01;
                if (r_re[k] == 0.0) r_re[k] = 0.005;
                r_im[k] = 0.0;
            }
            for (int k=0; k<degc; k++) {
                r_re[degr + k] = ((int)(-500 + 1000*Math.random()))*0.01;
                if (r_re[degr + k] == 0.0) r_re[degr + k] = 0.005;
                r_im[degr + k] = ((int)(-500 + 1000*Math.random()))*0.01;
                if (r_im[degr + k] == 0.0) r_im[degr + k] = 0.005;
            }
            init1(0.1+10*Math.random(), 0.1+10*Math.random(), r_re, r_im);
            
            double[] roots_re = new double[degr + degc];
            double[] roots_im = new double[degr + degc];
            
            int degree = solve(roots_re, roots_im);
            if (degree < degr + degc) {
                System.err.println("CONVERGENCE ERROR DEGREE = " + degree);
                for (int k=0; k<=degr+degc; k++) {
                    System.err.println("   " + coef_re[k] + "   " + coef_im[k]);
                }
                nErrors++;
                continue;
            }
            
            boolean hasInacc = false;
            for (int k=0; k<degr+degc; k++) {
                double diff = distanceFromClosestRoot(r_re[k], r_im[k], roots_re, roots_im);
                if (diff > 0.001) {
                    System.err.println("INACCURATE CPLX ROOT FOR POLY " + i + ": ROOT = " + r_re[k] + " + i*" + r_im[k] + "  DIFF = " + diff);
                    hasInacc = true;
                }
            }
            if (hasInacc) {
                for (int kk=0; kk<degr+degc; kk++) {
                    System.err.println("    root = " + r_re[kk] + " +/- i*" + r_im[kk]);
                }
                nInaccs++;
            }
        }
        
        
        for (int i=0; i<nIterMax; i++) {
            int deg = (int)(8 + 15*Math.random());
            double[] r_re = new double[deg];
            double[] r_im = new double[deg];
            for (int k=0; k<deg; k++) {
                double phi = 2*Math.PI*Math.random();
                double r = 0.9 + 0.2*Math.random();
                r_re[k] = r*Math.cos(phi);
                r_im[k] = r*Math.sin(phi);
            }
            init1(0.1+10*Math.random(), 0.1+10*Math.random(), r_re, r_im);
            
            double[] roots_re = new double[deg];
            double[] roots_im = new double[deg];
            
            int degree = solve(roots_re, roots_im);
            if (degree < deg) {
                System.err.println("CONVERGENCE CIRCLE ERROR DEGREE = " + degree);
                for (int k=0; k<=deg; k++) {
                    System.err.println("   " + coef_re[k] + "   " + coef_im[k]);
                }
                nErrors++;
                continue;
            }
            
            boolean hasInacc = false;
            for (int k=0; k<deg; k++) {
                double diff = distanceFromClosestRoot(r_re[k], r_im[k], roots_re, roots_im);
                if (diff > 0.001) {
                    System.err.println("INACCURATE CIRCLE ROOT FOR POLY " + i + ": ROOT = " + r_re[k] + " + i*" + r_im[k] + "  DIFF = " + diff);
                    hasInacc = true;
                }
            }
            if (hasInacc) {
                for (int kk=0; kk<deg; kk++) {
                    System.err.println("    root = " + r_re[kk] + " +/- i*" + r_im[kk]);
                }
                nInaccs++;
            }
        }
        
        System.err.println("ERRORS = " + nErrors);
        System.err.println("INACCURACIES = " + nInaccs);
    }

    
    
    
    public void mainWilco() {
        double[] rre = {-3.02, -3.86,-1.94, -0.41, -2.27, 4.97, -2.39, -2.37, -1.68, -2.17, -2.31, -3.13, -3.34, -3.03, 0.3};
        double[] rim = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.28, 1.12, -3.93, -4.87};
            
            init1(9.643264256234523, 0.856475436745643, rre, rim);
            double[] rootsr = new double[degree()];
            double[] rootsi = new double[degree()];
            int result = solve(rootsr, rootsi);
            for (int i=0; i<rootsr.length; i++) {
                System.out.println(rootsr[i] + "  i*" + rootsi[i]);
            }
}

    

    /**
     * A set of test polynomials, which demonstrate correct working of 
     * polynomial solver code. It does not contain particularly hard
     * polynomials, each solver should solve these without issue.
     */
    public void main2() {

        System.out.println();
        System.out.println("Testing Polynomial...");

        try {
            // Evaluate p(x) = x^2 + 2*x + 1.
            double[] coef = new double[3];
            double[] zeros_re = new double[coef.length-1];
            double[] zeros_im = new double[coef.length-1];
            coef[0] = 1.;
            coef[1] = 2.;
            coef[2] = 1.;
            CPoly poly = new CPoly(coef);

            System.out.println("    EXAMPLE 0:  A VERY SIMPLE POLYNOMIAL.");
            System.out.println("    p(x) = \n" + poly);
            System.out.println("    Zeros of p(x) are:");
            int deg = poly.solve(zeros_re, zeros_im);
            for (int i = 0; i < deg; ++i) {
                System.out.println("        " + zeros_re[i] + "    " + zeros_im[i]);
            }

            System.out.println("\n\n    EXAMPLE 1.  POLYNOMIAL WITH ZEROS 1,2,...,10.");
            coef = new double[11];
            zeros_re = new double[coef.length-1];
            zeros_im = new double[coef.length-1];
            coef[10] = 1;
            coef[9] = -55;
            coef[8] = 1320;
            coef[7] = -18150;
            coef[6] = 157773;
            coef[5] = -902055;
            coef[4] = 3416930;
            coef[3] = -8409500;
            coef[2] = 12753576;
            coef[1] = -10628640;
            coef[0] = 3628800;
            poly = new CPoly(coef);
            System.out.println("    p(x) = \n" + poly);
            deg = poly.solve(zeros_re, zeros_im);
            for (int i = 0; i < deg; ++i) {
                System.out.println("        " + zeros_re[i] + "    " + zeros_im[i]);
            }

            System.out.println("\n\n    EXAMPLE 2. ZEROS ON IMAGINARY AXIS DEGREE 3.");
            double[] coef_re = new double[4];
            double[] coef_im = new double[4];
            zeros_re = new double[coef_re.length-1];
            zeros_im = new double[coef_re.length-1];
            coef_re[0] = 0.0;  coef_im[0] = 1.0;
            coef_re[1] = -10001.0001;  coef_im[1] = 0.0;
            coef_re[2] = 0.0;  coef_im[2] = -10001.0001;
            coef_re[3] = 1.0;  coef_im[3] = 0.0;
            poly = new CPoly(coef_re, coef_im);
            System.out.println("    p(x) = \n" + poly);
            deg = poly.solve(zeros_re, zeros_im);
            for (int i = 0; i < deg; ++i) {
                System.out.println("        " + zeros_re[i] + "    " + zeros_im[i]);
            }

            System.out.println("\n\n    EXAMPLE 3. ZEROS AT 1+I,1/2*(1+I)....1/(2**-9)*(1+I).");
            coef_re = new double[] {0, -4.652065399568528E-10, 1.584803612786345E-7, -1.154642632172909E-5,
                                    0, 1.271507365163416E-2, -0.2002119533717632, 0.7567065954208374, 0, -1.998046875, 1};
            coef_im = new double[] {9.094947017729282E-13, -4.652065399568528E-10, 0, 1.154642632172909E-5, -7.820779428584501E-4,
                                    1.271507365163416E-2, 0., -7.567065954208374E-1, 2.658859252929688, -1.998046875, 0};
            poly = new CPoly(coef_re, coef_im);
            zeros_re = new double[coef_re.length-1];
            zeros_im = new double[coef_re.length-1];
            System.out.println("    p(x) = \n" + poly);
            deg = poly.solve(zeros_re, zeros_im);
            for (int i = 0; i < deg; ++i) {
                System.out.println("        " + zeros_re[i] + "    " + zeros_im[i]);
            }

            System.out.println("\n\n    EXAMPLE 4. MULTIPLE ZEROS.");
            coef_re = new double[] {288, -1344, 2204, -920, -1587, 2374, -1293, 284, 3, -10, 1};
            coef_im = new double[] {0, 504, -2352, 4334, -3836, 1394, 200, -334, 100, -10, 0};
            poly = new CPoly(coef_re, coef_im);
            zeros_re = new double[coef_re.length-1];
            zeros_im = new double[coef_re.length-1];
            System.out.println("    p(x) = \n" + poly);
            deg = poly.solve(zeros_re, zeros_im);
            for (int i = 0; i < deg; ++i) {
                System.out.println("        " + zeros_re[i] + "    " + zeros_im[i]);
            }

            System.out.println("\n\n    EXAMPLE 5. 12 ZEROS EVENLY DISTRIBUTED ON A CIRCLE OF RADIUS"
                    + " 1 CENTERED AT 0+2I.");
            coef_re = new double[] {4095, 0, -67584, 0, 126720, 0, -59136, 0, 7920, 0, -264, 0, 1};
            coef_im = new double[] {0, 24576, 0, -112640, 0, 101376, 0, -25344, 0, 1760, 0, -24, 0};
            poly = new CPoly(coef_re, coef_im);
            zeros_re = new double[coef_re.length-1];
            zeros_im = new double[coef_re.length-1];
            System.out.println("    p(x) = \n" + poly);
            deg = poly.solve(zeros_re, zeros_im);
            for (int i = 0; i < deg; ++i) {
                System.out.println("        " + zeros_re[i] + "    " + zeros_im[i]);
            }

            System.out.println("\n\n    EXAMPLE 6. Zeros of the polynomial 1 + 2x + 3x^2 + ... 51*x^50");
            coef = new double[51];
            for (int i=0; i<=50; i++) {
                coef[i] = 1+i;
            }
            poly = new CPoly(coef);
            zeros_re = new double[coef.length-1];
            zeros_im = new double[coef.length-1];
            System.out.println("    p(x) = \n" + poly);
            deg = poly.solve(zeros_re, zeros_im);
            for (int i = 0; i < deg; ++i) {
                System.out.println("        " + zeros_re[i] + "    " + zeros_im[i]);
            }

            System.out.println("\n\n    EXAMPLE 7. Zeros of the polynomial 1 + 2x + 3x^2 + ... 6*x^5");
            coef = new double[6];
            for (int i=0; i<=5; i++) {
                coef[i] = 1+i;
            }
            poly = new CPoly(coef);
            zeros_re = new double[coef.length-1];
            zeros_im = new double[coef.length-1];
            System.out.println("    p(x) = \n" + poly);
            deg = poly.solve(zeros_re, zeros_im);
            for (int i = 0; i < deg; ++i) {
                System.out.println("        " + zeros_re[i] + "    " + zeros_im[i]);
            }

        } 
        catch (Exception e) {
            e.printStackTrace();
        }    
    }
    
    
    
    
    
    
    public static double distanceFromClosestRoot(double re, double im, double[] rre, double[] rim) {
        double abs = Math.abs(re) + Math.abs(im);
        double minDiff = 1e10;
        for (int i=0; i<rre.length; i++) {
            double diff = (Math.abs(re - rre[i]) + Math.abs(im - rim[i]))/abs;
            if (diff < minDiff) 
                minDiff = diff;
        }
        return minDiff;
    }
}
