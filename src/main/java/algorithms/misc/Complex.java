package algorithms.misc;

/**
***********************************************************************
 *  Compilation:  javac Complex.java
 *  Execution:    java Complex
 *
 * Code is from http://introcs.cs.princeton.edu/java/stdlib/
 * "An Introduction to Programming in Java".
 * 
 * The code is released under the GNU General Public License, version 3 (GPLv3). 
 * If you wish to license the code under different terms, please contact our 
 * publisher to discuss.
 * 
 *  Data type for complex numbers.
 *
 *  The data type is "immutable" so once you create and initialize
 *  a Complex object, you cannot change it. The "final" keyword
 *  when declaring re and im enforces this rule, making it a
 *  compile-time error to change the .re or .im fields after
 *  they've been initialized.
 *
 *  % java Complex
 *  a            = 5.0 + 6.0i
 *  b            = -3.0 + 4.0i
 *  Re(a)        = 5.0
 *  Im(a)        = 6.0
 *  b + a        = 2.0 + 10.0i
 *  a - b        = 8.0 + 2.0i
 *  a * b        = -39.0 + 2.0i
 *  b * a        = -39.0 + 2.0i
 *  a / b        = 0.36 - 1.52i
 *  (a / b) * b  = 5.0 + 6.0i
 *  conj(a)      = 5.0 - 6.0i
 *  |a|          = 7.810249675906654
 *  tan(a)       = -6.685231390246571E-6 + 1.0000103108981198i
 *
 * NLK: added a function for nth roots, power, and natural log
 *************************************************************************/

public class Complex {
    private final double re;   // the real part
    private final double im;   // the imaginary part

    // create a new object with the given real and imaginary parts

    /**
     *
     @param real
     @param imag
     */
    public Complex(double real, double imag) {
        re = real;
        im = imag;
    }

    // return a string representation of the invoking Complex object
    public String toString() {
        if (im == 0) return re + "";
        if (re == 0) return im + "i";
        if (im <  0) return re + " - " + (-im) + "i";
        return re + " + " + im + "i";
    }

    // return abs/modulus/magnitude and angle/phase/argument

    /**
     *
     @return
     */
    public double abs()   { return Math.hypot(re, im); }  // Math.sqrt(re*re + im*im)

    /**
     *
     @return
     */
    public double phase() { return Math.atan2(im, re); }  // between -pi and pi

    // return a new Complex object whose value is (this + b)

    /**
     *
     @param b
     @return
     */
    public Complex plus(Complex b) {
        Complex a = this;             // invoking object
        double real = a.re + b.re;
        double imag = a.im + b.im;
        return new Complex(real, imag);
    }
    
    // return a new Complex object whose value is (this + b)

    /**
     *
     @param b
     @return
     */
    public Complex plus(double b) {
        Complex a = this;             // invoking object
        double real = a.re + b;
        double imag = a.im;
        return new Complex(real, imag);
    }

    // return a new Complex object whose value is (this - b)

    /**
     *
     @param b
     @return
     */
    public Complex minus(Complex b) {
        Complex a = this;
        double real = a.re - b.re;
        double imag = a.im - b.im;
        return new Complex(real, imag);
    }
    
    /**
     *
     @return
     */
    public Complex copy() {
        return new Complex(this.re, this.im);
    }

    // return a new Complex object whose value is (this * b)

    /**
     *
     @param b
     @return
     */
    public Complex times(Complex b) {
        Complex a = this;
        double real = a.re * b.re - a.im * b.im;
        double imag = a.re * b.im + a.im * b.re;
        return new Complex(real, imag);
    }

    // scalar multiplication
    // return a new object whose value is (this * alpha)

    /**
     *
     @param alpha
     @return
     */
    public Complex times(double alpha) {
        return new Complex(alpha * re, alpha * im);
    }

    // return a new Complex object whose value is the conjugate of this

    /**
     *
     @return
     */
    public Complex conjugate() {  return new Complex(re, -im); }

    // return a new Complex object whose value is the reciprocal of this

    /**
     *
     @return
     */
    public Complex reciprocal() {
        double scale = re*re + im*im;
        return new Complex(re / scale, -im / scale);
    }

    // return the real or imaginary part

    /**
     *
     @return
     */
    public double re() { return re; }

    /**
     *
     @return
     */
    public double im() { return im; }

    // return a / b

    /**
     *
     @param b
     @return
     */
    public Complex divided(Complex b) {
        Complex a = this;
        return a.times(b.reciprocal());
    }

    // return a new Complex object whose value is the complex exponential of this

    /**
     *
     @return
     */
    public Complex exp() {
        return new Complex(Math.exp(re) * Math.cos(im), Math.exp(re) * Math.sin(im));
    }

    // return a new Complex object whose value is the complex sine of this

    /**
     *
     @return
     */
    public Complex sin() {
        return new Complex(Math.sin(re) * Math.cosh(im), Math.cos(re) * Math.sinh(im));
    }

    // return a new Complex object whose value is the complex cosine of this

    /**
     *
     @return
     */
    public Complex cos() {
        return new Complex(Math.cos(re) * Math.cosh(im), -Math.sin(re) * Math.sinh(im));
    }

    // return a new Complex object whose value is the complex tangent of this

    /**
     *
     @return
     */
    public Complex tan() {
        return sin().divided(cos());
    }
    
    // a static version of plus

    /**
     *
     @param a
     @param b
     @return
     */
    public static Complex plus(Complex a, Complex b) {
        double real = a.re + b.re;
        double imag = a.im + b.im;
        Complex sum = new Complex(real, imag);
        return sum;
    }

    /**
     * get the (1/n) power of this instance
     * from Boas "mathematical methods in the physical sciences"
       chap 2, section 13
     @param n the number to take the (1/n) power of this instance
     @return 
     */
    public Complex nthRoot(double n) {
        // r^(1/n) * (cos(theta/n) + i*sin(theta/n))
        double rnth = Math.pow(abs(), 1./n);
        double theta = phase();
        return new Complex(rnth*Math.cos(theta/n), rnth*Math.sin(theta/n));
    }
    
    /**
     * get the natural log of this instance
     * from Boas "mathematical methods in the physical sciences"
       chap 2, section 13
     @return 
     */
    public Complex naturalLog() {
        double r = abs();
        double theta = phase();
        return new Complex(Math.log(r), theta);
    }
    
    /**
     * get this instance raised to the power b
     * from Boas "mathematical methods in the physical sciences"
       chap 2, section 13
     @param b the power to apply to this instance
     @return 
     */
    public Complex power(Complex b) {
        // a^b = e^(b*ln(a))
        Complex blna = b.times(naturalLog());
        return blna.exp();
    }
}
