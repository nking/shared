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
package thirdparty.net.oelen.polarith;

import java.io.Serializable;
import java.util.Arrays;

/**
 * Immutable, extended-precision floating-point numbers which maintain 105 bits
 * (approximately 30 decimal digits) of precision.
 * <p>
 * A DoubleDouble uses a representation containing two double-precision values.
 * A number x is represented as a pair of doubles, x.hi and x.lo, such that the
 * number represented by x is x.hi + x.lo, where
 * <pre>
 *    |x.lo| <= 0.5*ulp(x.hi)
 * </pre> and ulp(y) means "unit in the last place of y". The basic arithmetic
 * operations are implemented using convenient properties of IEEE-754
 * floating-point arithmetic.
 * <p>
 * The range of values which can be represented is the same as in IEEE-754. The
 * precision of the representable numbers is twice as great as IEEE-754 double
 * precision.
 * <p>
 * The correctness of the arithmetic algorithms relies on operations being
 * performed with standard IEEE-754 double precision and rounding. This is the
 * Java standard arithmetic model, but for performance reasons Java
 * implementations are not constrained to using this standard by default. Some
 * processors (notably the Intel Pentium architecure) perform floating point
 * operations in (non-IEEE-754-standard) extended-precision. A JVM
 * implementation may choose to use the non-standard extended-precision as its
 * default arithmetic mode. To prevent this from happening, this code uses the
 * Java <tt>strictfp</tt> modifier, which forces all operations to take place in
 * the standard IEEE-754 rounding model.
 * <p>
 * The API provides a value-oriented interface. DoubleDouble values are
 * immutable; operations on them return new objects carrying the result of the
 * operation. This provides a much simpler semantics for writing DoubleDouble
 * expressions, and Java memory management is efficient enough that this imposes
 * very little performance penalty.
 * <p>
 * This implementation uses algorithms originally designed variously by Knuth,
 * Kahan, Dekker, and Linnainmaa. Douglas Priest developed the first C
 * implementation of these techniques. Other more recent C++ implementation are
 * due to Keith M. Briggs and David Bailey et al.
 *
 * <h3>References</h3>
 * <ul>
 * <li>Priest, D., <i>Algorithms for Arbitrary Precision Floating Point
 * Arithmetic</i>, in P. Kornerup and D. Matula, Eds., Proc. 10th Symposium on
 * Computer Arithmetic, IEEE Computer Society Press, Los Alamitos, Calif., 1991.
 * <li>Yozo Hida, Xiaoye S. Li and David H. Bailey,
 * <i>Quad-Double Arithmetic: Algorithms, Implementation, and Application</i>,
 * manuscript, Oct 2000; Lawrence Berkeley National Laboratory Report BNL-46996.
 * <li>David Bailey, <i>High Precision Software Directory</i>;
 * <tt>http://crd.lbl.gov/~dhbailey/mpdist/index.html</tt>
 * </ul>
 *
 *
 * @author Martin Davis
 *
 */
public final strictfp class DoubleDouble implements Serializable, Comparable, Cloneable {
    
    public static final long serialVersionUID = Hash64.hash("DoubleDouble_v1.0");

    /**
     * The value nearest to the constant Pi.
     */
    public static final DoubleDouble PI = new DoubleDouble(
            3.141592653589793116e+00,
            1.224646799147353207e-16);

    /**
     * The value nearest to the constant 2 * Pi.
     */
    public static final DoubleDouble TWO_PI = new DoubleDouble(
            6.283185307179586232e+00,
            2.449293598294706414e-16);

    /**
     * The value nearest to the constant Pi / 2.
     */
    public static final DoubleDouble PI_2 = new DoubleDouble(
            1.570796326794896558e+00,
            6.123233995736766036e-17);

    /**
     * The value nearest to the constant e (the natural logarithm base).
     */
    public static final DoubleDouble E = new DoubleDouble(
            2.718281828459045091e+00,
            1.445646891729250158e-16);
    
    
    
    /**
     * A few simple numerical constants.
     */
    public static final DoubleDouble ZERO = new DoubleDouble();
    public static final DoubleDouble ONE = new DoubleDouble(1.0);
    public static final DoubleDouble TWO = new DoubleDouble(2.0);
    public static final DoubleDouble TEN = new DoubleDouble(10.0);
    public static final DoubleDouble MAXINT = new DoubleDouble(Integer.MAX_VALUE);
    public static final DoubleDouble MAXLONG = new DoubleDouble(Long.MAX_VALUE);
    public static final DoubleDouble MININT = new DoubleDouble(Integer.MIN_VALUE);
    public static final DoubleDouble MINLONG = new DoubleDouble(Long.MIN_VALUE);

    /**
     * A value representing the result of an operation which does not return a
     * valid number.
     */
    public static final DoubleDouble NaN = new DoubleDouble(Double.NaN, Double.NaN);

    /**
     * The smallest representable relative difference between two {link @
     * DoubleDouble} values
     */
    public static final double EPS = 1.23259516440783e-32;  /* = 2^-106 */


    
    
    /**
     * Converts the string argument to a DoubleDouble number.
     *
     * @param str a string containing a representation of a numeric value
     * @return the extended precision version of the value
     * @throws NumberFormatException if <tt>s</tt> is not a valid representation
     * of a number
     */
    public static DoubleDouble valueOf(String str)
            throws NumberFormatException {
        return parse(str);
    }

    
    
    
    /**
     * Converts the <tt>double</tt> argument to a DoubleDouble number.
     *
     * @param x a numeric value
     * @return the extended precision version of the value
     */
    public static DoubleDouble valueOf(double x) {
        return new DoubleDouble(x);
    }

    
    
    
    /**
     * Converts the <tt>long</tt> argument to a DoubleDouble number.
     *
     * @param x a numeric value
     * @return the extended precision version of the value
     */
    public static DoubleDouble valueOf(long x) {
        return new DoubleDouble(x);
    }

    
    
    
    /**
     * Converts the <tt>int</tt> argument to a DoubleDouble number.
     *
     * @param x a numeric value
     * @return the extended precision version of the value
     */
    public static DoubleDouble valueOf(int x) {
        return new DoubleDouble(x);
    }

    
    
    
    /**
     * The value to split a double-precision value on during multiplication
     */
    private static final double SPLIT = 0x08000001; // 2^27+1, for IEEE double

    
    
    
    /**
     * The high-order component of the double-double precision value.
     */
    double hi;  // Not private, must be accessible from same package.

    
    
    
    /**
     * The low-order component of the double-double precision value.
     */
    double lo;  // Not private, must be accessible from same package.

    
    
    
    /**
     * Creates a new DoubleDouble with value 0.0.
     */
    public DoubleDouble() {
        hi = lo = 0.0;
    }
    
    
    

    /**
     * Creates a new DoubleDouble with value x.
     *
     * @param x the value to initialize
     */
    public DoubleDouble(double x) {
        hi = x;
        lo = 0.0;
    }

    
    /**
     * Creates a new DoubleDouble with value i.
     *
     * @param i the value to initialize
     */
    public DoubleDouble(int i) {
        hi = i;
        lo = 0;
    }
    

    
    /**
     * Creates a new DoubleDouble with value l.
     *
     * @param l the value to initialize
     */
    public DoubleDouble(long l) { 
        hi = l & 0xfffffffffffff800l;
        lo = l & 0x7ffl;
        double s = hi + lo;
        double err = lo - (s - hi);
        hi = s;
        lo = err;
    }
    
    
    

    /**
     * Creates a new DoubleDouble with value (hi, lo).
     *
     * @param hi the high-order component
     * @param lo the high-order component
     */
    private DoubleDouble(double hi, double lo) {
        // This method does no normalization of a number, it assumes
        // that hi and lo are correct for further operations. For this
        // reason it only is made available for internal use in this
        // library.
        this.hi = hi;
        this.lo = lo;
    }
    
    
    

    /**
     * Creates a new DoubleDouble with value equal to the argument.
     *
     * @param str the value to initialize by
     * @throws NumberFormatException if <tt>str</tt> is not a valid
     * representation of a number
     */
    public DoubleDouble(String str)
            throws NumberFormatException {
        DoubleDouble val = parse(str);
        this.hi = val.hi;
        this.lo = val.lo;
    }
    
    
    

    /**
     * Creates and returns a copy of this value.
     *
     * @return a copy of this value
     * @throws CloneNotSupportedException but in practice this never
     * will occur for a DoubleDouble object
     */
    @Override
    public Object clone() throws CloneNotSupportedException {
        DoubleDouble cp = (DoubleDouble)super.clone();
        cp.hi = hi;
        cp.lo = lo;
        return cp;
    }
    
    
    

    private void RENORM() {
        double s = hi + lo;
        double err = lo - (s - hi);
        hi = s;
        lo = err;
    }
    
    
    
    
    /**
     * Returns a double[] array with approximate values of the supplied
     * DoubleDouble array. Null values are converted to 0.0.
     *
     * @param arr The DoubleDouble[] array to be converted to double[].
     * @return A double[] array with approximate values of the input array.
     */
    public static double[] toDouble(DoubleDouble[] arr) {
        double[] a = new double[arr.length];
        for (int i=0; i<arr.length; i++) {
            a[i] = arr[i] == null ? 0.0 : arr[i].hi + arr[i].lo;
        }
        return a;
    }
    
    
    
    
    /**
     * Returns a DoubleDouble[] array with values of the supplied
     * double array.
     *
     * @param arr The double[] array to be converted to DoubleDouble[].
     * @return A DoubleDouble[] array with values of the input array.
     */
    public static DoubleDouble[] toDoubleDouble(double[] arr) {
        DoubleDouble[] a = new DoubleDouble[arr.length];
        for (int i=0; i<arr.length; i++) {
            a[i] = new DoubleDouble(arr[i]);
        }
        return a;
    }
    
    
    
    
    /**
     * Converts a DoubleDouble[] array to a double[] array. The
     * output array can optionally be supplied. If no array is
     * supplied (null value is given), then the needed array is
     * created. Null values are converted to 0.0.
     *
     * @param a The output double[] array, in which the approximate converted
     * values are stored. A null-value may be supplied. In that case an output
     * array is created. If the output array is too short, then only part of
     * the DoubleDouble values is converted.
     * @param arr The DoubleDouble[] array to be converted to double[].
     * @return A double[] array with approximate values of the input array.
     * This can be the supplied input array or a newly allocated one.
     */
    public static double[] toDouble(double[] a, DoubleDouble[] arr) {
        if (a == null) a = new double[arr.length];
        for (int i=0; i<a.length && i<arr.length; i++) {
            a[i] = arr[i] == null ? 0.0 : arr[i].hi + arr[i].lo;
        }
        return a;
    }
    
    
    
    
    /**
     * Converts a double[] array to a DoubleDouble[] array. The
     * output array can optionally be supplied. If no array is
     * supplied (null value is given), then the needed array is
     * created.
     *
     * @param a The output DoubleDouble[] array, in which the converted
     * values are stored. A null-value may be supplied. In that case an
     * output array is created. If the output array is too short, then 
     * only part of the values is converted.
     * @param arr The double[] array to be converted to DoubleDouble[].
     * @return A DoubleDouble[] array with values of the input array.
     * This can be the supplied input array or a newly allocated one.
     */
    public static DoubleDouble[] toDoubleDouble(DoubleDouble[] a, double[] arr) {
        if (a == null) a = new DoubleDouble[arr.length];
        for (int i=0; i<a.length && i<arr.length; i++) {
            a[i] = new DoubleDouble(arr[i]);
        }
        return a;
    }
     
     
     
     
    /**
     * Returns a DoubleDouble whose value is <tt>(this + y)</tt>.
     *
     * @param y the addend
     * @return <tt>(this + y)</tt>
     */
    public DoubleDouble add(DoubleDouble y) {
        if (hi!=hi || y.hi!=y.hi) {
            return NaN;
        }
        double a, b, c;
        b = hi + y.hi;
        a = hi - b;
        c = ((hi - (a + b)) + (a + y.hi)) + (lo + y.lo);
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }
     
     
     
     
    /**
     * Returns a DoubleDouble whose value is <tt>(this + y)</tt>.
     *
     * @param y the addend
     * @return <tt>(this + y)</tt>
     */
    public DoubleDouble add(long y) {
        if (hi!=hi) {
            return NaN;
        }
        double yhi = y & 0xfffffffffffff800l;
        double ylo = y & 0x7ffl;
        double s = yhi + ylo;
        double err = ylo - (s - yhi);
        yhi = s;
        ylo = err;
        
        double a, b, c;
        b = hi + yhi;
        a = hi - b;
        c = ((hi - (a + b)) + (a + yhi)) + (lo + ylo);
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }
    
    

    
    /**
     * Returns a DoubleDouble whose value is <tt>(this + y)</tt>. This
     * method is slightly more accurate than the standard add() method,
     * but it also is slower. This method assures accuracy within one
     * ulp(), while the standard add() method can have an error of two
     * ulp().
     *
     * @param y the addend
     * @return <tt>(this + y)</tt>
     */
    public DoubleDouble addStrict(DoubleDouble y) {
        if (hi!=hi || y.hi!=y.hi) {
            return NaN;
        }
        double H, h, T, t, S, s, e, f;
        S = hi + y.hi;
        T = lo + y.lo;
        e = S - hi;
        f = T - lo;
        s = S - e;
        t = T - f;
        s = (y.hi - e) + (hi - s);
        t = (y.lo - f) + (lo - t);
        e = s + T;
        H = S + e;
        h = e + (S - H);
        e = t + h;

        double zhi = H + e;
        double zlo = e + (H - zhi);
        return new DoubleDouble(zhi, zlo);
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this + y)</tt>.
     *
     * @param y the addend
     * @return <tt>(this + y)</tt>
     */
    public DoubleDouble add(double y) {
        if (hi!=hi || y!=y) {
            return NaN;
        }

        double a, b, c;
        b = hi + y;
        a = hi - b;
        c = ((hi - (b + a)) + (y + a)) + lo;
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this + y)</tt>.
     *
     * @param y the addend
     * @return <tt>(this + y)</tt>
     */
    public DoubleDouble add(int y) {
        if (hi!=hi) {
            return NaN;
        }

        double a, b, c;
        b = hi + y;
        a = hi - b;
        c = ((hi - (b + a)) + (y + a)) + lo;
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this + 1)</tt>.
     *
     * @return <tt>(this + 1)</tt>
     */
    public DoubleDouble add1() {
        if (hi!=hi) {
            return NaN;
        }

        double a, b, c;
        b = hi + 1.0;
        a = hi - b;
        c = ((hi - (b + a)) + (1.0 + a)) + lo;
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this - y)</tt>.
     *
     * @param y the subtrahend
     * @return <tt>(this - y)</tt>
     */
    public DoubleDouble sub(DoubleDouble y) {
        if (hi!=hi || y.hi!=y.hi) {
            return NaN;
        }
        double a, b, c;
        b = hi - y.hi;
        a = hi - b;
        c = ((hi - (a + b)) + (a - y.hi)) + (lo - y.lo);
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this - y)</tt>.
     *
     * @param y the subtrahend
     * @return <tt>(this - y)</tt>
     */
    public DoubleDouble sub(long y) {
        if (hi!=hi) {
            return NaN;
        }
    
        double yhi = y & 0xfffffffffffff800l;
        double ylo = y & 0x7ffl;
        double s = yhi + ylo;
        double err = ylo - (s - yhi);
        yhi = s;
        ylo = err;
        
        double a, b, c;
        b = hi - yhi;
        a = hi - b;
        c = ((hi - (a + b)) + (a - yhi)) + (lo - ylo);
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }
    

    
    /**
     * Returns a DoubleDouble whose value is <tt>(this - y)</tt>. This
     * method is slightly more accurate than the standard sub() method,
     * but it also is slower. This method assures accuracy within one
     * ulp(), while the standard sub() method can have an error of two
     * ulp().
     *
     * @param y the addend
     * @return <tt>(this - y)</tt>
     */
    public DoubleDouble subStrict(DoubleDouble y) {
        if (hi!=hi || y.hi!=y.hi) {
            return NaN;
        }
        
        double H, h, T, t, S, s, e, f;
        S = hi - y.hi;
        T = lo - y.lo;
        e = S - hi;
        f = T - lo;
        s = S - e;
        t = T - f;
        s = (-y.hi - e) + (hi - s);
        t = (-y.lo - f) + (lo - t);
        e = s + T;
        H = S + e;
        h = e + (S - H);
        e = t + h;

        double zhi = H + e;
        double zlo = e + (H - zhi);
        return new DoubleDouble(zhi, zlo);
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this - y)</tt>.
     *
     * @param y the subtrahend
     * @return <tt>(this - y)</tt>
     */
    public DoubleDouble sub(double y) {
        if (hi!=hi || y!=y) {
            return NaN;
        }

        double a, b, c;
        b = hi - y;
        a = hi - b;
        c = ((hi - (b + a)) + (a - y)) + lo;
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this - y)</tt>.
     *
     * @param y the subtrahend
     * @return <tt>(this - y)</tt>
     */
    public DoubleDouble sub(int y) {
        if (hi!=hi) {
            return NaN;
        }

        double a, b, c;
        b = hi - y;
        a = hi - b;
        c = ((hi - (b + a)) + (a - y)) + lo;
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this - 1)</tt>.
     *
     * @return <tt>(this - 1)</tt>
     */
    public DoubleDouble sub1() {
        if (hi!=hi) {
            return NaN;
        }

        double a, b, c;
        b = hi - 1.0;
        a = hi - b;
        c = ((hi - (b + a)) + (a - 1.0)) + lo;
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>-this</tt>.
     *
     * @return <tt>-this</tt>
     */
    public DoubleDouble neg() {
        if (hi!=hi) {
            return NaN;
        }
        return new DoubleDouble(-hi, -lo);
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this * y)</tt>.
     *
     * @param y the multiplicand
     * @return <tt>(this * y)</tt>
     */
    public DoubleDouble mul(DoubleDouble y) {
        if (hi!=hi || y.hi!=y.hi) {
            return NaN;
        }
        
        double hx, tx, hy, ty, C, c;
        C = SPLIT * hi;
        hx = C - hi;
        c = SPLIT * y.hi;
        hx = C - hx;
        tx = hi - hx;
        hy = c - y.hi;
        C = hi * y.hi;
        hy = c - hy;
        ty = y.hi - hy;
        c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (hi * y.lo + lo * y.hi);
        double zhi = C + c;
        hx = C - zhi;
        double zlo = c + hx;
        return new DoubleDouble(zhi, zlo);
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this * y)</tt>.
     *
     * @param y the multiplicand
     * @return <tt>(this * y)</tt>
     */
    public DoubleDouble mul(long y) {
        if (hi!=hi) {
            return NaN;
        }
    
        double yhi = y & 0xfffffffffffff800l;
        double ylo = y & 0x7ffl;
        double s = yhi + ylo;
        double err = ylo - (s - yhi);
        yhi = s;
        ylo = err;
        
        double hx, tx, hy, ty, C, c;
        C = SPLIT * hi;
        hx = C - hi;
        c = SPLIT * yhi;
        hx = C - hx;
        tx = hi - hx;
        hy = c - yhi;
        C = hi * yhi;
        hy = c - hy;
        ty = yhi - hy;
        c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (hi * ylo + lo * yhi);
        double zhi = C + c;
        hx = C - zhi;
        double zlo = c + hx;
        return new DoubleDouble(zhi, zlo);
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this * y)</tt>.
     *
     * @param y the multiplicand
     * @return <tt>(this * y)</tt>
     */
    public DoubleDouble mul(double y) {
        if (hi!=hi || y!=y) {
            return NaN;
        }
        double hx, tx, hy, ty, C, c;
        C = SPLIT * hi;
        hx = C - hi;
        c = SPLIT * y;
        hx = C - hx;
        tx = hi - hx;
        hy = c - y;
        C = hi * y;
        hy = c - hy;
        ty = y - hy;
        c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (lo * y);
        double zhi = C + c;
        hx = C - zhi;
        double zlo = c + hx;
        return new DoubleDouble(zhi, zlo);
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this * y)</tt>.
     *
     * @param y the multiplicand
     * @return <tt>(this * y)</tt>
     */
    public DoubleDouble mul(int y) {
        if (hi!=hi) {
            return NaN;
        }
        double hx, tx, hy, ty, C, c;
        C = SPLIT * hi;
        hx = C - hi;
        c = SPLIT * y;
        hx = C - hx;
        tx = hi - hx;
        hy = c - y;
        C = hi * y;
        hy = c - hy;
        ty = y - hy;
        c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (lo * y);
        double zhi = C + c;
        hx = C - zhi;
        double zlo = c + hx;
        return new DoubleDouble(zhi, zlo);
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this * y)</tt>. This method
     * may ONLY be used when y is an exact power of two, e.g. 0.125, 0.25, 0.5,
     * 1, 2, 4, 8, and so on. A special much more efficient method for
     * multiplication is used in these cases.
     *
     * @param y the multiplicand
     * @return <tt>(this * y)</tt>
     */
    public DoubleDouble mulPowerOf2(double y) {
        if (hi!=hi) {
            return NaN;
        }
        return new DoubleDouble(hi * y, lo * y);
    }
    
    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this / y)</tt>.
     *
     * @param y the divisor
     * @return <tt>(this / y)</tt>
     */
    public DoubleDouble div(DoubleDouble y) {
        if (hi!=hi || y.hi!=y.hi) {
            return NaN;
        }
        double hc, tc, hy, ty, C, c, U, u;
        C = hi / y.hi;
        c = SPLIT * C;
        hc = c - C;
        u = SPLIT * y.hi;
        hc = c - hc;
        tc = C - hc;
        hy = u - y.hi;
        U = C * y.hi;
        hy = u - hy;
        ty = y.hi - hy;
        u = (((hc * hy - U) + hc * ty) + tc * hy) + tc * ty;
        c = ((((hi - U) - u) + lo) - C * y.lo) / y.hi;
        u = C + c;

        double zhi = u;
        double zlo = (C - u) + c;
        return new DoubleDouble(zhi, zlo);
    }
    
    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this / y)</tt>.
     *
     * @param y the divisor
     * @return <tt>(this / y)</tt>
     */
    public DoubleDouble div(long y) {
        if (hi!=hi) {
            return NaN;
        }
    
        double yhi = y & 0xfffffffffffff800l;
        double ylo = y & 0x7ffl;
        double s = yhi + ylo;
        double err = ylo - (s - yhi);
        yhi = s;
        ylo = err;
        
        double hc, tc, hy, ty, C, c, U, u;
        C = hi / yhi;
        c = SPLIT * C;
        hc = c - C;
        u = SPLIT * yhi;
        hc = c - hc;
        tc = C - hc;
        hy = u - yhi;
        U = C * yhi;
        hy = u - hy;
        ty = yhi - hy;
        u = (((hc * hy - U) + hc * ty) + tc * hy) + tc * ty;
        c = ((((hi - U) - u) + lo) - C * ylo) / yhi;
        u = C + c;

        double zhi = u;
        double zlo = (C - u) + c;
        return new DoubleDouble(zhi, zlo);
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this / y)</tt>.
     * This strict method is more precise than the normally used div()
     * method, but it also is slower. This strict method assures accuracy
     * within one ulp, the div() method sometimes may have one bit less
     * of accuracy, but in normal situations it works well.
     *
     * @param y the divisor
     * @return <tt>(this / y)</tt>
     */
    public DoubleDouble divStrict(DoubleDouble y) {
        if (hi!=hi || y.hi!=y.hi) {
            return NaN;
        }
        double a, b, c, d, e, f, g;
        f = hi / y.hi;
        a = 0x08000001 * y.hi;
        a += y.hi - a;
        b = y.hi - a;
        c = 0x08000001 * f;
        c += f - c;
        d = f - c;
        e = y.hi * f;
        c = (((a * c - e) + (a * d + b * c)) + b * d) + y.lo * f;
        b = lo - c;
        d = lo - b;
        a = hi - e;
        e = (hi - ((hi - a) + a)) + b;
        g = a + e;
        e += (a - g) + ((lo - (d + b)) + (d - c));
        a = g + e;
        b = a / y.hi;
        f += (e + (g - a)) / y.hi;
        a = f + b;
        return new DoubleDouble(a, b + (f - a));
    }

    
    
    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this / y)</tt>.
     *
     * @param y the divisor
     * @return <tt>(this / y)</tt>
     */
    public DoubleDouble div(double y) {
        if (hi!=hi || y!=y) {
            return NaN;
        }
        double hc, tc, hy, ty, C, c, U, u;
        C = hi / y;
        c = SPLIT * C;
        hc = c - C;
        u = SPLIT * y;
        hc = c - hc;
        tc = C - hc;
        hy = u - y;
        U = C * y;
        hy = u - hy;
        ty = y - hy;
        u = (((hc * hy - U) + hc * ty) + tc * hy) + tc * ty;
        c = (((hi - U) - u) + lo) / y;
        u = C + c;

        double zhi = u;
        double zlo = (C - u) + c;
        return new DoubleDouble(zhi, zlo);
    }

    
    
    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this / y)</tt>.
     *
     * @param y the divisor
     * @return <tt>(this / y)</tt>
     */
    public DoubleDouble div(int y) {
        if (hi!=hi) {
            return NaN;
        }
        double hc, tc, hy, ty, C, c, U, u;
        C = hi / y;
        c = SPLIT * C;
        hc = c - C;
        u = SPLIT * y;
        hc = c - hc;
        tc = C - hc;
        hy = u - y;
        U = C * y;
        hy = u - hy;
        ty = y - hy;
        u = (((hc * hy - U) + hc * ty) + tc * hy) + tc * ty;
        c = (((hi - U) - u) + lo) / y;
        u = C + c;

        double zhi = u;
        double zlo = (C - u) + c;
        return new DoubleDouble(zhi, zlo);
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this / y)</tt>. This method
     * may ONLY be used when y is an exact power of two, e.g. 0.125, 0.25, 0.5,
     * 1, 2, 4, 8, and so on. A special much more efficient method for
     * multiplication is used in these cases.
     *
     * @param y the divisor
     * @return <tt>(this * y)</tt>
     */
    public DoubleDouble divPowerOf2(double y) {
        if (hi!=hi) {
            return NaN;
        }
        y = 1.0 / y;
        return new DoubleDouble(hi * y, lo * y);
    }

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>1 / this</tt>.
     *
     * @return the reciprocal of this value
     */
    public DoubleDouble recip() {
        if (hi!=hi) {
            return NaN;
        }
        double hc, tc, hy, ty, C, c, U, u;
        C = 1.0 / hi;
        c = SPLIT * C;
        hc = c - C;
        u = SPLIT * hi;
        hc = c - hc;
        tc = C - hc;
        hy = u - hi;
        U = C * hi;
        hy = u - hy;
        ty = hi - hy;
        u = (((hc * hy - U) + hc * ty) + tc * hy) + tc * ty;
        c = ((((1.0 - U) - u)) - C * lo) / hi;

        double zhi = C + c;
        double zlo = (C - zhi) + c;
        return new DoubleDouble(zhi, zlo);
    }
    
    
    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(y / this)</tt>. Using this
     * method is much more efficient than creating a new DoubleDouble(y) and
     * calling div() on that value. This method allows a single double to be
     * divided by a DoubleDouble at DoubleDouble precision.
     *
     * @param y the divisor
     * @return <tt>(y / this)</tt>
     */
    public DoubleDouble recip(double y) {
        if (hi!=hi || y!=y) {
            return NaN;
        }
        double a, b, c, d, e, f;
        f = y / hi;
        a = 0x08000001 * hi;
        a += hi - a;
        b = hi - a;
        c = 0x08000001 * f;
        c += f - c;
        d = f - c;
        e = hi * f;
        b = ((y - e) - ((((a * c - e) + (a * d + b * c)) + b * d) + lo * f)) / hi;
        a = f + b;
        return new DoubleDouble(a, b + (f - a));
    }
    
    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(y / this)</tt>. Using this
     * method is much more efficient than creating a new DoubleDouble(y) and
     * calling div() on that value. This method allows a single long to be
     * divided by a DoubleDouble at DoubleDouble precision.
     *
     * @param y the divisor
     * @return <tt>(y / this)</tt>
     */
    public DoubleDouble recip(long y) {
        if (hi!=hi) {
            return NaN;
        }
    
        double yhi = y & 0xfffffffffffff800l;
        double ylo = y & 0x7ffl;
        double s = yhi + ylo;
        double err = ylo - (s - yhi);
        yhi = s;
        ylo = err;
        
        double hc, tc, hy, ty, C, c, U, u;
        C = yhi / hi;
        c = SPLIT * C;
        hc = c - C;
        u = SPLIT * hi;
        hc = c - hc;
        tc = C - hc;
        hy = u - hi;
        U = C * hi;
        hy = u - hy;
        ty = hi - hy;
        u = (((hc * hy - U) + hc * ty) + tc * hy) + tc * ty;
        c = ((((yhi - U) - u) + ylo) - C * lo) / hi;
        u = C + c;

        double zhi = u;
        double zlo = (C - u) + c;
        return new DoubleDouble(zhi, zlo);
    }
    
    
    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(y / this)</tt>. Using this
     * method is much more efficient than creating a new DoubleDouble(y) and
     * calling div() on that value. This method allows a single integer to be
     * divided by a DoubleDouble at DoubleDouble precision.
     *
     * @param y the divisor
     * @return <tt>(y / this)</tt>
     */
    public DoubleDouble recip(int y) {
        if (hi!=hi) {
            return NaN;
        }
        double a, b, c, d, e, f;
        f = y / hi;
        a = 0x08000001 * hi;
        a += hi - a;
        b = hi - a;
        c = 0x08000001 * f;
        c += f - c;
        d = f - c;
        e = hi * f;
        b = ((y - e) - ((((a * c - e) + (a * d + b * c)) + b * d) + lo * f)) / hi;
        a = f + b;
        return new DoubleDouble(a, b + (f - a));
    }
    
    


    
    /**
     * Computes the square of this value.
     *
     * @return the square of this value.
     */
    public DoubleDouble sqr() {
        if (hi!=hi) {
            return NaN;
        }
        double a, b, c;
        a = SPLIT * hi;
        a += hi - a;
        b = hi - a;
        c = hi * hi;
        b = ((((a * a - c) + a * b * 2) + b * b) + hi * lo * 2) + lo * lo;
        a = b + c;
        return new DoubleDouble(a, b + (c - a));
    }
    

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this * y + yy)</tt>.
     *
     * @param y the multiplicand
     * @param yy the addend
     * @return <tt>(this * y + yy)</tt>
     */
    public DoubleDouble muladd(DoubleDouble y, DoubleDouble yy) {
        if (hi!=hi || y.hi!=y.hi || yy.hi!=yy.hi) {
            return NaN;
        }
        
        double hx, tx, hy, ty, C, a, b, c;
        
        // Perform the multiplication.
        C = SPLIT * hi;
        hx = C - hi;
        c = SPLIT * y.hi;
        hx = C - hx;
        tx = hi - hx;
        hy = c - y.hi;
        C = hi * y.hi;
        hy = c - hy;
        ty = y.hi - hy;
        c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (hi * y.lo + lo * y.hi);
        double zhi = C + c;
        hx = C - zhi;
        double zlo = c + hx;
        
        // Perform the addition.
        b = zhi + yy.hi;
        a = zhi - b;
        c = ((zhi - (a + b)) + (a + yy.hi)) + (zlo + yy.lo);
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }
    

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this * y - yy)</tt>.
     *
     * @param y the multiplicand
     * @param yy the addend
     * @return <tt>(this * y - yy)</tt>
     */
    public DoubleDouble mulsub(DoubleDouble y, DoubleDouble yy) {
        if (hi!=hi || y.hi!=y.hi || yy.hi!=yy.hi) {
            return NaN;
        }
        
        double hx, tx, hy, ty, C, a, b, c;
        
        // Perform the multiplication.
        C = SPLIT * hi;
        hx = C - hi;
        c = SPLIT * y.hi;
        hx = C - hx;
        tx = hi - hx;
        hy = c - y.hi;
        C = hi * y.hi;
        hy = c - hy;
        ty = y.hi - hy;
        c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (hi * y.lo + lo * y.hi);
        double zhi = C + c;
        hx = C - zhi;
        double zlo = c + hx;
        
        // Perform the subtraction.
        b = zhi - yy.hi;
        a = zhi - b;
        c = ((zhi - (a + b)) + (a - yy.hi)) + (zlo - yy.lo);
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }
    

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this + y * yy)</tt>.
     *
     * @param y the first factor of the product
     * @param yy the second facvtor of the product
     * @return <tt>(this + y*yy)</tt>
     */
    public DoubleDouble addmul(DoubleDouble y, DoubleDouble yy) {
        if (hi!=hi || y.hi!=y.hi || yy.hi!=yy.hi) {
            return NaN;
        }
        
        double hx, tx, hy, ty, C, a, b, c;
        
        // Perform the multiplication and the negation.
        C = SPLIT * yy.hi;
        hx = C - yy.hi;
        c = SPLIT * y.hi;
        hx = C - hx;
        tx = yy.hi - hx;
        hy = c - y.hi;
        C = yy.hi * y.hi;
        hy = c - hy;
        ty = y.hi - hy;
        c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (yy.hi * y.lo + yy.lo * y.hi);
        double zhi = C + c;
        hx = C - zhi;
        double zlo = c + hx;
        
        // Perform the addition of this.
        b = zhi + hi;
        a = zhi - b;
        c = ((zhi - (a + b)) + (a + hi)) + (zlo + lo);
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }
    

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this - y * yy)</tt>.
     *
     * @param y the first factor of the product
     * @param yy the second facvtor of the product
     * @return <tt>(this - y*yy)</tt>
     */
    public DoubleDouble submul(DoubleDouble y, DoubleDouble yy) {
        if (hi!=hi || y.hi!=y.hi || yy.hi!=yy.hi) {
            return NaN;
        }
        
        double hx, tx, hy, ty, C, a, b, c;
        
        // Perform the multiplication and the negation.
        C = -SPLIT * yy.hi;
        hx = C + yy.hi;
        c = SPLIT * y.hi;
        hx = C - hx;
        tx = -yy.hi - hx;
        hy = c - y.hi;
        C = -yy.hi * y.hi;
        hy = c - hy;
        ty = y.hi - hy;
        c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) - (yy.hi * y.lo + yy.lo * y.hi);
        double zhi = C + c;
        hx = C - zhi;
        double zlo = c + hx;
        
        // Perform the addition of this.
        b = zhi + hi;
        a = zhi - b;
        c = ((zhi - (a + b)) + (a + hi)) + (zlo + lo);
        a = b + c;
        return new DoubleDouble(a, c + (b - a));
    }
    

    
    
    
    /**
     * Returns a DoubleDouble whose value is <tt>(this * y * yy)</tt>.
     *
     * @param y the multiplicand
     * @param yy the addend
     * @return <tt>(this * y * yy)</tt>
     */
    public DoubleDouble mulmul(DoubleDouble y, DoubleDouble yy) {
        if (hi!=hi || y.hi!=y.hi || yy.hi!=yy.hi) {
            return NaN;
        }
        
        double hx, tx, hy, ty, C, a, b, c;
        
        // Perform the first multiplication.
        C = SPLIT * hi;
        hx = C - hi;
        c = SPLIT * y.hi;
        hx = C - hx;
        tx = hi - hx;
        hy = c - y.hi;
        C = hi * y.hi;
        hy = c - hy;
        ty = y.hi - hy;
        c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (hi * y.lo + lo * y.hi);
        double zhi = C + c;
        hx = C - zhi;
        double zlo = c + hx;
        
        // Perform the second multiplication.
        C = SPLIT * zhi;
        hx = C - zhi;
        c = SPLIT * yy.hi;
        hx = C - hx;
        tx = zhi - hx;
        hy = c - yy.hi;
        C = zhi * yy.hi;
        hy = c - hy;
        ty = yy.hi - hy;
        c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (zhi * yy.lo + zlo * yy.hi);
        zhi = C + c;
        hx = C - zhi;
        zlo = c + hx;
        
        return new DoubleDouble(zhi, zlo);
    }

    
    
    
    /**
     * Returns the largest (closest to positive infinity) value that is not
     * greater than the argument and is equal to a mathematical integer. Special
     * cases:
     * <ul>
     * <li>If this value is NaN, returns NaN.
     * </ul>
     *
     * @return the largest (closest to positive infinity) value that is not
     * greater than the argument and is equal to a mathematical integer.
     */
    public DoubleDouble floor() {
        if (hi!=hi) {
            return NaN;
        }
        double fhi = StrictMath.floor(hi);
        if (fhi == hi) {
            // hi is already integral, floor the low word.
            double flo = StrictMath.floor(lo);
            double st = fhi + flo; flo = flo + (fhi - st); fhi = st;  
            return new DoubleDouble(fhi, flo);
        }
        
        return new DoubleDouble(fhi);
    }

    
    
    /**
     * Returns the smallest (closest to negative infinity) value that is not
     * less than the argument and is equal to a mathematical integer. Special
     * cases:
     * <ul>
     * <li>If this value is NaN, returns NaN.
     * </ul>
     *
     * @return the smallest (closest to negative infinity) value that is not
     * less than the argument and is equal to a mathematical integer.
     */
    public DoubleDouble ceil() {
        if (hi!=hi) {
            return NaN;
        }
        double fhi = StrictMath.ceil(hi);
        if (fhi == hi) {
            // hi is already integral, ceil the low word
            double flo = StrictMath.ceil(lo);
            double st = fhi + flo; flo = flo + (fhi - st); fhi = st;  
            return new DoubleDouble(fhi, flo);
        }
        return new DoubleDouble(fhi);
    }

    
    

    /**
     * Returns the integer which is largest in absolute value and not further
     * from zero than this value. Special cases:
     * <ul>
     * <li>If this value is NaN, returns NaN.
     * </ul>
     *
     * @return the integer which is largest in absolute value and not further
     * from zero than this value
     */
    public DoubleDouble trunc() {
        if (hi!=hi) {
            return NaN;
        }
        return (hi < 0.0 || (hi == 0.0 && lo < 0.0)) ? ceil() : floor();
    }

    
    
    
    /**
     * Returns an integer indicating the sign of this value.
     * <ul>
     * <li>if this value is > 0, returns 1
     * <li>if this value is < 0, returns -1 
     * <li>if this value is = 0, returns 0
     * <li>if this value is NaN, returns 0
     * </ul>
     *
     * @return an integer indicating the sign of this value
     */
    public int signum() {
        if (hi > 0.0 || (hi == 0.0 && lo > 0.0)) {
            return 1;
        }
        if (hi < 0.0 || (hi == 0.0 && lo < 0.0)) {
            return -1;
        }
        return 0;
    }

    
    
    
    /**
     * Rounds this value to the nearest integer. For positive integers, 0.5
     * is rounded upwards, e.g. 3.5 is rounded to 4. For negative integers,
     * 0.5 is rounded downwards (towards -infinity), e.g. -3.5 is rounded
     * to -4.
     * Special case:
     * <ul>
     * <li>If this value is NaN, returns NaN.
     * </ul>
     *
     * @return this value rounded to the nearest integer.
     */
    public DoubleDouble round() {
        if (hi!=hi) {
            return NaN;
        }
        
        return (hi > 0.0 || (hi == 0.0 && lo > 0.0)) ? add(0.5).floor() : sub(0.5).ceil();
    }

    
    
    
    
    /**
     * Returns the absolute value of this value. Special cases:
     * <ul>
     * <li>If this value is NaN, it is returned.
     * </ul>
     *
     * @return the absolute value of this value
     */
    public DoubleDouble abs() {
        if (hi!=hi) {
            return NaN;
        }
        if (hi < 0.0 || (hi == 0.0 && lo < 0.0)) {
            return new DoubleDouble(-hi, -lo);
        }
        return this;
    }

    
    


    
    /**
     * Computes the positive square root of this value. If the number is NaN or
     * negative, NaN is returned. 
     *
     * @return the positive square root of this number. If the argument is NaN
     * or less than zero, the result is NaN.
     */
    public DoubleDouble sqrt() {
        if (hi == 0.0 && lo == 0.0) {
            return this;
        }
        if (hi!=hi || hi < 0.0 || (hi == 0.0 && lo < 0.0)) {
            return NaN;
        }

        double a, b, c, d, e;
        d = 1 / StrictMath.sqrt(hi);
        e = hi * d;
        a = 0x08000001 * e;
        a += e - a;
        b = e - a;
        c = e * e;
        b = ((a * a - c) + a * b * 2) + b * b;
        a = hi - c;
        c = hi - a;
        //c = (a + ((((hi - (c + a)) + (c - c)) + lo) - b)) * d * 0.5;
        c = (a + ((((hi - (c + a))) + lo) - b)) * d * 0.5;
        a = e + c;
        b = e - a;
        return new DoubleDouble(a, (e - (b + a)) + (b + c));
    }
    
    
    
    
    
    /**
     * Computes the positive square root of this value. If the number is NaN or
     * negative, NaN is returned. This strict method is more precise than the
     * sqrt() method, but it also is slower. This strict method assures
     * accuracy within one ulp, the sqrt() method sometimes may have one bit
     * less of accuracy, but in normal situations it works well.
     *
     * @return the positive square root of this number. If the argument is NaN
     * or less than zero, the result is NaN.
     */
    
    public DoubleDouble sqrtStrict() {
        // Strategy:  Use Karp's trick:  if x is an approximation
        // to sqrt(a), then
        //
        // sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)
        //
        // The approximation is accurate to twice the accuracy of x.
        // Also, the multiplication (a*x) and [-]*x can be done with
        // only half the precision.
        if (hi == 0.0 && lo == 0.0) {
            return this;
        }
        if (hi!=hi || hi < 0.0 || (hi == 0.0 && lo < 0.0)) {
            return NaN;
        }

        double a, b, c, d, e, f, g, h;
        g = 1 / StrictMath.sqrt(hi);
        h = hi * g;
        g *= 0.5;
        a = 0x08000001 * h;
        a += h - a;
        b = h - a;
        c = h * h;
        b = ((a * a - c) + a * b * 2) + b * b;
        a = lo - b;
        f = lo - a;
        e = hi - c;
        d = hi - e;
        d = ((hi - (d + e)) + (d - c)) + a;
        c = e + d;
        b = (d + (e - c)) + ((lo - (f + a)) + (f - b));
        a = c + b;
        b += (c - a);
        c = 0x08000001 * a;
        c += a - c;
        d = a - c;
        e = 0x08000001 * g;
        e += g - e;
        f = g - e;
        a = a * g;
        e = ((c * e - a) + (c * f + d * e)) + d * f;
        e += b * g;
        b = a + e;
        e += a - b;
        f = b + h;
        c = b - f;
        return new DoubleDouble(f, e + ((b - (f + c)) + (h + c)));
    }
    
    // The function below is commented away, but kept here
    // for documenting purposes. This is Karp's trick made
    // explicit in working Java code. The code above is
    // used with everything written out such that no method
    // calls and instantiation overhead occurs.
    // The function below is equally well as the one above.
    /*
    public DoubleDouble sqrt() {

        if (hi == 0.0 && lo == 0.0) {
            return this;
        }
        if (hi!=hi || hi < 0.0 || (hi == 0.0 && lo < 0.0)) {
            return NaN;
        }

        double x = 1.0 / StrictMath.sqrt(hi);
        double ax = hi * x;

        DoubleDouble axdd = new DoubleDouble(ax);
        DoubleDouble diffSq = this.sub(axdd.sqr());
        double d2 = diffSq.hi * (x * 0.5);

        return axdd.add(new DoubleDouble(d2));
    }
    */

    
    
    
    /**
     * Computes the value of this number raised to an integral power. Follows
     * semantics of Java Math.pow as closely as possible.
     *
     * @param exp the integer exponent
     * @return x raised to the integral power exp
     */
    public DoubleDouble pow(int exp) {
        if (hi!=hi) {
            return NaN;
        }
        if (exp == 0) {
            return ONE;
        }

        // Variable r = this;
        double rhi = hi;
        double rlo = lo;
        
        // Variable s = 1.0;
        double shi = 1.0;
        double slo = 0.0;
        
        int n = (exp > 0) ? exp : -exp;
        if (n > 1) {
            /* Use binary exponentiation */
            while (n > 0) {
                if ((n & 1) == 1) {
                    // s = s*r
                    double hx, tx, hy, ty, C, c;
                    C = SPLIT * shi;
                    hx = C - shi;
                    c = SPLIT * rhi;
                    hx = C - hx;
                    tx = shi - hx;
                    hy = c - rhi;
                    C = shi * rhi;
                    hy = c - hy;
                    ty = rhi - hy;
                    c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (shi * rlo + slo * rhi);
                    shi = C + c;
                    hx = C - shi;
                    slo = c + hx;
                }
                n >>= 1;
                if (n > 0) {
                    // r = r^2
                    double a, b, c;
                    a = SPLIT * rhi;
                    a += rhi - a;
                    b = rhi - a;
                    c = rhi * rhi;
                    b = ((((a * a - c) + a * b * 2) + b * b) + rhi * rlo * 2) + rlo * rlo;
                    a = b + c;
                    rhi = a;
                    rlo = b + (c - a);
                }
            }
        } 
        else {
            shi = rhi;
            slo = rlo;
        }

        /* Compute the reciprocal if n is negative. */
        if (exp < 0) {
            // s = 1/s
            double a, b, c, d, e, f;
            f = 1.0 / shi;
            a = SPLIT * shi;
            a += shi - a;
            b = shi - a;
            c = SPLIT * f;
            c += f - c;
            d = f - c;
            e = shi * f;
            b = ((1.0 - e) - ((((a * c - e) + (a * d + b * c)) + b * d) + slo * f)) / shi;
            a = f + b;
            shi = a;
            slo = b + (f - a);
        }
        
        return new DoubleDouble(shi, slo);
    }
    
    
    
    
    /**
     * Computes y-th root of a DoubleDouble. The number y
     * must be positive.
     *
     * @param y the integer telling which root is taken
     * @return this, raised to the power 1/y.
     */
    public DoubleDouble root(int y) {
        if (y <=0 || hi != hi) {
            return NaN;
        }
        boolean isNeg = hi < 0.0 || (hi == 0.0 && lo < 0.0);

        if (isNeg && ((y & 1) == 0)) {
            return NaN;
        }
        
        if (hi == 0 && lo == 0) {
            return this;
        }
        
        switch (y) {
            case 1:
                return this;
            case 2:
                return sqrt();
            default:
                // fall through.
        }
        
        double a, b, c, d, e, f, g, h, i, j, k, l, m;
        int z;

        if (isNeg) {
            b = -hi;
            c = -lo;
        }
        else {
            b = hi;
            c = lo;
        }
        
        a = StrictMath.exp(StrictMath.log(b) / (-y));
        z = y;
        k = a;
        l = 0;
        g = 1;
        h = 0;
        while (z > 0) {
            if ((z & 1) > 0) {
                d = SPLIT * g;
                d += g - d;
                e = g - d;
                f = SPLIT * k;
                f += k - f;
                i = k - f;
                j = g * k;
                h = (((d * f - j) + (d * i + e * f)) + e * i) + (h * k + g * l);
                g = j + h;
                h += j - g;
            }
            f = SPLIT * k;
            f = f + (k - f);
            i = k - f;
            j = k * k;
            i = ((f * f - j) + f * i * 2) + i * i;
            i += k * l * 2;
            i += l * l;
            k = i + j;
            l = i + (j - k);
            z >>= 1;
        }

        l = SPLIT * b;
        l += b - l;
        m = b - l;
        d = SPLIT * g;
        d += g - d;
        e = g - d;
        f = b * g;
        d = (((l * d - f) + (l * e + m * d)) + m * e) + (c * g + b * h);
        e = 1 - f;
        l = e - d;
        m = (e - l) - d;
        d = SPLIT * l;
        d += l - d;
        e = l - d;
        f = SPLIT * a;
        f += a - f;
        g = a - f;
        l *= a;
        m *= a;
        m += (((d * f - l) + (d * g + e * f)) + e * g);
        d = l / y;
        e = SPLIT * d;
        e += d - e;
        f = d - e;
        g = SPLIT * y;
        g += y - g;
        h = y - g;
        i = d * y;
        j = l - i;
        k = l - j;
        m = (j + ((((l - (k + j)) + (k - i)) + m) - (((e * g - i) + (e * h + f * g)) + f * h))) / y;
        e = d + a;
        l = d - e;
        m += (d - (e + l)) + (a + l);
        if (isNeg) {
            e = -e;
            m = -m;
        }
        i = 1 / e;
        l = SPLIT * e;
        l += e - l;
        d = e - l;
        f = SPLIT * i;
        f += i - f;
        g = i - f;
        h = e * i;
        m = ((1 - h) - ((((l * f - h) + (l * g + d * f)) + d * g) + m * i)) / e;
        l = i + m;
        return new DoubleDouble(l, m + (i - l));
    }
    
    
    

    
    // Devil's values:
    // 0.693147180559945309417232121458174
    // 1.03972077083991796412584818218727
    // 1.03972077083991796312584818218727
    public DoubleDouble exp() {
        if (hi!=hi) {
            return NaN;
        }
        if (hi > 691.067739) {
            return new DoubleDouble(Double.POSITIVE_INFINITY);
        }

        double a, b, c, d, e, f, g = 0.5, h = 0, i, j, k, l, m, n, o, p, q = 2, r = 1;
        int s;

        a = SPLIT * hi;
        a += hi - a;
        b = a - hi;
        c = hi * 1.4426950408889634;
        b = (((a * 1.4426950514316559 - c) - (b * 1.4426950514316559 + a * 1.0542692496784412E-8)) + b * 1.0542692496784412E-8)
                + (lo * 1.4426950408889634 + hi * 2.0355273740931033E-17);
        s = (int) StrictMath.round(c);
        if (c == s) {
            s += (int) StrictMath.round(b);
        } 
        else if (StrictMath.abs(s - c) == 0.5 && b < 0.0) {
            s--;
        }
        e = 0.6931471805599453 * s;
        c = ((s * 0.6931471824645996 - e) - (s * 1.904654323148236E-9)) + 2.3190468138462996E-17 * s;
        b = lo - c;
        d = lo - b;
        e = hi - e;
        a = e + b;
        b = ((lo - (d + b)) + (d - c)) + (b + (e - a));
        e = a + 1;
        c = a - e;
        d = ((a - (e + c)) + (1 + c)) + b;
        c = e + d;
        d += e - c;
        e = SPLIT * a;
        e += a - e;
        f = a - e;
        i = a * a;
        f = ((e * e - i) + e * f * 2) + f * f;
        f += a * b * 2;
        f += b * b;
        e = f + i;
        f += i - e;
        i = e * g;
        j = f * g;
        do {
            k = d + j;
            l = d - k;
            m = c + i;
            n = c - m;
            n = ((c - (n + m)) + (n + i)) + k;
            o = m + n;
            d = (n + (m - o)) + ((d - (l + k)) + (l + j));
            c = o + d;
            d += o - c;
            k = SPLIT * e;
            k += e - k;
            l = e - k;
            m = SPLIT * a;
            m += a - m;
            n = a - m;
            o = e * a;
            f = (((k * m - o) + (k * n + l * m)) + l * n) + (f * a + e * b);
            e = o + f;
            f += o - e;
            n = g / ++q;
            k = SPLIT * n;
            k += n - k;
            l = n - k;
            m = n * q;
            o = g - m;
            p = g - o;
            h = (o + ((((g - (p + o)) + (p - m)) + h) - (((k * q - m) + l * q)))) / q;
            g = n;
            i = SPLIT * e;
            i += e - i;
            k = e - i;
            j = SPLIT * g;
            j += g - j;
            l = g - j;
            m = e * g;
            j = (((i * j - m) + (i * l + k * j)) + k * l) + (f * g + e * h);
            i = m + j;
            j += m - i;
        } 
        while (i > 1e-40 || i < -1e-40);

        if (s < 0) {
            s = -s;
            a = 0.5;
        } 
        else {
            a = 2;
        }

        while (s > 0) {
            if ((s & 1) > 0) {
                r *= a;
            }
            a *= a;
            s >>= 1;
        }
        a = d + j;
        b = d - a;
        e = c + i;
        f = c - e;
        f = ((c - (f + e)) + (f + i)) + a;
        c = e + f;
        d = (f + (e - c)) + ((d - (b + a)) + (b + j));
        return new DoubleDouble(c * r, d * r);
    }
    
    
    
    
    

    public DoubleDouble log() {
        if (hi!=hi || hi <= 0.0) {
            return NaN;
        }

        double a, b, c, d, e, f, g = 0.5, h = 0, i, j, k, l, m, n, o, p, q = 2, r = 1, s;
        int t;

        s = StrictMath.log(hi);

        a = SPLIT * s;
        a += s + a;
        b = s - a;
        c = s * -1.4426950408889634;
        b = (((a * -1.4426950514316559 - c) + (a * 1.0542692496784412E-8 - b * 1.4426950514316559)) + b * 1.0542692496784412E-8) - (s * 2.0355273740931033E-17);
        t = (int) StrictMath.round(c);
        if (a == t) {
            t += (int) StrictMath.round(b);
        } 
        else if (StrictMath.abs(t - a) == 0.5 && b < 0.0) {
            t--;
        }
        e = 0.6931471805599453 * t;
        c = ((t * 0.6931471824645996 - e) - (t * 1.904654323148236E-9)) + 2.3190468138462996E-17 * t;
        e += s;
        a = e + c;
        b = (a - e) - c;
        e = 1 - a;
        d = ((1 - e) - a) + b;
        c = e + d;
        d += e - c;
        e = SPLIT * -a;
        e -= a + e;
        f = a + e;
        i = a * a;
        f = ((e * e - i) - e * f * 2) + f * f;
        f += -a * b * 2;
        a = -a;
        f += b * b;
        e = f + i;
        f += i - e;
        l = SPLIT * e;
        l += e - l;
        //k = e - l;
        i = e * g;
        j = f * g;
        do {
            k = d + j;
            l = d - k;
            m = c + i;
            n = c - m;
            n = ((c - (n + m)) + (n + i)) + k;
            o = m + n;
            d = (n + (m - o)) + ((d - (l + k)) + (l + j));
            c = o + d;
            d += o - c;
            k = SPLIT * e;
            k += e - k;
            l = e - k;
            m = SPLIT * a;
            m += a - m;
            n = a - m;
            o = e * a;
            f = (((k * m - o) + (k * n + l * m)) + l * n) + (f * a + e * b);
            e = o + f;
            f += o - e;
            n = g / ++q;
            k = SPLIT * n;
            k += n - k;
            l = n - k;
            m = n * q;
            o = g - m;
            p = g - o;
            h = (o + ((((g - (p + o)) + (p - m)) + h) - (((k * q - m) + l * q)))) / q;
            g = n;
            i = SPLIT * e;
            i += e - i;
            k = e - i;
            j = SPLIT * g;
            j += g - j;
            l = g - j;
            m = e * g;
            j = (((i * j - m) + (i * l + k * j)) + k * l) + (f * g + e * h);
            i = m + j;
            j += m - i;
        } 
        while (i > 1e-40 || i < -1e-40);

        if (t < 0) {
            t = -t;
            k = 0.5;
        } 
        else {
            k = 2;
        }

        while (t > 0) {
            if ((t & 1) > 0) {
                r *= k;
            }
            k *= k;
            t >>= 1;
        }
        a = d + j;
        b = d - a;
        e = c + i;
        f = c - e;
        f = ((c - (f + e)) + (f + i)) + a;
        g = e + f;
        h = ((f + (e - g)) + ((d - (b + a)) + (b + j))) * r;
        g *= r;
        a = SPLIT * hi;
        a += hi - a;
        c = hi - a;
        b = SPLIT * g;
        b += g - b;
        d = g - b;
        e = hi * g;
        b = (((a * b - e) + (a * d + c * b)) + c * d) + (lo * g + hi * h);
        a = --e + b;
        b += e - a;
        c = a + s;
        d = a - c;
        b += ((a - (c + d)) + (s + d));
        a = c + b;
        return new DoubleDouble(a, b + (c - a));
    }

    
    
    
    /*------------------------------------------------------------
     *   Conversion Functions
     *------------------------------------------------------------
     */
    /**
     * Converts this value to the nearest double-precision number.
     *
     * @return the nearest double-precision number to this value
     */
    public double doubleValue() {
        return hi + lo;
    }

    
    
    
    /*------------------------------------------------------------
     *   Conversion Functions
     *------------------------------------------------------------
     */
    /**
     * Converts this value to the absolute value of the
     * nearest double-precision number.
     *
     * @return the absolute value of the nearest double-precision 
     * number to this value
     */
    public double dabs() {
        double a = hi + lo;
        return (a >= 0) ? a : -a;
    }

    
    
    
    /**
     * Converts this value to the nearest integer.
     *
     * @return the nearest integer to this value
     */
    public int intValue() {
        double fhi;
        if (hi > 0.0 || (hi == 0.0 && lo >= 0.0)) {
            // Non-negative number, use the floor function.
            fhi = StrictMath.floor(hi);
            if (fhi == hi) {
                fhi += StrictMath.floor(lo);  
            }
        }
        else {
            // Negative number, use the ceil function.
            fhi = StrictMath.ceil(hi);
            if (fhi == hi) {
                fhi += StrictMath.ceil(lo);  
            }
        }
        return (fhi < MININT.hi) ? Integer.MIN_VALUE : (fhi > MAXINT.hi ? Integer.MAX_VALUE : (int)fhi);
    }

    
    
    
    /**
     * Converts this value to the nearest integer.
     *
     * @return the nearest integer to this value
     */
    public long longValue() {
        double fhi, flo = 0.0;
        if (hi > 0.0 || (hi == 0.0 && lo >= 0.0)) {
            // Non-negative number, use the floor function.
            fhi = StrictMath.floor(hi);
            if (fhi == hi) {
                flo = StrictMath.floor(lo);
                double st = fhi + flo; flo = flo + (fhi - st); fhi = st;
            }
        }
        else {
            // Negative number, use the ceil function.
            fhi = StrictMath.ceil(hi);
            if (fhi == hi) {
                flo = StrictMath.ceil(lo);
                double st = fhi + flo; flo = flo + (fhi - st); fhi = st;  
            }
        }
        
        if (fhi > MAXLONG.hi || (fhi == MAXLONG.hi && flo >= MAXLONG.lo)) return Long.MAX_VALUE;
        if (fhi < MINLONG.hi || (fhi == MINLONG.hi && flo <= MINLONG.lo)) return Long.MIN_VALUE;
        
        return ((long)fhi) + ((long)flo);
    }

    
    
    
    /*------------------------------------------------------------
     *   Predicates
     *------------------------------------------------------------
     */
    
    /**
     * Tests whether this value is equal to 0.
     *
     * @return true if this value is equal to 0
     */
    public boolean isZero() {
        return hi == 0.0 && lo == 0.0;
    }

    
    
    
    /**
     * Tests whether this value is equal to 1.
     *
     * @return true if this value is equal to 1
     */
    public boolean isOne() {
        return hi == 1.0 && lo == 0.0;
    }

    
    
    
    /**
     * Tests whether this value is less than 0.
     *
     * @return true if this value is less than 0
     */
    public boolean isNegative() {
        return hi < 0.0 || (hi == 0.0 && lo < 0.0);
    }

    
    
    
    /**
     * Tests whether this value is greater than 0.
     *
     * @return true if this value is greater than 0
     */
    public boolean isPositive() {
        return hi > 0.0 || (hi == 0.0 && lo > 0.0);
    }

    
    
    
    /**
     * Tests whether this value is NaN.
     *
     * @return true if this value is NaN
     */
    public boolean isNaN() {
        return hi!=hi;
    }

    
    
    
    /**
     * Tests whether this value is close to another value.
     *
     * @param b a DoubleDouble value, to be compared with this
     * @param eps a specification of the relatieve precision
     * @return true if this value is close to the other value with relative
     * precision, given by eps.
     */
    public boolean isNear(DoubleDouble b, double eps) {
        if (hi!=hi || b.hi!=b.hi || eps!=eps) {
            return false;
        }
        double a = abs().doubleValue();
        double diff = sub(b).abs().doubleValue();
        return diff < a * eps;
    }

    
    
    
    /**
     * Tests whether this value is equal to another <tt>DoubleDouble</tt> value.
     *
     * @param y a DoubleDouble value
     * @return true if this value = y
     */
    public boolean equals(DoubleDouble y) {
        return hi == y.hi && lo == y.lo;
    }

    
    
    
    /**
     * Tests whether this value is equal to a given <tt>double</tt> value.
     *
     * @param y a double value
     * @return true if this value = y
     */
    public boolean equals(double y) {
        return hi == y && lo == 0;
    }

    
    
    
    /**
     * Tests whether this value is equal to a given <tt>long</tt> value.
     *
     * @param y a long value
     * @return true if this value = y
     */
    public boolean equals(long y) {
        double yhi = y & 0xfffffffffffff800l;
        double ylo = y & 0x7ffl;
        double s = yhi + ylo;
        double err = ylo - (s - yhi);
        yhi = s;
        ylo = err;
        return hi == yhi && lo == ylo;
    }

    
    
    
    /**
     * Tests whether this value is equal to a given <tt>int</tt> value.
     *
     * @param y an int value
     * @return true if this value = y
     */
    public boolean equals(int y) {
        return hi == y && lo == 0;
    }

    
    
    
    /**
     * Tests whether this value is greater than another <tt>DoubleDouble</tt>
     * value.
     *
     * @param y a DoubleDouble value
     * @return true if this value > y
     */
    public boolean gt(DoubleDouble y) {
        return (hi > y.hi) || (hi == y.hi && lo > y.lo);
    }

    
    
    
    /**
     * Tests whether this value is greater than a given <tt>double</tt>
     * value.
     *
     * @param y a double value
     * @return true if this value > y
     */
    public boolean gt(double y) {
        return (hi > y) || (hi == y && lo > 0.0);
    }

    
    
    
    /**
     * Tests whether this value is greater than a given
     * <tt>long</tt> value.
     *
     * @param y a long value
     * @return true if this value >= y
     */
    public boolean gt(long y) {
        double yhi = y & 0xfffffffffffff800l;
        double ylo = y & 0x7ffl;
        double s = yhi + ylo;
        double err = ylo - (s - yhi);
        yhi = s;
        ylo = err;
        return (hi > yhi) || (hi == yhi && lo > ylo);
    }

    
    
    
    /**
     * Tests whether this value is greater than a given <tt>int</tt>
     * value.
     *
     * @param y an int value
     * @return true if this value > y
     */
    public boolean gt(int y) {
        return (hi > y) || (hi == y && lo > 0.0);
    }

    
    
    
    /**
     * Tests whether this value is greater than or equals to another
     * <tt>DoubleDouble</tt> value.
     *
     * @param y a DoubleDouble value
     * @return true if this value >= y
     */
    public boolean ge(DoubleDouble y) {
        return (hi > y.hi) || (hi == y.hi && lo >= y.lo);
    }

    
    
    
    /**
     * Tests whether this value is greater than or equals to a given
     * <tt>double</tt> value.
     *
     * @param y a double value
     * @return true if this value >= y
     */
    public boolean ge(double y) {
        return (hi > y) || (hi == y && lo >= 0.0);
    }

    
    
    
    /**
     * Tests whether this value is greater than or equals to a given
     * <tt>long</tt> value.
     *
     * @param y a long value
     * @return true if this value >= y
     */
    public boolean ge(long y) {
        double yhi = y & 0xfffffffffffff800l;
        double ylo = y & 0x7ffl;
        double s = yhi + ylo;
        double err = ylo - (s - yhi);
        yhi = s;
        ylo = err;
        return (hi > yhi) || (hi == yhi && lo >= ylo);
    }

    
    
    
    /**
     * Tests whether this value is greater than or equals to a given
     * <tt>int</tt> value.
     *
     * @param y an int value
     * @return true if this value >= y
     */
    public boolean ge(int y) {
        return (hi > y) || (hi == y && lo >= 0.0);
    }

    
    
    
    /**
     * Tests whether this value is less than another <tt>DoubleDouble</tt>
     * value.
     *
     * @param y a DoubleDouble value
     * @return true if this value < y
     */
    public boolean lt(DoubleDouble y) {
        return (hi < y.hi) || (hi == y.hi && lo < y.lo);
    }

    
    
    
    /**
     * Tests whether this value is less than a given <tt>double</tt>
     * value.
     *
     * @param y a double value
     * @return true if this value < y
     */
    public boolean lt(double y) {
        return (hi < y) || (hi == y && lo < 0.0);
    }

    
    
    
    /**
     * Tests whether this value is less than a given
     * <tt>long</tt> value.
     *
     * @param y a long value
     * @return true if this value <= y
     */
    public boolean lt(long y) { 
        double yhi = y & 0xfffffffffffff800l;
        double ylo = y & 0x7ffl;
        double s = yhi + ylo;
        double err = ylo - (s - yhi);
        yhi = s;
        ylo = err;
        return (hi < yhi) || (hi == yhi && lo < ylo);
    }

    
    
    
    /**
     * Tests whether this value is less than a given <tt>int</tt>
     * value.
     *
     * @param y an int value
     * @return true if this value < y
     */
    public boolean lt(int y) {
        return (hi < y) || (hi == y && lo < 0.0);
    }

    
    
    
    /**
     * Tests whether this value is less than or equal to another
     * <tt>DoubleDouble</tt> value.
     *
     * @param y a DoubleDouble value
     * @return true if this value <= y
     */
    public boolean le(DoubleDouble y) {
        return (hi < y.hi) || (hi == y.hi && lo <= y.lo);
    }

    
    
    
    /**
     * Tests whether this value is less than or equal to a given
     * <tt>double</tt> value.
     *
     * @param y a double value
     * @return true if this value <= y
     */
    public boolean le(double y) {
        return (hi < y) || (hi == y && lo <= 0.0);
    }

    
    
    
    /**
     * Tests whether this value is less than or equal to a given
     * <tt>long</tt> value.
     *
     * @param y a long value
     * @return true if this value <= y
     */
    public boolean le(long y) { 
        double yhi = y & 0xfffffffffffff800l;
        double ylo = y & 0x7ffl;
        double s = yhi + ylo;
        double err = ylo - (s - yhi);
        yhi = s;
        ylo = err;
        return (hi < yhi) || (hi == yhi && lo <= ylo);
    }

    
    
    
    /**
     * Tests whether this value is less than or equal to a given
     * <tt>int</tt> value.
     *
     * @param y an int value
     * @return true if this value <= y
     */
    public boolean le(int y) {
        return (hi < y) || (hi == y && lo <= 0.0);
    }

    
    
    
    /**
     * Compares two DoubleDouble objects numerically.
     *
     * @return -1,0 or 1 depending on whether this value is less than, equal to
     * or greater than the value of <tt>o</tt>
     */
    @Override
    public int compareTo(Object o) {
        if (o instanceof DoubleDouble) {
            DoubleDouble other = (DoubleDouble)o;
            if (hi < other.hi || (hi == other.hi && lo < other.lo)) {return -1;}
            if (hi > other.hi || (hi == other.hi && lo > other.lo)) {return 1;}
            return 0;
        }
        
        if (o instanceof Double) {
            double other = (Double)o;
            if (hi < other || (hi == other && lo < 0)) {return -1;}
            if (hi > other || (hi == other && lo > 0)) {return 1;}
            return 0;
        }
        
        if (o instanceof Long) {
            DoubleDouble other = new DoubleDouble((Long)o);
            if (hi < other.hi || (hi == other.hi && lo < other.lo)) {return -1;}
            if (hi > other.hi || (hi == other.hi && lo > other.lo)) {return 1;}
            return 0;
        }
        
        if (o instanceof Integer) {
            double other = (Integer)o;
            if (hi < other || (hi == other && lo < 0)) {return -1;}
            if (hi > other || (hi == other && lo > 0)) {return 1;}
            return 0;
        }
        
        if (o instanceof Float) {
            double other = (Float)o;
            if (hi < other || (hi == other && lo < 0)) {return -1;}
            if (hi > other || (hi == other && lo > 0)) {return 1;}
            return 0;
        }
        
        if (o instanceof String) {
            DoubleDouble other = new DoubleDouble((String)o);
            if (hi < other.hi || (hi == other.hi && lo < other.lo)) {return -1;}
            if (hi > other.hi || (hi == other.hi && lo > other.lo)) {return 1;}
            return 0;
        }
        
        if (o instanceof Short) {
            double other = (Short)o;
            if (hi < other || (hi == other && lo < 0)) {return -1;}
            if (hi > other || (hi == other && lo > 0)) {return 1;}
            return 0;
        }
        
        if (o instanceof Byte) {
            double other = (Byte)o;
            if (hi < other || (hi == other && lo < 0)) {return -1;}
            if (hi > other || (hi == other && lo > 0)) {return 1;}
            return 0;
        }
        
        if (o instanceof Character) {
            double other = (Character)o;
            if (hi < other || (hi == other && lo < 0)) {return -1;}
            if (hi > other || (hi == other && lo > 0)) {return 1;}
            return 0;
        }
        
        throw new RuntimeException("Cannot compare class " + o.getClass().getName() + " to DoubleDouble value.");
    }

    
    
    
    
    /*------------------------------------------------------------
     *   Output
     *------------------------------------------------------------
     */
    
    /**
     * Dumps the components of this number to a string.
     *
     * @return a string showing the components of the number
     */
    public String dump() {
        return "DD<" + hi + ", " + lo + ">";
    }
    
        
    
    private static final char[] BASE_36_TABLE = { //
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', //
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', //
        'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', //
        'U', 'V', 'W', 'X', 'Y', 'Z'};
    
    private static final char[] ZEROES = { //
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', //
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', //
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', //
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', //
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', //
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', //
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', //
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', //
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', //
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', //
        '0', '0', '0', '0', '0'};
    
    
    private static final double[] DIGITS_PER_BIT = new double[37];
    static {
        double log2 = StrictMath.log(2);
        for (int base = 2; base <= 36; base++) {
            DIGITS_PER_BIT[base] = log2 / StrictMath.log(base);
        }
    }
    
    
    /**
     * Prints a DoubleDouble in the decimal number system to a String.
     * Normal notation is used for numbers in the range 0.001 to 1e10,
     * numbers outside of this range are printed in scientific notation.
     * @return The String, representing the number.
     */
    @Override
    public String toString() {
        return toString(10, 0, false);
    }

    
    
    /**
     * Prints this DoubleDouble to a String. This method provides a lot
     * of flexibility on how to print the number.
     * @param radix The radix of the number system in which this DoubleDouble
     * is printed. The exponent (if scientific notation is used) is always
     * printed in the decimal system. E.g. in the hexadecimal system, the
     * number A.BCDEFe20 stands for A.BCDEF*power, where power is 16^20 in
     * the decimal system.
     * @param maxPrecision The maximum number of digits in the given radix to
     * be printed. If the value 0 is supplied, then the number of digits is
     * such that the entire precision of the number is used (for radix 10 this
     * means that 31 digits are printed, for radix 16, this means that 26 
     * digits are printed). If a value is supplied, which uses more than the
     * precision of a DoubleDouble, then the number of digits is limited as if
     * the value 0 were aupplied.
     * @param sci A boolean, telling whether normal notation or scientific
     * notation with exponents needs to be used.
     * @return 
     */
    public String toString(int radix, int maxPrecision, boolean sci) {
        if (hi!=hi) {
            return "NaN";
        }
        
        // Get the precision. (The number of significant digits required
        // for an accurate representation of this number). 
        int precision = (int)(105 * DIGITS_PER_BIT[radix]);
        if (maxPrecision > 0 && maxPrecision < precision) {
            precision = maxPrecision;
        }

        // Get the raw digit representation.
        char[] chars = new char[precision + 1];
        int exp = to_digits(chars, precision, radix) + 1;

        // Get some properties.
        int left = StrictMath.max(0, -exp);
        int right = StrictMath.max(0, exp);
        if (chars[precision - 1] == 0) {
            precision--;
        }
        sci = sci || (exp<-3 || exp>10 || exp>precision);

        // Allocate exactly the right size string.
        StringBuilder out = new StringBuilder(precision + (sci ? 3 : left) + (exp > 0 ? 1 : 2));

        // Build the string.
        if (hi < 0) {
            out.append('-');
        }
        if (sci) {
            out.append(chars, 0, 1);
            out.append('.');
            out.append(chars, 1, precision - 1);
            out.append('e');
            out.append(exp - 1);
        } else {
            if (exp <= 0) {
                out.append('0');
            }
            if (right > 0) {
                out.append(chars, 0, right);
            }
            out.append('.');
            if (left > 0) {
                if (DoubleDouble.ZEROES.length < left) {
                    System.err.println(left);
                } else {
                    out.append(DoubleDouble.ZEROES, 0, left);
                }
            }
            out.append(chars, right, precision - right);
        }

        return out.toString();
    }

    
    
    
    private int to_digits(char[] s, int precision, int base) {
        int halfBase = (base + 1) >> 1;

        if (hi == 0.0 && lo == 0) {
            Arrays.fill(s, 0, precision, '0');
            return 0;
        }

        // First determine the (approximate) exponent.
        DoubleDouble temp = this.abs();
        int exp = (int) StrictMath.floor(StrictMath.log(temp.hi) / StrictMath.log(base));

        DoubleDouble p = new DoubleDouble(base);
        if (exp < -300) {
            temp = temp.mul(p.pow(150));
            p = p.pow(-exp - 150);
            temp = temp.mul(p);
        } 
        else {
            p = p.pow(-exp);
            temp = temp.mul(p);
        }

        // Fix roundoff errors. (eg. floor(log10(1e9))=floor(8.9999~)=8)
        if (temp.ge(base)) {
            exp++;
            temp = temp.div(base);
        } 
        else if (temp.lt(1.0)) {
            exp--;
            temp = temp.mul(base);
        }

        if (temp.ge(base) || temp.lt(1.0)) {
            throw new RuntimeException("Can't compute exponent.");
        }

        // Handle one digit more. Used afterwards for rounding.
        int numDigits = precision + 1;

        // Extract the digits.
        for (int i = 0; i < numDigits; i++) {
            int val = (int) temp.hi;
            temp = temp.sub(val);
            temp = temp.mul(base);
            s[i] = (char) val;
        }

        if (s[0] <= 0) {
            throw new RuntimeException("Negative leading digit.");
        }

        // Fix negative digits due to roundoff error in exponent.
        for (int i = numDigits - 1; i > 0; i--) {
            if (s[i] >= 32768) {
                s[i - 1]--;
                s[i] += base;
            }
        }

        // Round, handle carry.
        if (s[precision] >= halfBase) {
            s[precision - 1]++;
            int i = precision - 1;
            while (i > 0 && s[i] >= base) {
                s[i] -= base;
                s[--i]++;
            }
        }
        s[precision] = 0;

        // If first digit became too high, shift right and
        // replace first digit by "10".
        if (s[0] >= base) {
            exp++;
            for (int i = precision; i >= 2; i--) {
                s[i] = s[i-1];
            }
            s[1] = (char)(s[0] - base);
            s[0] = 1;
        }

        // If first digit became zero, shift left.
        if (s[0] == 0) {
            exp--;
            for (int i = 0; i < precision; i++) {
                s[i] = s[i+1];
            }
        }

        // Convert to ASCII.
        for (int i = 0; i < precision; i++) {
            s[i] = DoubleDouble.BASE_36_TABLE[s[i]];
        }

        return exp;
    }

    

    

    /*------------------------------------------------------------
     *   Input
     *------------------------------------------------------------
     */
    
    /**
     * Converts a string representation of a real number into a DoubleDouble
     * value. The format accepted is similar to the standard Java real number
     * syntax. It is defined by the following regular expression:
     * <pre>
     * [<tt>+</tt>|<tt>-</tt>] {<i>digit</i>} [ <tt>.</tt> {<i>digit</i>} ] [ ( <tt>e</tt> | <tt>E</tt> ) [<tt>+</tt>|<tt>-</tt>] {<i>digit</i>}+
     * <pre>
     *
     * @param str the string to parse
     * @return the value of the parsed number
     * @throws NumberFormatException if <tt>str</tt> is not a valid
     * representation of a number
     */
    private static DoubleDouble parse(String str) throws NumberFormatException {
        int i = 0;
        int strlen = str.length();

        // skip leading whitespace
        while (Character.isWhitespace(str.charAt(i))) {
            i++;
        }
        if (str.substring(i).equals("NaN")) {
            return NaN;
        }

        // check for sign
        boolean isNegative = false;
        if (i < strlen) {
            char signCh = str.charAt(i);
            if (signCh == '-' || signCh == '+') {
                i++;
                if (signCh == '-') {
                    isNegative = true;
                }
            }
        }

        // scan all digits and accumulate into an integral value
        // Keep track of the location of the decimal point (if any)
        // to allow scaling later.
        DoubleDouble val = new DoubleDouble();

        int numDigits = 0;
        int numBeforeDec = -1;
        int exp = 0;
        while (true) {
            if (i >= strlen) {
                break;
            }
            char ch = str.charAt(i);
            i++;
            if (Character.isDigit(ch)) {
                int d = ch - '0';
                val = val.mul(10);
                val = val.add(d);
                numDigits++;
                continue;
            }
            if (ch == '.') {
                if (numBeforeDec != -1) {
                    throw new NumberFormatException("Multiple decimal dots in number string " + str);
                }
                numBeforeDec = numDigits;
                continue;
            }
            if (ch == 'e' || ch == 'E') {
                String expStr = str.substring(i);
                // this should catch any format problems with the exponent
                try {
                    exp = Integer.parseInt(expStr);
                } catch (NumberFormatException ex) {
                    throw new NumberFormatException("Invalid exponent " + expStr + " in string " + str);
                }
                break;
            }
            throw new NumberFormatException("Unexpected character '" + ch
                    + "' at position " + i
                    + " in string " + str);
        }
        DoubleDouble val2 = val;
        
        // If there was no decimal dot, then set the number of digits before
        // the decimal dot to the number of digits.
        if (numBeforeDec == -1) {
            numBeforeDec = numDigits;
        }

        // scale the number correctly
        int numDecPlaces = numDigits - numBeforeDec - exp;
        if (numDecPlaces == 0) {
            val2 = val;
        } 
        else if (numDecPlaces > 0) {
            DoubleDouble scale = TEN.pow(numDecPlaces);
            val2 = val.div(scale);
        } 
        else if (numDecPlaces < 0) {
            DoubleDouble scale = TEN.pow(-numDecPlaces);
            val2 = val.mul(scale);
        }
        
        // apply leading sign, if any
        if (isNegative) {
            return val2.neg();
        }
        return val2;
    }
}
