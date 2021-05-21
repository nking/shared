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



public final class DoubleComplex implements Serializable {
    
    public static final long serialVersionUID = Hash64.hash("DoubleComplex_v1.0");
    
    public static final DoubleComplex ZERO = new DoubleComplex(DoubleDouble.ZERO, DoubleDouble.ZERO);
    public static final DoubleComplex ONE = new DoubleComplex(DoubleDouble.ONE, DoubleDouble.ZERO);
    public static final DoubleComplex I = new DoubleComplex(DoubleDouble.ZERO, DoubleDouble.ONE);
    
    public final DoubleDouble re;
    public final DoubleDouble im;
    
    
    
    public DoubleComplex(DoubleDouble a) {
        re = a;
        im = DoubleDouble.ZERO;
    }
    
    
    public DoubleComplex(DoubleDouble a, DoubleDouble b) {
        re = a;
        im = b;
    }
    
    
    
    public DoubleComplex(Complex a) {
        re = new DoubleDouble(a.re);
        im = new DoubleDouble(a.im);
    }
    
    
    
    
    public DoubleComplex(int a) {
        re = new DoubleDouble(a);
        im = DoubleDouble.ZERO;
    }
    
    
    
    
    public DoubleComplex(int a, int b) {
        re = new DoubleDouble(a);
        im = new DoubleDouble(b);
    }
    
    
    
    
    public DoubleComplex(long a) {
        re = new DoubleDouble(a);
        im = DoubleDouble.ZERO;
    }
    
    
    
    
    public DoubleComplex(long a, long b) {
        re = new DoubleDouble(a);
        im = new DoubleDouble(b);
    }
    
    
    
    
    public DoubleComplex(double a) {
        re = new DoubleDouble(a);
        im = DoubleDouble.ZERO;
    }
    
    
    
    
    public DoubleComplex(double a, double b) {
        re = new DoubleDouble(a);
        im = new DoubleDouble(b);
    }
    
    
    
    public DoubleComplex(String a) {
        re = new DoubleDouble(a);
        im = DoubleDouble.ZERO;
    }
    
    
    
    public DoubleComplex(String a, String b) {
        re = new DoubleDouble(a);
        im = new DoubleDouble(b);
    }
    
    
    
    
    /**
     * Returns a Complex[] array with approximate values of the supplied
     * DoubleComplex array.
     *
     * @param arr The DoubleComplex[] array to be converted to Complex[].
     * @return A Complex[] array with approximate values of the input array.
     */
    public static Complex[] toComplex(DoubleComplex[] arr) {
        Complex[] a = new Complex[arr.length];
        for (int i=0; i<arr.length; i++) {
            a[i] = new Complex(arr[i].re.hi + arr[i].re.lo, arr[i].im.hi + arr[i].im.lo);
        }
        return a;
    }
    
    
    
    /**
     * Returns a DoubleDouble[] array with all real values of the supplied
     * DoubleComplex array.
     *
     * @param arr The DoubleComplex[] array to be converted to real values.
     * @return A DoubleDouble[] array with real values from the input array.
     */
    public static DoubleDouble[] toReal(DoubleComplex[] arr) {
        DoubleDouble[] arr_re = new DoubleDouble[arr.length];
        for (int i=0; i<arr.length; i++) {
            arr_re[i] = arr[i].re;
        }
        return arr_re;
    }
    
    
    
    /**
     * Returns a DoubleDouble[] array with all imaginary values of the supplied
     * DoubleComplex array.
     *
     * @param arr The DoubleComplex[] array to be converted to imaginary values.
     * @return A DoubleDouble[] array with imaginary values from the input array.
     */
    public static DoubleDouble[] toImag(DoubleComplex[] arr) {
        DoubleDouble[] arr_im = new DoubleDouble[arr.length];
        for (int i=0; i<arr.length; i++) {
            arr_im[i] = arr[i].im;
        }
        return arr_im;
    }
    
    
    
    public Complex complexValue() {
        return new Complex(re.hi+re.lo, im.hi+im.lo);
    }



    public boolean isZero() {
        return re.isZero() && im.isZero();
    }



    public boolean isOne() {
        return re.isOne() && im.isZero();
    }



    public boolean isNear(DoubleComplex b, double eps) {
        double a = abs1().doubleValue();
        double diff = sub(b).abs1().doubleValue();
        return diff < a*eps;
    }



    public boolean isNear(DoubleDouble b, double eps) {
        double a = abs1().doubleValue();
        double diff = sub(b).abs1().doubleValue();
        return diff < a*eps;
    }



    public DoubleComplex add(DoubleComplex a) {
        return new DoubleComplex(re.add(a.re), im.add(a.im));
    }

    
    
    public DoubleComplex add(DoubleDouble a) {
        return new DoubleComplex(re.add(a), im);
    }

    
    
    public DoubleComplex add1() {
        return new DoubleComplex(re.add1(), im);
    }

    
    
    public DoubleComplex add(DoubleDouble a, DoubleDouble b) {
        return new DoubleComplex(a == null ? re : re.add(a), 
                                 b == null ? im : im.add(b));
    }



    public DoubleComplex sub(DoubleComplex a) {
        return new DoubleComplex(re.sub(a.re), im.sub(a.im));
    }

    
    
    public DoubleComplex sub(DoubleDouble a) {
        return new DoubleComplex(re.sub(a), im);
    }

    
    
    public DoubleComplex sub1() {
        return new DoubleComplex(re.sub1(), im);
    }

    
    
    public DoubleComplex sub(DoubleDouble a, DoubleDouble b) {
        return new DoubleComplex(a == null ? re : re.sub(a), 
                                 b == null ? im : im.sub(b));
    }
    
    
    
    public DoubleComplex mul(DoubleComplex b) {
        DoubleDouble r = re.mul(b.re).sub(im.mul(b.im));
        DoubleDouble i = re.mul(b.im).add(im.mul(b.re));
        return new DoubleComplex(r, i);
    }
    
    
    
    public DoubleComplex mul(DoubleDouble b) {
        DoubleDouble r = re.mul(b);
        DoubleDouble i = im.mul(b);
        return new DoubleComplex(r, i);
    }
    
    
    
    public DoubleComplex mul(double b) {
        DoubleDouble r = re.mul(b);
        DoubleDouble i = im.mul(b);
        return new DoubleComplex(r, i);
    }
    
    
    
    public DoubleComplex mul(DoubleDouble bre, DoubleDouble bim) {
        if (bre != null) {
            if (bim != null) {
                DoubleDouble r = re.mul(bre).sub(im.mul(bim));
                DoubleDouble i = re.mul(bim).add(im.mul(bre));
                return new DoubleComplex(r, i);
            }
            else {
                DoubleDouble r = re.mul(bre);
                DoubleDouble i = im.mul(bre);
                return new DoubleComplex(r, i);
            }
        }
        else {
            if (bim != null) {
                DoubleDouble r = im.mul(bim).neg();
                DoubleDouble i = re.mul(bim);
                return new DoubleComplex(r, i);
            }
            else {
                return DoubleComplex.ZERO;
            }
        }
    }
    
    
    
    
    
    
    public DoubleComplex div(DoubleComplex b) {
        DoubleDouble abre = b.re.abs();
        DoubleDouble abim = b.im.abs();
        if (abre.compareTo(abim) > 0) {
            DoubleDouble r = b.im.div(b.re);
            DoubleDouble den = b.re.add(r.mul(b.im));
            DoubleDouble cre = re.add(r.mul(im)).div(den);
            DoubleDouble cim = im.sub(r.mul(re)).div(den);
            return new DoubleComplex(cre, cim);
        }
        else {
            DoubleDouble r = b.re.div(b.im);
            DoubleDouble den = b.im.add(r.mul(b.re));
            DoubleDouble cre = re.mul(r).add(im).div(den);
            DoubleDouble cim = im.mul(r).sub(re).div(den);
            return new DoubleComplex(cre, cim);
        }
    }
    
    
    
    public DoubleComplex div(DoubleDouble b) {
        return new DoubleComplex(re.div(b), im.div(b));
    }
    
    
    
    public DoubleComplex div(DoubleDouble bre, DoubleDouble bim) {
        if (bre != null) {
            if (bim != null) {
                DoubleDouble abre = bre.abs();
                DoubleDouble abim = bim.abs();
                if (abre.compareTo(abim) > 0) {
                    DoubleDouble r = bim.div(bre);
                    DoubleDouble den = bre.add(r.mul(bim));
                    DoubleDouble cre = re.add(r.mul(im)).div(den);
                    DoubleDouble cim = im.sub(r.mul(re)).div(den);
                    return new DoubleComplex(cre, cim);
                }
                else {
                    DoubleDouble r = bre.div(bim);
                    DoubleDouble den = bim.add(r.mul(bre));
                    DoubleDouble cre = re.mul(r).add(im).div(den);
                    DoubleDouble cim = im.mul(r).sub(re).div(den);
                    return new DoubleComplex(cre, cim);
                }
            }
            else {
                return new DoubleComplex(re.div(bre), im.div(bre));
            }
        }
        else {
            if (bim != null) {
                return new DoubleComplex(im.div(bim), re.div(bim).neg());
            }
            else {
                // Throws exception!
                return div(new DoubleDouble(0));
            }
        }
    }
    
    
    
    
    
    
    public DoubleComplex recip() {
        DoubleDouble are = re.abs();
        DoubleDouble aim = im.abs();
        if (are.compareTo(aim) > 0) {
            DoubleDouble r = im.div(re);
            DoubleDouble den = re.add(r.mul(im));
            
            DoubleDouble cre = den.recip();
            DoubleDouble cim = r.div(den).neg();
            return new DoubleComplex(cre, cim);
        }
        else {
            DoubleDouble r = re.div(im);
            DoubleDouble den = im.add(r.mul(re));
            DoubleDouble cre = r.div(den);
            DoubleDouble cim = den.recip().neg();
            return new DoubleComplex(cre, cim);
        }
    }
    


    public DoubleComplex neg() {
        return new DoubleComplex(re.neg(), im.neg());
    }
    


    public DoubleComplex conj() {
        return new DoubleComplex(re, im.neg());
    }
    
    
    
    
    public DoubleDouble abs() {
        if (re.isZero()) {
            return im.abs();
        }
        if (im.isZero()) {
            return re.abs();
        }
        
        DoubleDouble x = re.abs();
        DoubleDouble y = im.abs();
        
        if (x.compareTo(y) > 0) {
            DoubleDouble temp = y.div(x);
            return temp.sqr().add1().sqrt().mul(x);
        }
        else {
            DoubleDouble temp = x.div(y);
            return temp.sqr().add1().sqrt().mul(y);
        }
    }
    
    
    
    
    public DoubleDouble abs1() {
        DoubleDouble x = re.abs();
        DoubleDouble y = im.abs();
        return x.add(y);
    }
    
    
    
    
    public DoubleDouble absinf() {
        DoubleDouble x = re.abs();
        DoubleDouble y = im.abs();
        return x.compareTo(y) > 0 ? x : y;
    }
    
    
    
    
    public DoubleDouble sqrabs() {
        return re.sqr().add(im.sqr());
    }
    
    
    
    public DoubleDouble real() {
        return re;
    }
    
    
    
    public DoubleDouble imag() {
        return im;
    }

    
    
    
    public Complex ComplexValue() {
        return new Complex(re.doubleValue(), im.doubleValue());
    }


    
    
    @Override
    public String toString() {
        return toString(10, 0, false);
    }



    
    public String toString(int radix, int maxPrecision, boolean sci) {
        if (im.isZero()) {
            return "" + re.toString(radix, maxPrecision, sci);
        }
        if (re.isZero()) {
            if (im.isNegative()) {
                return "-i*" + im.neg().toString(radix, maxPrecision, sci);
            }
            return "i*" + im.toString(radix, maxPrecision, sci);
        }
        
        DoubleDouble r = re;
        if (r.isNegative()) {
            r = r.neg();
        }
        
        DoubleDouble i = im;
        boolean imIsNeg = false;
        if (i.isNegative()) {
            i = i.neg();
            imIsNeg = true;
        }
        
        double threshold = (double)(1l << 52);
        threshold *= threshold;
        if (maxPrecision > 0) {
            double precthresh = Math.pow(radix, maxPrecision);
            if (precthresh < threshold) {
                threshold = precthresh;
            }
        }
        
        
        DoubleDouble r_i = r.div(i);
        if (r_i.gt(threshold)) {
            return "" + re.toString(radix, maxPrecision, sci);
        }
        if (r_i.lt(1.0/threshold)) {
            return (imIsNeg ? "i*" : "-i*") + i.toString(radix, maxPrecision, sci);
        }
        
        return "" + re.toString(radix, maxPrecision, sci) + (imIsNeg ? " + i*" : " - i*") + i.toString(radix, maxPrecision, sci);
    }
    
    
/*

FP arg(const cplx &a)
{
  return atan2(a.im, a.re);
}



cplx euler(FP phi)
{
  FP re, im;
  sincos(phi, re, im);
  return cplx(re, im);
}

*/
    
}
