package thirdparty.net.oelen.polarith;

import java.io.Serializable;



// An immutable Complex class, supporting basic arithmetic
// on complex numbers. The public members re, im are final
// and hence cannot be modified.

public final class Complex implements Serializable {
    private static final long serialVersionUID = Hash64.hash("Complex_v1.0");
    
    public static final Complex ZERO = new Complex(0.0, 0.0);
    public static final Complex ONE = new Complex(1.0, 0.0);
    public static final Complex I = new Complex(0.0, 1.0);
    
    public final double re;
    public final double im;

    public Complex(double re, double im) {
        this.re = re;
        this.im = im;
    }

    public Complex(double re) {
        this.re = re;
        this.im = 0.0;
    }

    public Complex() {
        this.re = 0.0;
        this.im = 0.0;
    }



    public Complex add(double re) {
        return new Complex(this.re + re, this.im);
    }


    public Complex add(double re, double im) {
        return new Complex(this.re + re, this.im + im);
    }


    public Complex add(Complex a) {
        return new Complex(this.re + a.re, this.im + a.im);
    }



    
    public Complex add1() {
        return new Complex(this.re + 1.0, this.im);
    }




    public Complex sub(double re) {
        return new Complex(this.re - re, this.im);
    }


    public Complex sub(double re, double im) {
        return new Complex(this.re - re, this.im - im);
    }


    public Complex sub(Complex a) {
        return new Complex(this.re - a.re, this.im - a.im);
    }




    public Complex sub1() {
        return new Complex(this.re - 1.0, this.im);
    }




    public Complex mul(double re) {
        return new Complex(this.re * re, this.im * re);
    }

    public Complex mul(double re, double im) {
        return new Complex(this.re * re - this.im * im, this.im * re + this.re * im);
    }

    public Complex mul(Complex a) {
        return new Complex(this.re * a.re - this.im * a.im, this.im * a.re + this.re * a.im);
    }





    public Complex sqr() {
        return new Complex(this.re * this.re - this.im * this.im, 2 * this.im * this.re);
    }




    public Complex div(double re) {
        return new Complex(this.re/re, this.im/re);
    }

    
    
    
    private static final double THRESH_MAX = 0.5 * Math.sqrt(Double.MAX_VALUE);
    private static final double THRESH_MIN = 1.0/THRESH_MAX;
    
    public Complex div(double re, double im) {
        double are = Math.abs(re);
        double aim = Math.abs(im);

        if (are + aim > THRESH_MAX) {
            // Use a special algorithm, which assures that
            // no overflow occurs, but which is slower.
            if (are > aim) {
                double im_re = im/re;
                double rr = 1.0/(re + im_re*im);
                return new Complex((this.re + this.im*im_re)*rr, (this.im - this.re*im_re)*rr);
            }
            else {
                double re_im = re/im;
                double rr = 1.0/(re*re_im + im);
                return new Complex((this.re*re_im + this.im)*rr, (this.im*re_im - this.re)*rr);
            }
        }

        if (are + aim < THRESH_MIN) {
            // Use a special algorithm, which assures that
            // no underflow occurs, but which is slower.
            if (are < aim) {
                double im_re = im/re;
                double rr = 1.0/(re + im_re*im);
                return new Complex((this.re + this.im*im_re)*rr, (this.im - this.re*im_re)*rr);
            }
            else {
                double re_im = re/im;
                double rr = 1.0/(re*re_im + im);
                return new Complex((this.re*re_im + this.im)*rr, (this.im*re_im - this.re)*rr);
            }
        }

        double rr = 1.0/(re*re + im*im);
        return new Complex((this.re*re + this.im*im)*rr, (this.im*re - this.re*im)*rr);
    }

    
    
    
    public Complex div(Complex a) {
        double are = Math.abs(a.re);
        double aim = Math.abs(a.im);

        if (are + aim > THRESH_MAX) {
            return div(a.re, a.im);
        }

        if (are + aim < THRESH_MIN) {
            return div(a.re, a.im);
        }

        double rr = 1.0/(a.re*a.re + a.im*a.im);
        return new Complex((this.re*a.re + this.im*a.im)*rr, (this.im*a.re - this.re*a.im)*rr);
    }
    
    
    
    
    
    
    public Complex recip() {
        double are = Math.abs(re);
        double aim = Math.abs(im);
        if (are > aim) {
            double r = im / re;
            double den = re + r * im;
            
            double cre = 1.0/den;
            double cim = -r/den;
            return new Complex(cre, cim);
        }
        else {
            double r = re / im;
            double den = im + r * re;
            double cre = r / den;
            double cim = -1.0/den;
            return new Complex(cre, cim);
        }
    }
    
    



    public double abs() {
        double are = Math.abs(re);
        double aim = Math.abs(im);
        
        if (aim == 0.0) {
            return are;
        }
        if (are == 0.0) {
            return aim;
        }
        
        // We do not simply use sqrt(re*re+im*im). If we used that
        // then there would be an intermediate overflow for values
        // in the order of magnitude of 10^154. The code below does
        // a little smarter arithmetic and only overflows if the
        // number really cannot be represented anymore in a double
        // (which is in the order of magnitude of 10^308).
        
        if (are > aim) {
            double im_re = aim/are;
            return are * Math.sqrt(1.0 + im_re*im_re);
        }
        else {
            double re_im = are/aim;
            return aim * Math.sqrt(re_im*re_im + 1.0);
        }
    }



    public double abs1() {
        double are = Math.abs(re);
        double aim = Math.abs(im);
        return are + aim;
    }



    public double real() {
        return this.re;
    }



    public double imag() {
        return this.im;
    }
    
    
    
    public Complex neg() {
        return new Complex(-re, -im);
    }
    
    
    
    public Complex conj() {
        return new Complex(re, -im);
    }
    
    
    
    public boolean isZero() {
        return re==0.0 && im==0.0;
    }
    
    
    
    @Override
    public String toString() {
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
    
    
    
    public static double[] toReal(Complex[] arr) {
        double[] arr_re = new double[arr.length];
        for (int i=0; i<arr.length; i++) {
            arr_re[i] = arr[i].re;
        }
        return arr_re;
    }
    
    
    
    public static double[] toImag(Complex[] arr) {
        double[] arr_im = new double[arr.length];
        for (int i=0; i<arr.length; i++) {
            arr_im[i] = arr[i].im;
        }
        return arr_im;
    }
}
