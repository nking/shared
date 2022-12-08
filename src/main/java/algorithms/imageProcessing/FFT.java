package algorithms.imageProcessing;

import algorithms.misc.Complex;
import algorithms.misc.MiscMath0;

/**
 * 
 * adapted from  Cormen, Leiserson, Rivest, and Stein pseudocode
 * 
 * first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.
   * 
 * @author nichole
 */
public class FFT {
    
    private boolean performNormalization = true;
    
    /**
     *
     */
    public void setToNotNormalize() {
        performNormalization = false;
    }
    
    /**
     * compute the FFT of x[], assuming its length is a power of 2.
     * runtime complexity is O(N*log(N)).
     * (adapted from  Cormen, Leiserson, Rivest, and Stein pseudocode)
     * ability added
     @param x
     @return 
     */ 
    public Complex[] fft(Complex[] x) {
        
        return fft(x, true);
    }

    /**
     * FFT from  Cormen, Leiserson, Rivest, and Stein pseudocode for iterative FFT w/ an inverse
     * ability added.  Note that the length of x must be a power of 2.
     * runtime complexity is O(N*log(N)).
     * 
     @param x
     @param forward run the transform in forward if true, else perform inverse 
     * transform.
     @return 
     */
    protected Complex[] fft(Complex[] x, boolean forward) {
       
        if (x == null || x.length == 0) {
            throw new IllegalArgumentException("xReal cannot be null or empty");
        }
     
        int n = x.length;

        if (n == 1) {
            return x;
        }
        
        if (!MiscMath0.isAPowerOf2(n)) {
            throw new IllegalArgumentException("x's length has to be a power of 2");
        }
        
        Complex[] a = bitReverseCopy(x);
        
        //TODO: could improve the speed at expense of space by caching
        // norm, end, m, eCoeff, and wn
        
        double norm = 1./Math.sqrt(n);
        
        int end = (int)(Math.log(n)/Math.log(2));
        Complex t;
        Complex u;
        Complex w;
        int m;
        double eCoeff;
        Complex wn;
        int k;
        int j;
        for (int s = 1; s <= end; s++) {
            
            m = 1 << s;
            
            eCoeff = 2. * Math.PI/(double)m;
            
            wn = forward ?
                new Complex(Math.cos(eCoeff), -Math.sin(eCoeff)) :
                new Complex(Math.cos(eCoeff), Math.sin(eCoeff));
            
            for (k = 0; k < n; k+=m) {
                
                w = new Complex(1, 0);
                
                for (j = 0; j < (m/2); j++) {
                    
                    t = w.times(a[k + j + (m/2)]);
                    
                    u = a[k + j];
                    a[k + j + (m/2)] = u.minus(t);
                    a[k + j] = u.plus(t);
                    
                    w = w.times(wn);
                }
            }
        }
        
        if (performNormalization) {
            for (int i = 0; i < a.length; i++) {
                a[i] = a[i].times(norm);
            }
        }
        
        /*
        bit-reverse-copy(a,A)
        for s=1 to lg n {
            m = 2^s
            wm = exp^(i*2*PI/m)
            for k=0 to n-1 by m {
                w=1
                for j=0 to ((m/2)-1) {
                    do t=w*A[k + j + (m/2)]
                    u = A[k + j]
                    A[k + j] = u + t
                    A[k + j + (m/2)] = u - t
                    w = w * wm
                }
            }
        }
        */
        
        return a;
    }
    
    /**
     * compute the inverse FFT of x[], assuming its length is a power of 2
     * runtime complexity is O(N).
     * (adapted from  Cormen, Leiserson, Rivest, and Stein pseudocode for forward transform)
     @param x
     @return 
     */
    public Complex[] ifft(Complex[] x) {
        
       return fft(x, false);
        
    }

    /**
     *
     @param x
     @return
     */
    protected double[] bitReverseCopy(double[] x) {
        
        int n = x.length;
        
        int nBits = MiscMath0.numberOfBits(n - 1);
                        
        double[] r = new double[n];
        int idx;
        for (int k = 0; k < n; k++) {
            
            idx = MiscMath0.bitReverse(k, nBits);
            
            r[idx] = x[k];
        }
        
        return r;
    }
    
    /**
     *
     @param x
     @return
     */
    protected double[] bitReverseCopy(int[] x) {
        
        int n = x.length;
        
        int nBits = MiscMath0.numberOfBits(n - 1);
                        
        double[] r = new double[n];
        int idx;
        for (int k = 0; k < n; k++) {
            
            idx = MiscMath0.bitReverse(k, nBits);
            
            r[idx] = x[k];
        }
        
        return r;
    }
    
    /**
     *
     @param x
     @return
     */
    protected Complex[] bitReverseCopy(Complex[] x) {
        
        int n = x.length;
        
        int nBits = MiscMath0.numberOfBits(n - 1);
                        
        Complex[] r = new Complex[n];
        int idx;
        for (int k = 0; k < n; k++) {
            
            idx = MiscMath0.bitReverse(k, nBits);
            
            r[idx] = x[k].copy();
        }
        
        return r;
    }
    
    /**
     * perform FFT on x
     * runtime complexity is O(N*log(N)).
     * (adapted from  Cormen, Leiserson, Rivest, and Stein pseudocode for forward transform).
     * Note that the result cannot be inverted because only the real portion is 
     * returned and inverse FFT needs the real and complex components.
     @param x
     @return 
     */
    public double[] fft(int[] x) {

        if (x == null || x.length == 0) {
            throw new IllegalArgumentException("xReal cannot be null or empty");
        }
     
        int n = x.length;

        if (n == 1) {
            return new double[]{x[0]};
        }
        
        if (!MiscMath0.isAPowerOf2(n)) {
            throw new IllegalArgumentException("xReal's length has to be a power of 2");
        }
        
        double[] a = bitReverseCopy(x);
        
        return _fft_body(a);
    }

    /**
     * perform FFT on x
     * runtime complexity is O(N).
     * (adapted from  Cormen, Leiserson, Rivest, and Stein pseudocode for forward transform).
     * Note that the result cannot be inverted because only the real portion is
     * returned and inverse FFT needs the real and complex components.
     @param x
     @return 
     */
    public double[] fft(double[] x) {

        if (x == null || x.length == 0) {
            throw new IllegalArgumentException("xReal cannot be null or empty");
        }

        int n = x.length;

        if (n == 1) {
            return new double[]{x[0]};
        }

        if (!MiscMath0.isAPowerOf2(n)) {
            throw new IllegalArgumentException("xReal's length has to be a power of 2");
        }

        double[] a = bitReverseCopy(x);

        return _fft_body(a);
    }

    /**
     * the main body of FFT
     @param bitReversedX
     @return 
     */
    double[] _fft_body(double[] bitReversedX) {

        if (bitReversedX == null || bitReversedX.length == 0) {
            throw new IllegalArgumentException("xReal cannot be null or empty");
        }

        int n = bitReversedX.length;

        if (!MiscMath0.isAPowerOf2(n)) {
            throw new IllegalArgumentException("xReal's length has to be a power of 2");
        }

        double[] a = bitReversedX;

        double wReal;
        double wImag;
        double tReal;
        double tImag;
        double tAbs;
        double u;
        double eCoeff;
        double wnReal;
        double wnImag;
        double norm = 1./Math.sqrt(n);
        int end = (int)(Math.log(n)/Math.log(2));
        int m;

        for (int s = 1; s <= end; s++) {

            m = 1 << s;
            eCoeff = 2. * Math.PI/(double)m;
            wnReal = Math.cos(eCoeff);
            wnImag = Math.sin(eCoeff);

            for (int k = 0; k < n; k+=m) {
                wReal = 1;
                wImag = 0;
                for (int j = 0; j < (m/2); j++) {
                    //complex multiplication:
                    tReal = wReal * a[k + j + (m/2)];
                    tImag = wImag * a[k + j + (m/2)];
                    tAbs = Math.hypot(tReal, tImag);;

                    u = a[k + j];
                    a[k + j] = (u + tAbs);
                    a[k + j + (m/2)] = (u - tAbs);

                    //complex multiplication:
                    wReal = wReal * wnReal - (wImag * wnImag);
                    wImag = wReal * wnImag + (wImag * wnReal);
                }
            }
        }

        if (performNormalization) {
            for (int i = 0; i < a.length; i++) {
                a[i] *= norm;
            }
        }

        /*
        bit-reverse-copy(a,A)
        for s=1 to lg n {
            m = 2^s
            wm = exp^(i*2*PI/m)
            for k=0 to n-1 by m {
                w=1
                for j=0 to ((m/2)-1) {
                    do t=w*A[k + j + (m/2)]
                    u = A[k + j]
                    A[k + j] = u + t
                    A[k + j + (m/2)] = u - t
                    w = w * wm
                }
            }
        }
        */

        return a;
    }
}
