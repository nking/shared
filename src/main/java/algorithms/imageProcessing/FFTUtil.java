package algorithms.imageProcessing;

import algorithms.misc.Complex;
import java.util.Arrays;
import thirdparty.ca.uol.aig.fftpack.Complex1D;
import thirdparty.ca.uol.aig.fftpack.ComplexDoubleFFT;

/**
 *
 * first implemented in project
     http://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)
   then moved here to share with other projects.

   NOTE: consider adding Sparse FFT:
     http://groups.csail.mit.edu/netmit/sFFT/

 * @author nichole
 */
public class FFTUtil {
   
    /**
     *
     @param input
     @param forward
     @return
     */
    public Complex[][] create2DFFT(double[][] input, boolean forward) {

        // performs normalization by default
        return create2DFFT(input, true, forward);
    }

    /**
     * perform fft on input.
     @param input
     @param doNormalize
     @param forward
     @return
     */
    public Complex[][] create2DFFT(final double[][] input, boolean doNormalize,
        boolean forward) {

        Complex[][] input2 = new Complex[input.length][];
        for (int i = 0; i < input.length; ++i) {
            input2[i] = new Complex[input[0].length];
            for (int j = 0; j < input[0].length; ++j) {
                input2[i][j] = new Complex(input[i][j], 0);
            }
        }

        return create2DFFT(input2, doNormalize, forward);
    }

    /**
     *
     @param input
     @param forward
     @return
     */
    public Complex[][] create2DFFT(Complex[][] input, boolean forward) {

        // performs normalization by default
        return create2DFFT(input, true, forward);
    }

    /**
     * runtime complexity: is O(N*lg_2(N)) for N not power of 2,
     * else is
     *
     * perform fft on input.
     @param input
     @param doNormalize
     @param forward
     @return
     */
    public Complex[][] create2DFFT(final Complex[][] input, boolean doNormalize,
        boolean forward) {

        final int n0 = input.length;
        final int n1 = input[0].length;

        int nn0 = 1 << (int)(Math.ceil(Math.log(n0)/Math.log(2)));
        int nn1 = 1 << (int)(Math.ceil(Math.log(n1)/Math.log(2)));

        if (nn0 > n0 || nn1 > n1) {
            Complex1D[] input2 = copyToComplex1D(input);
            Complex1D[] output = create2DFFT2(input2, doNormalize, forward);
            Complex[][] output2 = copyToComplex(output);
            return output2;
        }

        Complex[][] output = copy(input);

        // padding is at front of cols and rows

        FFT fft = new FFT();
        if (!doNormalize) {
            fft.setToNotNormalize();
        }

        // ----- perform FFT by dimension 0 -----
        for (int i0 = 0; i0 < nn0; i0++) {
            if (forward) {
                output[i0] = fft.fft(output[i0]);
            } else {
                output[i0] = fft.ifft(output[i0]);
            }
        }

        // re-use array for the FFT by dimension 1
        Complex[] tmp = new Complex[nn0];

        /*
        nn0
         |
        \|/
        [0]  ..........nn1-1
        [1]  ..........nn1-1
        */

        // ----- perform the FFT on dimension 1 ------
        for (int i1 = 0; i1 < nn1; ++i1) {

            // store each column in tmp array and perform fft on it then
            // recopy values back into columns
            for (int i0 = 0; i0 < nn0; ++i0) {
                tmp[i0] = output[i0][i1];
            }

            if (forward) {
                tmp = fft.fft(tmp);
            } else {
                tmp = fft.ifft(tmp);
            }

            for (int i0 = 0; i0 < nn0; ++i0) {
                output[i0][i1] = tmp[i0];
            }
        }

        return output;
    }
    
    /**
     * runtime complexity: is O(N*lg_2(N)) for N not power of 2,
     * else is
     *
     * perform fft on input.
     @param input
     @param forward
     @return
     */
    public Complex[] create1DFFT(final Complex[] input, boolean forward) {

        Complex1D input2 = copyToComplex1D(input);
        Complex1D output = create1DFFT2(input2, forward);
        Complex[] output2 = copyToComplex(output);
        
        return output2;
    }

    /**
     * runtime complexity: is O(N*lg_2(N)) for N not power of 2,
     * else is
     *
     * perform fft on input.
     @param input
     @param forward
     @return
     */
    public Complex[] create1DFFTNormalized(final Complex[] input, boolean forward) {

        Complex1D input2 = copyToComplex1D(input);
        Complex1D output = create1DFFT2(input2, forward);
        Complex[] output2 = copyToComplex(output);
        normalize(output2);
        return output2;
    }

    private void normalize(Complex[] a) {
        double norm = 1./Math.sqrt(a.length);
        for (int i = 0; i < a.length; ++i) {
            a[i] = a[i].times(norm);
        }
    }

    /**
     * runtime complexity: is O(N*lg_2(N)) for N not power of 2,
     * else is
     *
     * perform fft on input.
     @param input
     @param forward
     @return
     */
    public Complex[] create1DFFT(final double[] input, boolean forward) {

        Complex1D input2 = copyToComplex1D(input);
        Complex1D output = create1DFFT2(input2, forward);
        Complex[] output2 = copyToComplex(output);
        
        return output2;
    }

    /**
     * runtime complexity: is O(N*lg_2(N))
     *
     * perform fft on input.
     @param input
     @param forward
     @return
     */
    public Complex[] create1DFFTNormalized(final double[] input, boolean forward) {

        Complex1D input2 = copyToComplex1D(input);
        Complex1D output = create1DFFT2(input2, forward);
        Complex[] output2 = copyToComplex(output);
        normalize(output2);
        return output2;
    }

    /**
     * perform a 2-dimension FFT using the JFFTPack library.
     *
     @param input double array of complex data in format double[nRows][2*nColumns]
     * where the column elements are alternately the complex real number and the
     * complex imaginary number.
     @param forward
     @return two dimensional complex array of size Complex[nRows][input.nCols/2)
     */
    public Complex1D[] create2DFFT2(Complex1D[] input, boolean forward) {

        // perform normalization by default
        return create2DFFT2(input, true, forward);
    }

    /**
     * perform a 2-dimension FFT using the JFFTPack library.
     *
     * runtime complexity: is O(N*lg_2(N)) for N not power of 2.
     *
     @param input double array of complex data in format double[nRows][2*nColumns]
     * where the column elements are alternately the complex real number and the
     * complex imaginary number.
     @param performNormalization
     @param forward
     @return two dimensional complex array of size Complex[nRows][input.nCols/2)
     */
    public Complex1D[] create2DFFT2(Complex1D[] input, boolean performNormalization,
        boolean forward) {

        final int n0 = input.length;
        final int n1 = input[0].x.length;

        Complex1D[] output = Arrays.copyOf(input, input.length);

        ComplexDoubleFFT fft1 = new ComplexDoubleFFT(n1);

        final double norm1 = 1./Math.sqrt(n1);

        // ----- perform FFT by dimension 0 -----
        for (int i0 = 0; i0 < n0; i0++) {

            if (forward) {
                fft1.ft(output[i0]);
            } else {
                fft1.bt(output[i0]);
            }

            // normalize the data
            if (performNormalization) {
                Complex1D a = output[i0];
                for (int idx = 0; idx < a.x.length; ++idx) {
                    a.x[idx] *= norm1;
                    a.y[idx] *= norm1;
                }
            }
        }

        // re-use array for the FFT by dimension 1 (across rows)
        Complex1D tmp = new Complex1D();
        tmp.x = new double[n0];
        tmp.y = new double[n0];

        ComplexDoubleFFT fft0 = new ComplexDoubleFFT(n0);

        final double norm0 = performNormalization ? (1./Math.sqrt(n0)) : 1.;

        // ----- perform the FFT on dimension 1 ------
        for (int i1 = 0; i1 < n1; ++i1) {

            // store each column in tmp array and perform fft on it then
            // recopy values back into columns
            for (int i0 = 0; i0 < n0; ++i0) {
                tmp.x[i0] = output[i0].x[i1];
                tmp.y[i0] = output[i0].y[i1];
            }

            if (forward) {
                fft0.ft(tmp);
            } else {
                fft0.bt(tmp);
            }

            for (int i0 = 0; i0 < n0; ++i0) {
                output[i0].x[i1] = tmp.x[i0] * norm0;
                output[i0].y[i1] = tmp.y[i0] * norm0;
            }
        }

        return output;
    }
    
    /**
     * perform a 1-dimension FFT using the JFFTPack library.
     *
     * runtime complexity: is O(N*lg_2(N)) for N not power of 2.
     *
     @param input double array of complex data in format double[nRows][2*nColumns]
     * where the column elements are alternately the complex real number and the
     * complex imaginary number.
     @param forward
     @return two dimensional complex array of size Complex[nRows][input.nCols/2)
     */
    public Complex1D create1DFFT2(Complex1D input, boolean forward) {

        final int n0 = input.x.length;

        Complex1D output = new Complex1D();
        output.x = Arrays.copyOf(input.x, input.x.length);
        output.y = Arrays.copyOf(input.y, input.y.length);

        ComplexDoubleFFT fft1 = new ComplexDoubleFFT(n0);

        if (forward) {
            fft1.ft(output);
        } else {
            fft1.bt(output);
        }

        return output;
    }

    /**
     *
     @param input
     @return
     */
    public Complex1D[] copyToComplex1D(Complex[][] input) {

        int n0 = input.length;
        int n1 = input[0].length;

        Complex1D[] output = new Complex1D[n0];
        for (int i = 0; i < n0; ++i) {
            output[i] = new Complex1D();
            output[i].x = new double[n1];
            output[i].y = new double[n1];
            for (int j = 0; j < n1; ++j) {
                output[i].x[j] = input[i][j].re();
                output[i].y[j] = input[i][j].im();
            }
        }

        return output;
    }
    
    /**
     *
     @param input
     @return
     */
    public Complex1D copyToComplex1D(Complex[] input) {

        int n0 = input.length;

        Complex1D output = new Complex1D();
        output.x = new double[n0];
        output.y = new double[n0];
        for (int i = 0; i < n0; ++i) {
            output.x[i] = input[i].re();
            output.y[i] = input[i].im();
        }

        return output;
    }
    
    /**
     *
     @param input
     @return
     */
    public Complex1D copyToComplex1D(double[] input) {

        int n0 = input.length;

        Complex1D output = new Complex1D();
        output.x = Arrays.copyOf(input, input.length);
        output.y = new double[n0];

        return output;
    }
    
    /**
     *
     @param input
     @return
     */
    public Complex[] copyToComplex(Complex1D input) {

        int n0 = input.x.length;

        Complex[] output = new Complex[n0];
        for (int i = 0; i < n0; ++i) {
            output[i] = new Complex(input.x[i], input.y[i]);
        }

        return output;
    }
    
    /**
     *
     @param input
     @return
     */
    public Complex[][] copyToComplex(Complex1D[] input) {

        int n0 = input.length;
        int n1 = input[0].x.length;

        Complex[][] output = new Complex[n0][];
        for (int i = 0; i < n0; ++i) {
            output[i] = new Complex[n1];
            for (int j = 0; j < n1; ++j) {
                output[i][j] = new Complex(input[i].x[j], input[i].y[j]);
            }
        }

        return output;
    }

    /**
     *
     @param input
     @return
     */
    public Complex[][] copy(Complex[][] input) {

        int n0 = input.length;

        Complex[][] output = new Complex[n0][];
        for (int i = 0; i < n0; ++i) {
            output[i] = Arrays.copyOf(input[i], input[i].length);
        }

        return output;
    }

    /**
     *
     @param fftData
     @return
     */
    public double[] extractAbs(Complex[] fftData) {
        double[] a = new double[fftData.length];
        for (int i = 0; i < fftData.length; ++i) {
            a[i] = fftData[i].abs();
        }
        return a;
    }
}
