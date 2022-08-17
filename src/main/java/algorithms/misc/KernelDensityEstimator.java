package algorithms.misc;

import algorithms.imageProcessing.FFT;

import java.util.Arrays;

/**
 * TODO: consider implementing the Improved Sheather - Jones algorithm for estimating bandwidth:.
 * can see implementations
 * https://github.com/tommyod/KDEpy/blob/master/KDEpy/bw_selection.py
 * and
 * https://kdepy.readthedocs.io/en/latest/_modules/KDEpy/FFTKDE.html#FFTKDE
 *
 * TODO: consider implementing
 * "Fast & Accurate Gaussian Kernel Density Estimation", Jeffrey Heer, University of Washiington.
 * https://idl.cs.washington.edu/files/2021-FastKDE-VIS.pdf
 *
 * TODO: consider using MSER in a 2-D KDE estimator
 *
 * TODO: explore https://github.com/TasCL/cpda/ one day.
 * Probability Density Approximation.
 * The authors are Yi-Shin Lin yishin.lin@utas.edu.au, Andrew Heathcote, William Holmes
 * The repository implements some of
 * Holmes, W. (2015). A practical guide to the Probability Density Approximation (PDA) with improved implementation
 * and error characterization. Journal of Mathematical Psychology, 68-69, 13--24, doi:
 * http://dx.doi.org/10.1016/j.jmp.2015.08.006.
 */
public class KernelDensityEstimator {

    final static double BIG = 0.8 * Double.MAX_VALUE;

    /**
     * assuming a zero-centered mean gaussian kernel, return the bandwidth which minimizes the
     * mean integrated squared error (MISE).
     * <pre>
     *     Silverman 1981, "Kernel Density Estimation Using theFast Fourier Transform"
     * </pre>
     * @param standardDeviation
     * @param nSample
     * @return
     */
    public static double optimalBandwidthGaussianKernel(double standardDeviation, int nSample) {
        return 1.06 * standardDeviation * Math.pow(nSample, -1./5);
    }

    /**
     * assuming a zero-centered mean gaussian kernel, return Silverman's "rule of thumb" bandwidth.
     * This is considered an improvement over the optimal bandwidth for gaussian kernels.
     * <pre>
     *     Silverman's ‘rule of thumb’, Silverman (1986, page 48, eqn (3.31))).
     *     and
     *     https://rdrr.io/r/stats/bandwidth.html
     </pre>
     * @param standardDeviation
     * @param nSample
     * @param IQR the interquartile range which is the ordered statistic divided into 4 equal parts, each summed,
     *            then iqr = the 3rd quartile minus the first quartile = sum[3] - sum[0].
     *            see MiscMath0.calcMedianAndIQR()
     * @return
     */
    public static double ruleOfThumbBandwidthGaussianKernel(double standardDeviation, double IQR, int nSample) {
        return 0.9 * Math.min(standardDeviation, IQR/1.34) * Math.pow(nSample, -1./5);
    }

    public static class KDE {

        /**
         * FFT of the histogram of the data.
         u(s) = (1./sqrt(2*pi)) * (1/n) * sum j=1 to n of ( exp(i*s*X_j) ) where i is imaginary
         and n is the length of the data.
         Note that the data are a histogram of the original data.
         <pre>
         reference:
         Silverman, B. W. (1982). Algorithm as 176: Kernel density estimation using the fast Fourier transform. Journal of the Royal Statistical Society. Series C (Applied Statistics), 31(1), 93-99. https://dx.doi.org/10.2307/2347084.
         https://www.jstor.org/stable/2347084#metadata_info_tab_contents
         </pre>
         */
        public Complex[] u;

        /**
         * the x bins of the histogram of the data used in creating u.
         */
        public double[] hx;

        /**
         * kernel bandwidth
         */
        public double h;

        /**
         * the kde calculated from u and h.
         * kde = inverse FFT(f_n(s))
         *     = inverse FFT( exp(-(0.5) * h^2 * s^2) * u(s)  )
         */
        public double[] kde;
    }

    /**
     * <pre>
     *      reference:
     *      Silverman, B. W. (1982). Algorithm as 176: Kernel density estimation using the fast Fourier transform. Journal of the Royal Statistical Society. Series C (Applied Statistics), 31(1), 93-99. https://dx.doi.org/10.2307/2347084.
     *      https://www.jstor.org/stable/2347084#metadata_info_tab_contents
     *      </pre>
     * @param x zero-centered data
     * @param h bandwidth
     * @return
     */
    public static KDE viaFFTGaussKernel(double[] x, double h, int histNBins, double histBinWidth, double histMinBin,
                                        double histMaxBin) {

        // the histogram length (hist[1].length) is a power of 2 and the range is enlarged by at least 2*3*h
        double[][] hist = createFineHistogram(x, h, histNBins, histBinWidth, histMinBin, histMaxBin);

        //assert(assertHistRange(hist[0], MiscMath0.getMinMax(x), h) == true);

        Complex[] yHist = convertToComplex(hist[1]);

        // normalization is performed by default:
        FFT fft = new FFT();

        // u is the portion that can be re-used on subsequent iterations.  e.g. when iterating
        //   to find minimum bandwidth h.
        //double[] u = fft.fft(hist[1]);
        Complex[] u = fft.fft(yHist);

        assert(u.length == hist[0].length);

        return viaFFTGaussKernel(u, hist[0], h);
    }

    /**
     * <pre>
     *      reference:
     *      Silverman, B. W. (1982). Algorithm as 176: Kernel density estimation using the fast Fourier transform. Journal of the Royal Statistical Society. Series C (Applied Statistics), 31(1), 93-99. https://dx.doi.org/10.2307/2347084.
     *      https://www.jstor.org/stable/2347084#metadata_info_tab_contents
     *      </pre>
     * @param x data observed
     * @param h bandwidth
     * @return
     */
    public static KDE viaFFTGaussKernel(double[] x, double h) {

        //TODO: consider using PeriodicFFT from the curvature scale space project instead of FFT
        //     to avoid edge effects.

        // the histogram length (hist[1].length) is a power of 2 and the range is enlarged by at least 2*3*h
        double[][] hist = createFineHistogram(x, h);

        //assert(assertHistRange(hist[0], MiscMath0.getMinMax(x), h) == true);

        Complex[] yHist = convertToComplex(hist[1]);

        // normalization is performed by default:
        FFT fft = new FFT();

        // u is the portion that can be re-used on subsequent iterations.  e.g. when iterating
        //   to find minimum bandwidth h.
        //double[] u = fft.fft(hist[1]);
        Complex[] u = fft.fft(yHist);

        assert(u.length == hist[0].length);

        return viaFFTGaussKernel(u, hist[0], h);

        /*
        summary of the paper:
        note that the f_n(x) is highly inefficient to use directly on a grid of points.
             where f_n(x) = (1./(n*h)) * sum i=1 to n of (K((x - X_i)/h)   EQN (1)

          first, discretize the data to a fine grid
          then use FFT to convolve the data w/ the kernel to calculate .
            *Take FFT of eqn (1):  FFT(f_n(s)) = (sqrt(2*pi)) * FFT(K(h*s))*u(s)
               where u(s) is the fourier transform of the data
               u(s) = (1./sqrt(2*pi)) * (1/n) * sum j=1 to n of ( exp(i*s*X_j) ) where i is imaginary
               ** A discrete approx of u(s) is found by constructing a histogram on a grid of 2^k cells
                  and then apply the FFT to it.
                  (follow Gentelman and Sande as imple. by Monro 1976
               ** The Munroe 1966 FFT takes an array of X(M) of M=2^k real values and returns
                  their discrete fourier transform Y stored with the real parts of Y_0 to Y_{M/2} in
                  locations X(1) to X((M/2)_1) and the imaginary parts of Y_1 to Y_{(M/2)-1} stored
                  M/2 locations above their corresponding real parts.
            *Substitute the fourier transform of the Gaussian Kernel K(h*s)
               FFT(f_n(s)) = exp(-(0.5) * h^2 * s^2) * u(s)   EQN (3)
            *Then f_n(s) = inverse FFT(f_n(s))
             negative values of f_n are set to 0.
               ** note that if several different uses of this algorithm for different h are employed,
               that the discrete FFT has to be calculated only once and can be reused.
               ** The algorithm also avoids exponential underflow by setting
                  FFT(f_n(s)) equal to 0 if 0.5*h*h*s*s is > BIG which is large for the machine.
          The discrete fourier transform should be performed on an interval larger than the interval of
             interest because of wrap around edge conditions.  they recommend enlarging the interval
             by 3*h at each end.  Note that the enlargement is not needed for circular data.
         */
    }

    private static Complex[] convertToComplex(double[] a) {
        Complex[] c = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            c[i] = new Complex(a[i], 0);
        }
        return c;
    }

    /**
     * calculate a fine resolution histogram for x, for a larger data range than x's.
     * (NOTE that if the data are circular, a method can be created to calculate the histogram
     * using the same data range as x, not larger)
     * @param x
     * @return a 2-dimensional array of the histogram where hist[0] holds the centers of the histogram bins,
     * and hist[1] holds the counts within the bins.
     */
    protected static double[][] createFineHistogram(double[] x, double h) {
        // the number of bins need to be a power of 2, and larger than x.length.
        int nBins = (int)Math.pow(2, Math.ceil(Math.log(x.length * 2)/Math.log(2)));

        // unless the data are circular, the range has to be larger than the range of x in order to
        // avoid wrap around edge conditions
        double[] minMaxX = MiscMath0.getMinMax(x);
        double range = minMaxX[1] - minMaxX[0];
        double dr;
        if (0.5 * range > 3.*h) {
            dr = 0.5 * range;
        } else {
            dr = 4. * h;
        }
        double min = minMaxX[0] - dr;
        double max = minMaxX[1] + dr;

        // nBins = (maxX - minX)/binWidth
        double binWidth = nBins/(max - min);

        return createFineHistogram(x, h, nBins, binWidth, min, max);
    }

    /**
     * calculate a fine resolution histogram for x, for a larger data range than x's.
     * (NOTE that if the data are circular, a method can be created to calculate the histogram
     * using the same data range as x, not larger)
     * @param x
     * @return a 2-dimensional array of the histogram where hist[0] holds the centers of the histogram bins,
     * and hist[1] holds the counts within the bins.
     */
    protected static double[][] createFineHistogram(double[] x, double h, int nBins, double binWidth,
                                                    double minBin, double maxBin) {

        if (!MiscMath0.isAPowerOf2(nBins)) {
            throw new IllegalArgumentException("the number of histogram bins must be a power of 2");
        }

        System.out.printf("nBins=%d, binWidth=%.4f min=%.4f, max=%.4f\n", nBins, binWidth, minBin, maxBin);

        double[][] hist = new double[2][];
        hist[0] = new double[nBins];
        hist[1] = new double[nBins];

        for (int i = 0; i < nBins; i++) {
            hist[0][i] = minBin + i*binWidth + (binWidth/2.);
        }

        int bin;
        for (int i = 0; i < x.length; i++) {
            bin = (int) ((x[i] - minBin)/binWidth);
            if ((bin > -1) && (bin < nBins)) {
                hist[1][bin]++;
            }
        }

        // divide by x.length
        for (int i = 0; i < hist[1].length; ++i) {
            hist[1][i] /= x.length;
        }
        return hist;
    }

    /**
     *
     <pre>
     reference:
     Silverman, B. W. (1982). Algorithm as 176: Kernel density estimation using the fast Fourier transform. Journal of the Royal Statistical Society. Series C (Applied Statistics), 31(1), 93-99. https://dx.doi.org/10.2307/2347084.
     https://www.jstor.org/stable/2347084#metadata_info_tab_contents
     </pre>
     * @param u the FFT of the histogram of the data.
     * @param histBins the x bins of the histogram of the data that were used to create u
     * @return kernel density estimate
     */
    public static KDE viaFFTGaussKernel(Complex[] u, double[] histBins, double h) {

        // perform the fourier transform of the Gaussian Kernel K(h*s)
        //    FFT( K(h*s) ) = exp(-(0.5) * h^2 * s^2)
        // s = (x - xTilde[i])/h;  K(s*h)
        double[] fftKernel = new double[histBins.length];
        // note, this is not normalized.
        double zh;
        int i;
        for (i = 0; i < histBins.length; ++i) {
            zh = histBins[i] * h;
            if (zh < BIG) {
                fftKernel[i] = Math.exp(-0.5 * zh * zh);
            }
        }

        // FFT(f_n(s)) = exp(-(0.5) * h^2 * s^2) * u(s)   EQN (3)
        // Then f_n(s) = inverse FFT(f_n(s))
        //     negative values of f_n are set to 0.
        //     note that if several different uses of this algorithm for different h are employed,
        //     that the discrete FFT has to be calculated only once and can be reused.
        // The algorithm also avoids exponential underflow by setting
        // FFT(f_n(s)) equal to 0 if 0.5*h*h*s*s is > BIG which is large for the machine.

        // element-wise multiplication
        Complex[] eqn3 = new Complex[fftKernel.length];
        for (i = 0; i < fftKernel.length; ++i) {
            eqn3[i] = u[i].times(fftKernel[i]);
            //if (eqn3[i].re() < 0) {
            //    eqn3[i] = new Complex(0, 0);
            //}
        }

        FFT fft = new FFT();

        Complex[] kdeC = fft.ifft(eqn3);
        double[] kd = new double[kdeC.length];
        for (i = 0; i < kdeC.length; ++i) {
            kd[i] = kdeC[i].abs();
            //TODO: follow up on the math. presumably want magnitude instead of just the real component.
        }

        KDE kde = new KDE();
        kde.u = Arrays.copyOf(u, u.length);
        kde.h = h;
        kde.hx = Arrays.copyOf(histBins, histBins.length);
        kde.kde = Arrays.copyOf(kd, kd.length);

        return kde;
    }

    protected static boolean assertHistRange(double[] histBins, double[] dataMinMax, double h) {
        if (histBins.length < 3) {
            throw new IllegalArgumentException("histBins.length must be >= 3");
        }
        double dataRange = dataMinMax[1] - dataMinMax[0];
        double histBinSize = histBins[1] - histBins[0];
        double histRange = histBins[histBins.length - 1] - histBins[0] + histBinSize;

        // assert histRange + 6h >= dataRange
        return (histRange + 6.*h) >= dataRange;
    }

    /**
     NOTE: this method is replaced by the Silverman FFT approach.

    Wasserman chap 20:
    let X_0,...X_{N-1) denote the observed data which is a sample from f (= probability distribution).
    K is the kernel.
    the KDE estimator f_hat(x) = (1/n) * sum_i=0_to_{N-1} ( (1/h^d) * K((x - X_i)/h)
       where d is the dimensionality

     once the kde is estimated, that is, f_hat(x), the risk can be estimated with Wasserman eqn (20.24)
     J_hat(h) = integral( f_hat^2(x)*dz ) - (2/n) * summation_i=1_to_n( f_hat_{-i)(X_i) }
     where -i denotes "leave one out" of the sample.
    */

    /**
     estimate the KDE for a single value x.
     NOTE that this method is not efficent if used to calculate many x.  one should instead use the FFT method.
     <pre>
     Wasserman chap 20:
     let X_0,...X_{N-1) denote the observed data which is a sample from f (= probability distribution).
     K is the kernel.
     the1-Dimensional  KDE estimator f_hat(x) = (1/n) * sum_i=0_to_{N-1} ( (1/h) * K((x - X_i)/h).
     runtime complexity is O(xTilde.length).
     </pre>
     * @param kernel
     * @param x
     * @param xTilde the observed data
     * @param h the kernel bandwidth
     * @return the kde at x
     */
    public static double univariateKernelDensityEstimate(IKernel kernel, double x, double[] xTilde, double h) {
        return kernel.kernel(x, xTilde, h);
    }

}
