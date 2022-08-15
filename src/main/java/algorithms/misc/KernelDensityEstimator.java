package algorithms.misc;

import algorithms.imageProcessing.FFT;

/**
 * TODO: consider implementing the Improved Sheather - Jones algorithm for estimating bandwidth:.
 * can see implementations
 * https://github.com/tommyod/KDEpy/blob/master/KDEpy/bw_selection.py
 * and
 * https://kdepy.readthedocs.io/en/latest/_modules/KDEpy/FFTKDE.html#FFTKDE
 *
 * TODO: consider implementing Silverman Fast Fourier Transform (FFT): bin the data, map it to the frequency domain
 * using FFT, convolution, then inverse FFT.
 * The runtime complexity is O(n + m log m), with binning of n points followed by FFT calls on m-sized grids.
 * references:
 * https://github.com/tommyod/KDEpy/blob/master/KDEpy
 * and
 * "Fast & Accurate Gaussian Kernel Density Estimation", Jeffrey Heer, University of Washiington.
 * https://idl.cs.washington.edu/files/2021-FastKDE-VIS.pdf
 *
 * TODO: consider using MSER in a 2-D KDE estimator
 */
public class KernelDensityEstimator {

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
         * FFT of the data.
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
     * @param x
     * @param h
     * @return
     */
    public static KDE viaFFTGaussKernel(double[] x, double h) {

        //TODO: consider using PeriodicFFT from the curvature scale space project instead of FFT
        //     to avoid edge effects.

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
                  (follow Gentelman and Sande as imple. by Monro 1976 (? keeping to a power of 2?
                  the code below multiplies n by factor of 10 then takes the ceiling of the log base 2 of n...
                  math.ceil( math.log(n*10)/math.log(2) )
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

        throw new UnsupportedOperationException("not yet implemented");
    }

    /**
     * calculate a fine resolution histogram for x, for a larger data range than x's.
     * (NOTE that if the data are circular, a method can be created to calculate the histogram
     * using the same data range as x, not larger)
     * @param x
     * @return a 2-dimensional array of the histogram where hist[0] holds the centers of the histogram bins,
     * and hist[1] holds the counts within the bins.
     */
    protected static double[][] createFineHistogram(double[] x) {
        // the number of bins need to be a power of 2, and larger than x.length.
        int nBins = (int)Math.ceil(Math.log(x.length*2)/Math.log(2));

        // unless the data are circular, the range has to be larger than the range of x in order to
        // avoid wrap around edge conditions
        double[] minMaxX = MiscMath0.getMinMax(x);
        double range = minMaxX[1] - minMaxX[0];
        double min = minMaxX[0] - 0.25*range;
        double max = minMaxX[1] + 0.25*range;

        // nBins = (maxX - minX)/binWidth
        double binWidth = nBins/(max - min);

        double[][] hist = new double[2][];
        hist[0] = new double[nBins];
        hist[1] = new double[nBins];

        for (int i = 0; i < nBins; i++) {
            hist[0][i] = min + i*binWidth + (binWidth/2.);
        }

        int bin;
        for (int i = 0; i < x.length; i++) {
            bin = (int) ((x[i] - min)/binWidth);
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
     * @param u
     * @param h
     * @return
     */
    public static double[] viaFFTGaussKernel(Complex[] u, double h) {
        throw new UnsupportedOperationException("not yet implemented");
    }

    /**
    Wasserman chap 20:
    let X_0,...X_{N-1) denote the observed data which is a sample from f (= probability distribution).
    K is the kernel.
    the KDE estimator f_hat(x) = (1/n) * sum_i=0_to_{N-1} ( (1/h^d) * K((x - X_i)/h)
       where d is the dimensionality

    public static double fEst(IKernel kernel, double[][] x, double[] h) {
        throw new UnsupportedOperationException("not yet implemented");
    }*/

    /**
     Wasserman chap 20:
     let X_0,...X_{N-1) denote the observed data which is a sample from f (= probability distribution).
     K is the kernel.
     the1-Dimensional  KDE estimator f_hat(x) = (1/n) * sum_i=0_to_{N-1} ( (1/h) * K((x - X_i)/h).
     runtime complexity is O(xTilde.length).
     */
    public static double univariateKernelDensityEstimate(IKernel kernel, double x, double[] xTilde, double h) {
        return kernel.kernel(x, xTilde, h);
    }

    /*
    using FFT to create KDE/PDA:

    https://github.com/TasCL/cpda/

    refine these:

    construct a simulated histogram,
    use FFT to transform the histogram to spectral domain,
    applies standard Gaussian kernels to smooth it,
    uses inverse FFT to get simulated PDF
     */
}
