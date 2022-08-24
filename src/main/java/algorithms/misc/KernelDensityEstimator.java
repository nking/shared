package algorithms.misc;

import algorithms.imageProcessing.FFTUtil;
import algorithms.matrix.MatrixUtil;

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
     * NOTE, the method ruleOfThumbBandwidthGaussianKernel() is preferred.
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
     * estimate the density by using properties of exponentials to make a faster algorithm.
     <pre>
     reference:
     Silverman, B. W. (1982). Algorithm as 176: Kernel density estimation using the fast Fourier transform. Journal of the Royal Statistical Society. Series C (Applied Statistics), 31(1), 93-99. https://dx.doi.org/10.2307/2347084.
     https://www.jstor.org/stable/2347084#metadata_info_tab_contents
     </pre>
     calculate the KDE in a fast manner by using discrete FFTs and using a grid of data points
     the kde = f_hat(x) ~ summation over i=1 to n of (K(|| x - X_i ||/h).
     the kernel estimate is a convolution of the data with the kernel.
     naive implementation is O(n^2).
     if the kernel K is chosen to be a Gaussian, one can use a property of exponentials to rewrite:
     exp(x - X_i) = exp(x) * exp(-X_i) (neglecting details)
     also note that the discrete FFT is a summation of exponentials.
     convolution theorem: one can use the elementwise multiplication between the fourier paired functions.
     (Chap 15.5, Boas "Mathematical Methods in the Physical Sciences")
     convolution to FFT:
     FFT(kde) ~ FFT(K(h*s)) * FFT( hist(X) )
     one can choose s to be the same spatial intervals (=grid) in the histogram
     and in the kernel, to avoid interpolation.  the multiplication is element-wise.
     then inverse FFT of FFT(kde) = kde.
     the runtime complexity is then O(n_s*log(n_s))
     * @param x observed data
     * @param h bandwidth
     * @return
     */
    public static KDE viaFFTGaussKernel(double[] x, double h, int histNBins, double histBinWidth, double histMinBin,
                                        double histMaxBin) {

        // the histogram length (hist[1].length) is a power of 2 and the range is enlarged by at least 2*3*h
        double[][] hist = createFineHistogram(x, h, histNBins, histBinWidth, histMinBin, histMaxBin);

        //assert(assertHistRange(hist[0], MiscMath0.getMinMax(x), h) == true);

        // normalization is performed by default:
        FFTUtil fft = new FFTUtil();

        // u is the portion that can be re-used on subsequent iterations.  e.g. when iterating
        //   to find minimum bandwidth h.
        Complex[] u = fft.create1DFFTNormalized(hist[1], true);

        assert(u.length == hist[0].length);

        return viaFFTGaussKernel(u, hist[0], h);
    }

    /**
     *
     * estimate the density by using properties of exponentials to make a faster algorithm.
     <pre>
     reference:
     Silverman, B. W. (1982). Algorithm as 176: Kernel density estimation using the fast Fourier transform. Journal of the Royal Statistical Society. Series C (Applied Statistics), 31(1), 93-99. https://dx.doi.org/10.2307/2347084.
     https://www.jstor.org/stable/2347084#metadata_info_tab_contents
     </pre>
     calculate the KDE in a fast manner by using discrete FFTs and using a grid of data points
     the kde = f_hat(x) ~ summation over i=1 to n of (K(|| x - X_i ||/h).
     the kernel estimate is a convolution of the data with the kernel.
     naive implementation is O(n^2).
     if the kernel K is chosen to be a Gaussian, one can use a property of exponentials to rewrite:
     exp(x - X_i) = exp(x) * exp(-X_i) (neglecting details)
     also note that the discrete FFT is a summation of exponentials.
     convolution theorem: one can use the elementwise multiplication between the fourier paired functions.
     (Chap 15.5, Boas "Mathematical Methods in the Physical Sciences")
     convolution to FFT:
     FFT(kde) ~ FFT(K(h*s)) * FFT( hist(X) )
     one can choose s to be the same spatial intervals (=grid) in the histogram
     and in the kernel, to avoid interpolation.  the multiplication is element-wise.
     then inverse FFT of FFT(kde) = kde.
     the runtime complexity is then O(n_s*log(n_s))
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

        // normalization is performed by default:
        FFTUtil fft = new FFTUtil();

        // u is the portion that can be re-used on subsequent iterations.  e.g. when iterating
        //   to find minimum bandwidth h.
        Complex[] u = fft.create1DFFTNormalized(hist[1], true);

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

    /**
     * estimate the density of x using a Gaussian Kernel of bandwidth h.
     * note that the grid to which the kernel is applied is also x in this case.
     <pre>
     f_hat(x) ~ (1/n) * summation over i=1 to n of ( (1/h) * K(|| x - X_i ||/h) ).
     double[] kde = iter from j=1-1 to xGrid.length-1 estimating f_hat(xGrid[j]
     reference: Wasserman's "All of Statistics", eqn (20.21)
     </pre>
     runtime complexity is O(n^2) where n = x.length.
     * @param x data observed
     * @param h bandwidth
     * @return
     */
    public static KDE viaGaussKernel(double[] x, double h) {
        int n = x.length;
        int i;
        KDE kde = new KDE();
        kde.kde = new double[n];
        GaussianKernel kernel = new GaussianKernel();
        double sum = 0;
        for (i = 0; i < n; ++i) {
            kde.kde[i] = kernel.kernel(x[i], x, h);
            sum += kde.kde[i];
        }
        System.out.printf("SUM KDE=%.4e\n", sum);
        kde.hx = Arrays.copyOf(x, x.length);
        kde.h = h;
        kde.u = null;
        return kde;
    }

    /**
     * estimate the density of x using a Gaussian Kernel of bandwidth h.
     * the grid to which the kernel is applied is xGrid in this case.
     <pre>
     f_hat(x) ~ (1/n) * summation over i=1 to n of ( (1/h) * K(|| x - X_i ||/h) ).
     double[] kde = iter from j=1-1 to xGrid.length-1 estimating f_hat(xGrid[j]
     reference: Wasserman's "All of Statistics", eqn (20.21)
     </pre>
     runtime complexity is O(n*m) where n = x.length and m=xGrid.length
     * @param x data observed
     * @param xGrid the grid of data at which to apply the kernel.
     * @param h bandwidth
     * @return
     */
    public static KDE viaGaussKernel(double[] x, double[] xGrid, double h) {
        int n = xGrid.length;
        int i;
        KDE kde = new KDE();
        kde.kde = new double[n];
        GaussianKernel kernel = new GaussianKernel();
        double sum = 0;
        for (i = 0; i < n; ++i) {
            kde.kde[i] = kernel.kernel(xGrid[i], x, h);
            sum += kde.kde[i];
        }
        System.out.printf("SUM KDE=%.4e\n", sum);
        kde.hx = Arrays.copyOf(xGrid, xGrid.length);
        kde.h = h;
        kde.u = null;

        return kde;
    }

    protected static Complex[] convertToComplex(double[] a) {
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
    public static double[][] createFineHistogram(double[] x, double h) {
        // the number of bins need to be a power of 2, and larger than x.length.
        int nBins;
        //if (x.length < 50) {
        //    nBins = x.length;
        //} else {
            nBins = (int) Math.pow(2, Math.ceil(Math.log(x.length * 3) / Math.log(2)));
        //}

        // the power of 10 was inspired by method vbwkde, variable n_dct from documentation:
        //https://user-web.icecube.wisc.edu/~peller/pisa_docs/_modules/pisa/utils/vbwkde.html

        //TODO: improve this with an estimate for available memory of machine
        if (nBins > 16384) {
            nBins = 16384;
        }

        // unless the data are circular, the range has to be larger than the range of x in order to
        // avoid wrap around edge conditions
        double[] minMaxX = MiscMath0.getMinMax(x);
        double range = minMaxX[1] - minMaxX[0];
        double dr;
        double he = 3.*h;
        if (0.5 * range > he) {
            dr = 0.5 * range;
        } else {
            dr = he;
        }

        double min = minMaxX[0] - dr;
        double max = minMaxX[1] + dr;

        // nBins = (maxX - minX)/binWidth
        double binWidth = (max - min)/nBins;

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

        System.out.printf("nBins=%d, binWidth=%.4f min=%.4f, max=%.4f\n", nBins, binWidth, minBin, maxBin);
        System.out.flush();

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
            hist[1][i] /= ((double)x.length);
        }
        return hist;
    }

    /**
     * calculate a substitute for the fine resolution histogram for x for use in the
     * risk estimator that uses cross-validation.
     * @param x
     * @return a 2-dimensional array of the histogram where hist[0] holds the centers of the histogram bins,
     * and hist[1] holds the counts within the bins.
     */
    protected static double[] createFineHistogramSubstitute(double[] x) {
        double[] yHist = new double[x.length];
        Arrays.fill(yHist, 1.);
        // divide by x.length
        for (int i = 0; i < yHist.length; ++i) {
            yHist[i] /= x.length;
        }
        return yHist;
    }

    /**
     estimate the density by using properties of exponentials to make a faster algorithm.
     <pre>
           reference:
           Silverman, B. W. (1982). Algorithm as 176: Kernel density estimation using the fast Fourier transform. Journal of the Royal Statistical Society. Series C (Applied Statistics), 31(1), 93-99. https://dx.doi.org/10.2307/2347084.
           https://www.jstor.org/stable/2347084#metadata_info_tab_contents
     </pre>
     calculate the KDE in a fast manner by using discrete FFTs and using a grid of data points
     the kde = f_hat(x) ~ summation over i=1 to n of (K(|| x - X_i ||/h).
     the kernel estimate is a convolution of the data with the kernel.
     A naive implementation such as the method viaGaussKernel() is O(n^2).
     If the kernel K is chosen to be a Gaussian, one can use a property of exponentials to rewrite:
     exp(x - X_i) = exp(x) * exp(-X_i) (neglecting details)
     also note that the discrete FFT is a summation of exponentials.
     convolution theorem: one can use the element-wise multiplication between the fourier paired functions.
     (Chap 15.5, Boas "Mathematical Methods in the Physical Sciences")
     convolution to FFT:
     FFT(kde) ~ FFT(K(h*s)) * FFT( hist(X) )
     one can choose s to be the same spatial intervals (=grid) in the histogram
     and in the kernel, to avoid interpolation.  the multiplication is element-wise.
     then inverse FFT of FFT(kde) = kde.
     the runtime complexity is then O(n_s*log(n_s))
     * The Gaussian Kernel with discrete fast fourier transforms is O(s*log(s)).
     * @param u the FFT of the histogram of the data.
     * @param histBins the x bins of the histogram of the data that were used to create u
     * @return kernel density estimate
     */
    public static KDE viaFFTGaussKernel(Complex[] u, double[] histBins, double h) {

        //K(x) ≥ 0, ∫K(x)dx = 1, ∫xK(x)dx = 0

        // perform the fourier transform of the Gaussian Kernel K(h*s)
        //    FFT( K(h*s) ) = exp(-(0.5) * h^2 * s^2)
        // s = (x - xTilde[i])/h;  K(s*h)

        // FFT(f_n(s)) = exp(-(0.5) * h^2 * s^2) * u(s)   EQN (3)
        // Then f_n(s) = inverse FFT(f_n(s))
        //     negative values of f_n are set to 0.
        //     note that if different uses of this algorithm for different h are employed,
        //     that the discrete FFT has to be calculated only once and can be reused.
        // The algorithm also avoids exponential underflow by setting
        // FFT(f_n(s)) equal to 0 if 0.5*h*h*s*s is > BIG which is large for the machine.
        int i;
        double zh;
        double c = 1./Math.sqrt(2.*Math.PI);
        // element-wise multiplication
        Complex[] eqn3 = new Complex[histBins.length];
        for (i = 0; i < histBins.length; ++i) {
            zh = histBins[i] * h;
            if (zh < BIG) {
                eqn3[i] = u[i].times(c * Math.exp(-0.5 * zh * zh));
            }
        }

        FFTUtil fft = new FFTUtil();

        Complex[] kdeC = fft.create1DFFTNormalized(eqn3, false);
        double[] kd = new double[kdeC.length];
        //double d = 1./(kd.length - 1.);
        for (i = 0; i < kdeC.length; ++i) {
            kd[i] = kdeC[i].abs();
        }

        //
        double sum = sumKDE(kd, histBins[1] - histBins[0]);
        System.out.printf("h=%.4e sumKDE=%.4e\n", h,  sum);

        KDE kde = new KDE();
        kde.u = Arrays.copyOf(u, u.length);
        kde.h = h;
        kde.hx = Arrays.copyOf(histBins, histBins.length);
        kde.kde = Arrays.copyOf(kd, kd.length);

        return kde;
    }

    static double sumKDE(double[] kde, double dk) {
        double sum = 0;
        for (int i = 0; i < kde.length; ++i) {
            sum += kde[i];
        }
        return sum;
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

    /**
     * a fast cross-validation method which can be used in choosing the bandwidth for the kernel density estimator.
     * the runtime complexity is max( O(n), O(n_s*log_2(ns)))
     *  where n is the data length, and n_s is the number of elements in the data histogram
     * (which is equal to u.length).
     * This method is implemented for 1-D only, but the references (eqn 27) provide a formula for multiple dimensions.
     <pre>
     References:
     Wasserman, "All of Statistics", eqn (20.25)
     and
     https://www.stat.cmu.edu/~larry/=sml/densityestimation.pdf
     36-708 Statistical Methods for Machine Learning by Larry Wasserman, CMU
     eqn (27) (which is (27)+(28))
     and
     Indirect Cross-validation for Density Estimation
     Olga Y. Savchuk, Jeffrey D. Hart, Simon J. Sheather
     https://doi.org/10.48550/arxiv.0812.0051,
     https://arxiv.org/abs/0812.0051,
     eqn (2).
     </pre>
     * @param u fft of the histogram of the data.
     * @param histBins the x axis of the histogram of the data.
     * @param h the bandwidth to use
     * @return the cross validation score
     */
    public static double crossValidationScore(Complex[] u, double[] histBins, double h) {

        /*
        let K(z, sigma) = normal (gaussian) kernel with mean 0 and variance sigma^2.

        r_hat(h) = (K(0, sqrt(2)*h)/(n-1))
                   + ((n-2)/(n*((n-1)^2))) * summation over i where i!=j of (K(Xi-Xj, sqrt(2)*h))
                   - (2/(n*(n-1))) * summation over i where i!=j of (K(Xi-Xj, h))

        The 2nd and 3rd terms can be estimated using the properties of convolution and discrete FFTs,
        similar to the way the kernel density is estimated using FFTs.
        */

        int n = histBins.length;

        // perform the fourier transform of the Gaussian Kernel K(h*s)
        //    FFT( K(h*s) ) = exp(-(0.5) * h^2 * s^2)
        // s = (x - xTilde[i])/h;  K(s*h)  <== multiply h differently for term1, term2, term3
        double sh = Math.sqrt(2) * h;
        double term1 = 0;

        int i;
        double zh;
        double m;
        // term2: ((n-2)/(n*((n-1)^2))) * summation over i where i!=j of (K(Xi-Xj, sqrt(2)*h))
        //term3: (2/(n*(n-1))) * summation over i where i!=j of (K(Xi-Xj, h))
        double[] f2 = new double[n];
        double[] f3 = new double[n];
        for (i = 0; i < n; ++i) {
            zh = histBins[i] / sh;
            if (zh < BIG) {
                m = -0.5 * zh * zh;
                term1 += Math.exp(m);
                f2[i] = u[i].times(Math.exp(m) ).abs();
                // if (m >= 0) {
            }
            zh = histBins[i] / h;
            if (zh < BIG) {
                f3[i] = u[i].times(Math.exp(-0.5 * zh * zh)).abs();
                // if (m >= 0) {
            }
        }
        term1 /= (n*(n - 1.));

        FFTUtil fft = new FFTUtil();
        Complex[] t2 = fft.create1DFFTNormalized(f2, false);
        Complex[] t3 = fft.create1DFFTNormalized(f3, false);

        double term2 = 0;
        double term3 = 0.;
        for (i = 0; i < n; ++i) {
            term2 += t2[i].abs();
            term3 += t3[i].abs();
        }

        term2 *= ((n-2.)/(n*(n-1.)*(n-1.)));
        term3 *= (2./(n*(n-1.)));

        // there is an error in my implementation
        double r = term1 + term2 - term3;

        // fudge, to be removed when the bug is found
        r = term1 + (term2 - term3)*(n-1);

        System.out.printf("h=%.4f terms=%.4e %.4e %.4e  r=%.4e\n", h, term1, term2, term3, r);

        return r;

        /*
        try again, but with the Wasserman eqn (20.25)

        r_hat(h) ~ (2/(n*h)) * K(0)
                    + (1/(h*n*n)) * sum_i( sum_j(
                         K_ast((X_i-X_j)/h) ))

        where
        K_2(z) = integral( K(z-y)*K(y)*dy )
        K_ast(x) = K_2(x) - 2*K(x)

        z = (X_i-X_j)/h;
        K_ast((X_i-X_j)/h) = K_2(z) - 2*K(z) = N(0,2) - 2*K(z)
                           = integral( K(z-y)*K(y)*dy ) - 2*K(z)

        r_hat(h) ~ (2/(n*h)) * N(0,1)
                 + (1/(h*n*n)) * sum_i( sum_j(
                     integral( K(z-y)*K(y)*dy ) - 2*K(z)
                 ))
        */
    }

    /**
     * calculate the risk estimator for the "leave-one-out" method of cross validation.
     * This method can be used to select a kde bandwidth by minimizing the risk estimator for a range
     * of bandwidths.  It is a data based method as the true probability distribution is unknown.
     <pre>
     see eqn (20.24) of Wasserman's "All of Statistics"
     see eqn (26) of
     https://www.stat.cmu.edu/~larry/=sml/densityestimation.pdf
     36-708 Statistical Methods for Machine Learning by Larry Wasserman, CMU
     </pre>
     * @return
     */
    static double riskEstimatorLeaveOneOut() {

        // r_hat(h) = integral( (p_hat(x))^2*dx - (2/n) * summation_i=1_to_n( p_hat(X_i) )
        //    the 3rd compoonent on the right hand side is dropped as it's a constant (and so cancels out in comparisons
        //    for bandwidth selection) and it is nearly negligible as n increases.
        //
        //     for the "leave-one-out" algorithm, p_hat is the density estimator after removing
        //        the "i-th" observation.
        //     for the "data-splitting" algorithm,
        //        X_i is split into half randomly, and that instead of the entire X_i is used
        //        in the summation above.
        //        Note that splitting in half is V-fold or k-fold of 2, but a larger number can be used instead.

        throw new UnsupportedOperationException("not yet implemented");
    }

    /**
     * calculate the risk estimator for the "data-splitting" V-fold or k-fold method of cross validation.
     * This method can be used to select a kde bandwidth by minimizing the risk estimator for a range
     * of bandwidths.  It is a data based method as the true probability distribution is unknown.
     <pre>
     see eqn (20.24) of Wasserman's "All of Statistics".
     see eqn (26) of
     https://www.stat.cmu.edu/~larry/=sml/densityestimation.pdf
     36-708 Statistical Methods for Machine Learning by Larry Wasserman, CMU
     </pre>
     * @return
     */
    static double riskEstimatorDataSplitting() {

        // r_hat(h) = integral( (p_hat(x))^2*dx - (2/n) * summation_i=1_to_n( p_hat(X_i) )
        //    the 3rd compoonent on the right hand side is dropped as it's a constant (and so cancels out in comparisons
        //    for bandwidth selection) and it is nearly negligible as n increases.
        //
        //     for the "leave-one-out" algorithm, p_hat is the density estimator after removing
        //        the "i-th" observation.
        //     for the "data-splitting" algorithm,
        //        X_i is split into half randomly, and that instead of the entire X_i is used
        //        in the summation above.
        //        Note that splitting in half is V-fold or k-fold of 2, but a larger number can be used instead.

        throw new UnsupportedOperationException("not yet implemented");
    }

    /**
     * calculate the risk estimator for the "leave-one-out" method of cross validation.
     * This method can be used to select a kde bandwidth by minimizing the risk estimator for a range
     * of bandwidths.  It is a data based method as the true probability distribution is unknown.
     <pre>
     see eqn (27), (28) of
     https://www.stat.cmu.edu/~larry/=sml/densityestimation.pdf
     36-708 Statistical Methods for Machine Learning by Larry Wasserman, CMU.
     see eqn (20.24) of Wasserman's "All of Statistics"
     </pre>
     * @return
     */
    static double riskEstimatorGaussianKernel() {

        // r_hat(h) = integral( (p_hat(x))^2*dx - (2/n) * summation_i=1_to_n( p_hat(X_i) )
        //    the 3rd compoonent on the right hand side is dropped as it's a constant (and so cancels out in comparisons
        //    for bandwidth selection) and it is nearly negligible as n increases.
        //
        //     for the "leave-one-out" algorithm, p_hat is the density estimator after removing
        //        the "i-th" observation.
        //     for the "data-splitting" algorithm,
        //        X_i is split into half randomly, and that instead of the entire X_i is used
        //        in the summation above.
        //        Note that splitting in half is V-fold or k-fold of 2, but a larger number can be used instead.

        throw new UnsupportedOperationException("not yet implemented");
    }

    static double sumHistogram(double[][] hist) {
        double sum = 0;
        double dh = hist[0][1] - hist[0][0];
        for (int i = 0; i < hist[0].length; ++i) {
            sum += (hist[1][i]);
        }
        return sum;
    }

}
