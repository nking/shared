package algorithms.statistics;

import algorithms.misc.MiscMath0;

import java.security.SecureRandom;

/**
 *  the Gumbel distribution (Generalized Extreme Value distribution Type-I)
 *  is used to model the distribution of the maximum (or the minimum) of a number
 *  of samples of various distributions.
 *
 * <pre>
 * parameters:
 *    mu is the location
 *    sigma is the scale
 *    k is the shape
 *
 * For Gumbel distribution, k = 0, sigma .gt. 0
 *
 * Let z = (x-mu)/sigma
 *
 *  PDF = f(x) = (1/sigma)*exp(-z)*exp(-exp(-z) )
 *  CDF = F(x) = exp( -exp(-z) )
 *  inverse of CDF = G(alpha) = mu - sigma*log( log(1/alpha)
 *         where alpha is the random variate drawn from U(0,1))
 *
 * references:
 *     "Statistical Distibutions", 2011, ed 4
 *         Merran Evans, Nicholas Hastings, Brian Peacock, and Catherine Forbes
 *
 *      https://www.statisticshowto.com/gumbel-distribution/
 *
 *     https://en.m.wikipedia.org/wiki/Gumbel_distribution
 * </pre>
 */
public class Gumbel {

    /**
     * generate a point in the Gumbel cumulative distribution
     * @param x a general element, that is a quantile, of the larger range of general variate X
     * @param location parameter of the distribution function
     * @param scale parameter of the distribution function
     * @return
     */
    public static double cdf(double x, double location, double scale) {
        double z = (x - location) / scale;
        return Math.exp(-Math.exp(-z));
    }

    /**
     * calculate the probability of a point in the Gumbel discrete probability density function
     * This uses the minimum Gumbel distribution and has a right leaning distribution.
     * @param x a general element, that is a quantile, of the larger range of general variate X
     * @param location parameter of the distribution function
     * @param scale parameter of the distribution function
     * @return
     */
    public static double pdfForMinimum(double x, double location, double scale) {
        double z = (x - location) / scale;
        return (1. / scale) * Math.exp(z) * Math.exp(-Math.exp(z));
    }

    /**
     * calculate the probability of a point in the Gumbel discrete probability density function
     * This uses the maximum Gumbel distribution and has a left leaning distribution.
     * @param x a general element, that is a quantile, of the larger range of general variate X
     * @param location parameter of the distribution function
     * @param scale parameter of the distribution function
     * @return
     */
    public static double pdf(double x, double location, double scale) {
        double z = (x - location) / scale;
        return (1. / scale) * Math.exp(-z) * Math.exp(-Math.exp(-z));
    }

    /**
     * calculate the inverse CDF of the Gumbel, that is, a random variate x given the
     * probability alpha.
     * @param alpha random variate drawn from U(0,1) where U is the uniform distribution.
     * @param location parameter of the distribution function
     * @param scale parameter of the distribution function
     * @return
     */
    public static double inverseCdf(double alpha, double location, double scale) {
        return location - scale*Math.log(-Math.log(alpha));
    }

    /**
     * sample from a Gumbel distribution G(location, scale).
     * @param location
     * @param scale
     * @param nDraws number of random draws to make
     * @return a fair sampling from a Gumbel distribution G(location, scale).
     */
    public static double[] sampleRandomlyFrom(double location, double scale,
        int nDraws, SecureRandom rand) {

        double[] out = new double[nDraws];
        int i;
        double u;
        double eps = 1e-320;
        for (i = 0; i < nDraws; ++i) {
            u = rand.nextDouble();
            while (u < eps) {
                // to stay within math domain of Math.log(u)
                u = rand.nextDouble();
            }
            out[i] = inverseCdf(u, location, scale);
        }

        return out;
    }

    /**
     * @param x ordered statistic of an observed Gumbel distribution.
     */
    public static double[] empiricalCdf(double[] x) {
        throw new UnsupportedOperationException("not yet implemented");
    }

    public static double[] generateCurve(double[] x, double location, double scale) {
        double[] y = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            y[i] = pdf(x[i], location, scale);
        }
        return y;
    }

    /**
     * calculate a rough estimate of Gumbel distribution parameters for the given x.
     * A more precise estimate can be obtained from fitGumbelUsingML or fitGumbelUsingBayesian
     * when they are implemented.
     * <pre>
     * references:
     *     Chap 19 of "Statistical Distributions" by Evans et al.
     *     and
     *     https://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm
     * </pre>
     * @param x ordered statistic of an observed Gumbel distribution.
     * @return
     */
    public static double[] fitUsingMethodOfMoments(double[] x) {

        int n = x.length;

        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(x);
        double kMAD = 1.4826;
        double s = kMAD*mADMinMax[0];

        double[] meanStdv = MiscMath0.getAvgAndStDev(x);

        double sigma0 = meanStdv[1] * Math.sqrt(6.)/Math.PI;
        double sigma1 = s * Math.sqrt(6.)/Math.PI;
        double mu0 = meanStdv[0] - sigma0 * MiscMath0.eulerMascheroniConstant();
        double mu1 = meanStdv[0] - sigma1 * MiscMath0.eulerMascheroniConstant();

        //double r0 = mADMinMax[1] - 3*s;
        //double r1 = mADMinMax[1] + 3*s;

        System.out.println("mu0,1=" + mu0 + ", " + mu1);
        System.out.println("sigma0,1=" + sigma0 + ", " + sigma1);

        return new double[]{mu1, sigma1, 0};
    }

    /**
     * estimate the parameters mu and sigma (location and scale, respectively) using
     * method of maximum likelihood simultaneous solution.
     * <pre>
     *    reference is Chap 19 of "Statistical Distributions" by Evans et al.
     *    and
     *    https://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm
     * </pre>
     * @param x ordered values for which to find the best fitting Gumbel distribution parameters
     * @return
     */
    public static double[] fitUsingMaximumLikelihood(double[] x) {

        double[] params = fitUsingMethodOfMoments(x);
        double[] avgAndStdev = MiscMath0.getAvgAndStDev(x);
        /*
        method of maximum likelihood simultaneous solutions:
        scaleEst = x_avg - ( sum_over_i(x_i * exp(-x_i/scaleEst) ) / sum_over_i(exp(-x_i/scaleEst) ) )
        locEst = -scaleEst * Math.log( (1/n) * sum_over_i( exp(-x_i/scaleEst) ))

        initial estimates from method of moments

        assert that locEst < x_avg (due to skew)
        */

        double locEst;
        double scaleEst;

        // simple loop until scale estimate change is small is not efficient
        // so prefer the optimization
        boolean useOptimization = false;
        if (!useOptimization) {
            final int nIterMax = 100;
            int nIter = 0;
            // TODO: allow tolerance and nIterMax to be method arguments if keep the simple loop iteration
            final double tol = 1e-3;
            double prevScaleEst = params[1];
            double diffScale;
            do {
                scaleEst = estimateScaleML(x, avgAndStdev[0], prevScaleEst);
                locEst = estimateLocML(x, scaleEst);
                diffScale = scaleEst - prevScaleEst;
                System.out.printf("locEst=%11.6e, prevScaleEst=%11.6e, scaleEst=%11.6e, diffScale=%11.6e\n",
                        locEst, prevScaleEst, scaleEst, diffScale);
                prevScaleEst = scaleEst;
                nIter++;
            } while ((nIter < nIterMax) && Math.abs(diffScale) > tol);
            System.out.printf("nIter=%d\n", nIter);
        } else {
            // non-linear optimization.
            // solving for 2 parameter Gumbel distribution should be convex, so one could use the
            //   first derivative, the second derivative, and or the gradient.

            throw new UnsupportedOperationException("not yet implemented");
        }
        return new double[]{locEst, scaleEst};
    }

    private static double estimateLocML(double[] x, double scaleEst) {
        int n = x.length;
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += Math.exp(-x[i]/scaleEst);
        }
        return -scaleEst * Math.log((1./n)*sum);
    }

    private static double estimateScaleML(double[] x, double xAvg, double scaleEst) {
        int n = x.length;
        double sum0 = 0;
        double sum1 = 0;
        for (int i = 0; i < n; ++i) {
            sum0 += (x[i] * Math.exp(-x[i]/scaleEst));
            sum1 += (Math.exp(-x[i]/scaleEst));
        }
        scaleEst = xAvg - ( sum0 / sum1 );
        return scaleEst;
    }

    /**
     * generate the generalized Gumbel probability density curve (GEV Type I).
     * <pre>
     *     reference is Chap 19 of "Statistical Distributions" by Evans et al.
     * </pre>
     * @param x1
     * @param mu
     * @param sigma
     * @return
     */
    public static double[] generateGumbelCurve(double[] x1, double mu, double sigma) {
        if (sigma <= 0) {
            throw new IllegalArgumentException("sigma must be > 0");
        }
        double[] yGEV = new double[x1.length];
        double z;
        double a;
        for (int i = 0; i < x1.length; i++) {
            z = (x1[i] - mu)/sigma;
            a = Math.exp(-z);
            yGEV[i] = (1./sigma) * Math.exp(-z) * Math.exp(-a);
            //yGEV[i] = (1./sigma) * Math.exp(-(Math.exp(-z) + z)); // same result
        }
        return yGEV;
    }
}
