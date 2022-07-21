package algorithms.statistics;

import algorithms.matrix.MatrixUtil;
import no.uib.cipr.matrix.NotConvergedException;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Random;

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
 * For Gumbel distribution, k = 0, sigma > 0
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
public class GumbelCDF {

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
        double eps = 1e-19;
        for (i = 0; i < nDraws; ++i) {
            u = rand.nextDouble();
            while (u < eps) {
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
}
