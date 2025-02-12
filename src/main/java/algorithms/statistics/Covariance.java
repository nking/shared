package algorithms.statistics;

import algorithms.matrix.MatrixUtil;

import java.util.Arrays;

/**
 * class holding methods for covariance and correlations.
 *
 * NOTE: see MultivariateDistance for fast distance covariance.
 *
 * @author nichole
 */
public class Covariance {

    /**
     * enumeration of data preparation.
     * @param UNCENTERED: the data have had no transformations
     * @param MEAN_CENTERED : the data have had their means subtracted
     * @param UNIT_STAND_MEAN0_STD1 the data have had their means subtracted then were divided
     *                              by their standard deviations
     */
    public static enum STAND_TYPE {
        UNCENTERED,
        MEAN_CENTERED,
        UNIT_STAND_MEAN0_STD1
    }

    /**
     * calculating Pearson correlation as
     * (E(X*Y) - E(X)*E(Y))/ sqrt((E[X^2] - E[X]^2)*(E[Y^2] - E[Y]^2))
     * @param x 1st dimension points of xy dataset.  the points do not need to be mean centered.
     * @param y 2nd dimension points of xy dataset.  the points do not need to be mean centered.
     * @return the Pearson correlation.
     */
    public static double correlationPearson(double[] x, double[] y) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }
        int n = x.length;

        double[] moments = Util.caldc2DMomentsX2Y2(x, y);

        // divide by n to make them averages:
        for (int i = 0; i < moments.length; ++i) {
            moments[i] /= n;
        }

        /*
        = (avg( moment[xy]) - avg(moment[x])*avg(moment[y]))
           / (sqrt(
              avg(moment[x^2] - (moment[x]*moment[x]))
              * avg(moment[y^2] - (moment[y]*moment[y]) )
             ))
        */
        //x    y    x^2    y^2    xy   = moments
        //0    1    2      3       4   = indexes in output

        double rxy = (moments[4] - moments[0] * moments[1]) /
                (Math.sqrt(moments[2] - moments[0] * moments[0])
                * Math.sqrt(moments[3] - moments[1] * moments[1]) );

        return rxy;
    }

    /**
     * calculate the sample correlation of the xy dataset given the type of standardiation that
     * has been performed on them.
     @param x 1st dimension points of xy dataset.  the points do not need to be mean centered.
     @param y 2nd dimension points of xy dataset.  the points do not need to be mean centered.
     @return the Pearson correlation.
     @param type type of dataprocessing already performed on the data.  if type is null, it is
     assigned UNCENTERED.
     * @return the sample correlation
     */
    public static double correlationSample(double[] x, double[] y, STAND_TYPE type) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }
        int n = x.length;

        if (type == null) {
            type = STAND_TYPE.UNCENTERED;
        }

        //x    y    x^2    y^2    xy   = moments
        //0    1    2      3       4   = indexes in output
        double[] moments = Util.caldc2DMomentsX2Y2(x, y);

        switch (type) {
            case UNCENTERED: {
                return  (n * moments[4] - moments[0] * moments[1])
                        / (
                        Math.sqrt(n*moments[2] - moments[0] * moments[0])
                        * Math.sqrt(n*moments[3] - moments[1] * moments[1])
                );
            }
            case MEAN_CENTERED: {
                return moments[4]/(Math.sqrt(moments[2]) * Math.sqrt(moments[3]));
            }
            case UNIT_STAND_MEAN0_STD1: {
                return (1/(n-1)) * (moments[0] * moments[1]);
            }
        }
        throw new IllegalArgumentException("error in alg");
    }

    /**
     * calculate the sample covariance of the xy dataset given the type of standardiation that
     * has been performed on them.
     @param x 1st dimension points of xy dataset.  the points do not need to be mean centered.
     @param y 2nd dimension points of xy dataset.  the points do not need to be mean centered.
     @return the sample covariance.
     @param type type of dataprocessing already performed on the data.  if type is null, it is
     assigned UNCENTERED.
      * @return the sample covariance.
     <pre>
     if type==UNCENTERED, it calculates the covarianve using brute force method .
     if type==MEAN_CENTERED, uses the sum of moments to calculate result quickly.
    if type==UNIT_STAND_MEAN0_STD1, uses the sum of mements, but the returned value is
    actually the correlation.  To transform that number to covariance with respect to the
    original dataset reference frame, multiply this result by the standard deviation of
    x and the standard deviation of x.
    </pre>
     */
    public static double covarianceSample(double[] x, double[] y, STAND_TYPE type) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }
        int n = x.length;

        if (type == null) {
            type = STAND_TYPE.UNCENTERED;
        }

        double[] meanX = new double[1];
        double[] meanY = new double[1];
        double[] stdX = new double[]{1};
        double[] stdY = new double[]{1};

        switch (type) {
            case UNCENTERED: {
                double[][] a = new double[n][];
                for (int i = 0; i < n; ++i) {
                    a[i] = new double[]{x[i], y[i]};
                }
                double[][] cov = BruteForce.covariance(a);
                return cov[0][1];
            }
            case MEAN_CENTERED: {
                double[] moments = Util.caldc2DMomentsX2Y2(x, y);
                //    Cov(X, Y) = (1 / (n - 1)) * Σ ( (Xi - mean(X)) / std(X) ) * ( (Yi - mean(Y)) / std(Y) )
                double cov = (1./(n - 1.)) * moments[4];

                /* agrees with _cov[0][1]
                double[][] a = new double[n][];
                for (int i = 0; i < n; ++i) {
                    a[i] = new double[]{x[i], y[i]};
                }
                double[][] _cov = BruteForce.covariance(a);*/

                return cov;
            }
            default : break;
            // no processing needed for unit_stanc...
        }

        //x    y    x^2    y^2    xy   = moments
        //0    1    2      3       4   = indexes in output
        double[] moments = Util.caldc2DMomentsX2Y2(x, y);

        //    Cov(X, Y) = (1 / (n - 1)) * Σ ( (Xi - mean(X)) / std(X) ) * ( (Yi - mean(Y)) / std(Y) )

        // needs to be multiplied by stddev(x)*stdev(y) to transform to original data reference frame
        double cor = (1./(n - 1.)) * moments[4];

        return cor;
    }

    /**
     * calculate the sample mean of x as (1/(n-1)) * X^T*X
     * where the bias term (1/(n-1)) corrects for having the sample mean instead of the true population mean.
     @param x and n x m array of data where n is the number of samples and m is the number of variables (== dimensions).
     @param isZeroCentered true if x has already been zero-centered.
     @return the sample covariance
     */
    public static double[][] calcSampleCovariance(double[][] x, boolean isZeroCentered) {

        //NOTE: same result as BruteForce.covariance

        if (!isZeroCentered) {
            x = Standardization.zeroCenterMean(x);
        }

        double[][] s = MatrixUtil.createATransposedTimesA(x);
        MatrixUtil.multiply(s, 1./(x.length - 1));

        return s;
    }
}
