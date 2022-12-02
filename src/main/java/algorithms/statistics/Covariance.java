package algorithms.statistics;

import algorithms.correlation.MultivariateDistance;
import algorithms.matrix.MatrixUtil;

public class Covariance {

    /**
     * calculate the sample mean of x as (1/(n-1)) * X^T*X
     * where the bias term (1/(n-1)) corrects for having the sample mean instead of the true population mean.
     * @param x and n x m array of data where n is the number of samples and m is the number of variables (== dimensions).
     * @param isZeroCentered true if x has already been zero-centered.
     * @return the sample covariance
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
