package algorithms.misc;

import algorithms.util.FormatArray;
import junit.framework.TestCase;

public class KernelDensityEstimatorTest extends TestCase {

    /**
     * get a very small test data set. the values are already sorted, ascending.
     * @return
     */
    private double[] getData0() {
        // test data from https://en.wikipedia.org/wiki/Kernel_density_estimation#Example
        return new double[]{-2.1, -1.3, -0.4, 1.9, 5.1, 6.2};
    }

    public void testOptimalBandwidthGaussianKernel() {
        double[] data = getData0();

        double[] meanStDev = MiscMath0.getAvgAndStDev(data);

        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(data);
        double kMAD = 1.4826;
        double s = kMAD*mADMinMax[0];
        double r0 = mADMinMax[1] - 3*s;
        double r1 = mADMinMax[1] + 3*s;

        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(data);

        double bw0 = KernelDensityEstimator.optimalBandwidthGaussianKernel(meanStDev[1], data.length);
        double bw1 = KernelDensityEstimator.optimalBandwidthGaussianKernel(s, data.length);

        double bwRT0 = KernelDensityEstimator.ruleOfThumbBandwidthGaussianKernel(meanStDev[1], medianAndIQR[1],
                data.length);
        double bwRT1 = KernelDensityEstimator.ruleOfThumbBandwidthGaussianKernel(s, medianAndIQR[1],
                data.length);

        assertTrue(bw0 > bwRT0);
        assertTrue(bw1 > bwRT1);
    }

    public void testViaFFTGaussKernel() {

        double[] data = getData0();

        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(data);
        double kMAD = 1.4826;
        double s = kMAD*mADMinMax[0];
        double r0 = mADMinMax[1] - 3*s;
        double r1 = mADMinMax[1] + 3*s;

        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(data);

        double h = KernelDensityEstimator.ruleOfThumbBandwidthGaussianKernel(s, medianAndIQR[1],
                data.length);

        KernelDensityEstimator.KDE kde = KernelDensityEstimator.viaFFTGaussKernel(data, h);

        System.out.printf("kde=%s", FormatArray.toString(kde.kde, "%14.5e"));
    }

    public void testCreateFineHistogram() {
    }

    public void testTestViaFFTGaussKernel() {
    }

    public void testUnivariateKernelDensityEstimate() {
    }
}