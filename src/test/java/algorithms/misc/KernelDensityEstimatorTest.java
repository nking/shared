package algorithms.misc;

import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

import java.io.IOException;
import java.util.Arrays;

public class KernelDensityEstimatorTest extends TestCase {

    /**
     * get a very small test data set. the values are already sorted, ascending.
     * @return
     */
    private double[] getData0() {
        // test data from https://en.wikipedia.org/wiki/Kernel_density_estimation#Example
        return new double[]{-2.1, -1.3, -0.4, 1.9, 5.1, 6.2};
    }

    private double[] getCenteredData0() {
        double[] data = Arrays.copyOf(getData0(), getData0().length);
        double[] meanStDev = MiscMath0.getAvgAndStDev(data);
        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(data);
        double m = meanStDev[0];
        m = medianAndIQR[0];
        for (int i = 0; i < data.length; ++i) {
            data[i] -= m;
        }
        return data;
    }

    public void testOptimalBandwidthGaussianKernel() {
        double[] data;
        data = getData0();

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

    public void testViaFFTGaussKernel() throws IOException {

        double[] data;
        //data = getData0();
        data = getCenteredData0();

        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(data);
        double kMAD = 1.4826;
        double s = kMAD*mADMinMax[0];
        double r0 = mADMinMax[1] - 3*s;
        double r1 = mADMinMax[1] + 3*s;

        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(data);

        double h = KernelDensityEstimator.ruleOfThumbBandwidthGaussianKernel(s, medianAndIQR[1],
                data.length);

        KernelDensityEstimator.KDE kde;

        kde = KernelDensityEstimator.viaFFTGaussKernel(data, h);
        //kde = KernelDensityEstimator.viaFFTGaussKernel(data, h, 16, 2,
        //        data[0] - 4*h, data[data.length - 1] + 4*h);

        //System.out.printf(" x =%s", FormatArray.toString(kde.hx, "%14.5e"));
        //System.out.printf("kde=%s", FormatArray.toString(kde.kde, "%14.5e"));

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));
        for (int i = 0; i < 10; ++i) {
            h /= 2.;
            kde = KernelDensityEstimator.viaFFTGaussKernel(kde.u, kde.hx, h);
            plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));
        }

        plotter.writeFile("kde_0");
    }

    public void testCreateFineHistogram() {
    }

    public void testTestViaFFTGaussKernel() {
    }

    public void testUnivariateKernelDensityEstimate() {
    }
}