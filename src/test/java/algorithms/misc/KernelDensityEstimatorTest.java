package algorithms.misc;

import algorithms.correlation.UnivariateDistance;
import algorithms.imageProcessing.FFTUtil;
import algorithms.statistics.SMCFileReader;
import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import junit.framework.TestCase;
import thirdparty.ca.uol.aig.fftpack.Complex1D;

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

    private double[] getSMCData() throws IOException {
        // 8_435_778 lines, in ascending order
        // 8_435_778 lines of data.  n*log_2(n) = 194_091_137 so not limited by most processors for algorithms in KDE.
        // For 1 GB available memory in the JVM, can hold 14 copies of the data as double (double is 8 Bytes).
        return SMCFileReader.readDiffFile("smc118.1_diffs.txt");
    }
    private double[] getCenteredSMCData() throws IOException {
        double[] data = getSMCData();
        data = Arrays.copyOf(data, data.length);
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

    private double calcH(double[] data) {

        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(data);
        double kMAD = 1.4826;
        double s = kMAD*mADMinMax[0];
        double r0 = mADMinMax[1] - 3*s;
        double r1 = mADMinMax[1] + 3*s;

        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(data);

        return KernelDensityEstimator.ruleOfThumbBandwidthGaussianKernel(s, medianAndIQR[1],
                data.length);
    }

    public void est0() throws IOException {

        double[] data = getData0();

        double h = calcH(data);

        // h *= 10;

        KernelDensityEstimator.KDE kde = null;

        //kde = KernelDensityEstimator.cVTerm2(data, h);

        //System.out.printf(" x =%s", FormatArray.toString(kde.hx, "%14.5e"));
        //System.out.printf("kde=%s", FormatArray.toString(kde.kde, "%14.5e"));

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));
        for (int i = 0; i < 10; ++i) {
            h /= 2.;
            kde = KernelDensityEstimator.viaFFTGaussKernel(kde.u, kde.hx, h);
            plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));
        }
        plotter.writeFile("kde_4");

    }

    public void testViaFFTGaussKernel0() throws IOException {

        double[] data;
        //data = getData0();
        data = getCenteredData0();

        double h = calcH(data);

       // h *= 10;

        KernelDensityEstimator.KDE kde;

        kde = KernelDensityEstimator.viaFFTGaussKernel(data, h);
        //kde = KernelDensityEstimator.viaFFTGaussKernel(data, h, 16, 2,
        //        data[0] - 4*h, data[data.length - 1] + 4*h);

        //System.out.printf(" x =%s", FormatArray.toString(kde.hx, "%14.5e"));
        //System.out.printf("kde=%s", FormatArray.toString(kde.kde, "%14.5e"));

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));
        PolygonAndPointPlotter plotter2 = new PolygonAndPointPlotter();
        plotter2.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));

        for (int i = 0; i < 10; ++i) {
            kde = KernelDensityEstimator.viaFFTGaussKernel(data, h);
            plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));

            kde = KernelDensityEstimator.viaFFTGaussKernel(kde.u, kde.hx, h);
            plotter2.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));
            h /= 2.;
        }

        plotter.writeFile("kde_0");
        plotter.writeFile("kde_0_1");
    }

    private static class HistAndFFT {
        double[] hx;
        double[] hy;
        Complex[] u;
        double h;
    }
    private HistAndFFT getSMCHist() throws IOException {
        double[] data;
        data = getSMCData();
        //data = getCenteredSMCData();

        double h = calcH(data);

        // the histogram length (hist[1].length) is a power of 2 and the range is enlarged by at least 2*3*h
        double[][] hist = KernelDensityEstimator.createFineHistogram(data, h);

        HistAndFFT hf = new HistAndFFT();
        hf.hx = Arrays.copyOf(hist[0], hist[0].length);
        hf.hy = Arrays.copyOf(hist[1], hist[1].length);
        hf.h = h;

        return hf;
    }

    private HistAndFFT getSMCFFTHist() throws IOException {

        HistAndFFT hf = getSMCHist();

        System.gc();

        // normalization is performed by default:
        FFTUtil fft = new FFTUtil();
        hf.u = fft.create1DFFT(hf.hy, true);

        return hf;
    }

    public void estViaFFTGaussKernel2() throws IOException {

        HistAndFFT hf = getSMCFFTHist();
        System.gc();

        KernelDensityEstimator.KDE kde;

        double h = hf.h;

        kde = KernelDensityEstimator.viaFFTGaussKernel(hf.u, hf.hx, h);

        //System.out.printf(" x =%s", FormatArray.toString(kde.hx, "%14.5e"));
        //System.out.printf("kde=%s", FormatArray.toString(kde.kde, "%14.5e"));

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));
        for (int i = 0; i < 10; ++i) {
            h /= 2.;
            kde = KernelDensityEstimator.viaFFTGaussKernel(hf.u, hf.hx, h);
            plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));
        }
        plotter.writeFile("kde_2");
    }

    public void testCreateFineHistogram() {
    }

    public void testUnivariateKernelDensityEstimate() {
    }
}