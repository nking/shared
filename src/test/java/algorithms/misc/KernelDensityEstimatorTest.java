package algorithms.misc;

import algorithms.correlation.UnivariateDistance;
import algorithms.imageProcessing.FFTUtil;
import algorithms.statistics.SMCFileReader;
import algorithms.statistics.UnivariateNormalDistribution;
import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import junit.framework.TestCase;
import thirdparty.ca.uol.aig.fftpack.Complex1D;

import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
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

    public void estOptimalBandwidthGaussianKernel() {
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

    private double[] calcH(double[] data) {

        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(data);
        double kMAD = 1.4826;
        double s = kMAD*mADMinMax[0];
        double r0 = mADMinMax[1] - 3*s;
        double r1 = mADMinMax[1] + 3*s;

        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(data);

        double rt = KernelDensityEstimator.ruleOfThumbBandwidthGaussianKernel(s, medianAndIQR[1],
                data.length);
        double ob = KernelDensityEstimator.optimalBandwidthGaussianKernel(s, data.length);

        System.out.printf("stdev=%.4e,  s=%.4e,  ob=%.4e,  rt=%.4e\n",
                MiscMath0.getAvgAndStDev(data)[1], s, ob, rt);
        return new double[]{rt, ob};
    }

    public double[] getRandomGaussianData0() throws NoSuchAlgorithmException {

        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //seed = 175696292403861L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);

        double mean1 = 0.3;
        double mean2 = 0.8;
        double sigma1 = 0.13;
        double sigma2 = 0.13;
        int n1 = 100;
        int n2 = 3*n1;
        double h = 0.5;

        double[] g1 = UnivariateNormalDistribution.randomSampleOf(mean1, sigma1, rand, n1);
        double[] g2 = UnivariateNormalDistribution.randomSampleOf(mean2, sigma2, rand, n2);
        double[] data = Arrays.copyOf(g1, g1.length + g2.length);
        System.arraycopy(g2, 0, data, g1.length, g2.length);

        return data;
    }

    public void est0() throws IOException {

        double[] data = getData0();

        double h = calcH(data)[0];

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

    public void testViaFFTGaussKernel0() throws IOException, NoSuchAlgorithmException {

        double[] data;
        //data = getData0();
        //data = getCenteredData0();
        //data = getSMCData();
        data = getRandomGaussianData0();

        double h = calcH(data)[0];
        h *= 10;

        System.out.printf("data.length=%d h=%.4e\n", data.length, h);

        double[][] hist = KernelDensityEstimator.createFineHistogram(data, h);

        double sum = KernelDensityEstimator.sumHistogram(hist);
        assertTrue(Math.abs(sum - 1.) < 1e-7);
        System.out.printf("dx=%.4e  sumHist=%.4e  data.length=%d  hist[0].length=%d\n",
                hist[0][1] - hist[0][0], sum, data.length, hist[0].length);

        KernelDensityEstimator.KDE kde;

        kde = KernelDensityEstimator.viaFFTGaussKernel(data, h);

        double cv = KernelDensityEstimator.crossValidationScore(kde.u, kde.hx, h);

        //System.out.printf(" x =%s", FormatArray.toString(kde.hx, "%14.5e"));
        //System.out.printf("kde=%s", FormatArray.toString(kde.kde, "%14.5e"));

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f, cv=%.4f", h, cv));
        PolygonAndPointPlotter plotter2 = new PolygonAndPointPlotter();
        //PolygonAndPointPlotter plotter2 = new PolygonAndPointPlotter(0, 1.f, 0, 1.f);
        plotter2.addPlot(hist[0], hist[1], hist[0], hist[1], String.format("histogram"));

        for (int i = 0; i < 40; ++i) {
            //h /= 2.;
            h *= 0.85;
            kde = KernelDensityEstimator.viaFFTGaussKernel(data, h);
            cv = KernelDensityEstimator.crossValidationScore(kde.u, kde.hx, h);
            plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f, cv=%.4f", h, cv));
        }

        plotter.writeFile("kde_0");
        plotter2.writeFile("kde_0_1");
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

        double h = calcH(data)[0];

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

        double cv = KernelDensityEstimator.crossValidationScore(hf.u, hf.hx, h);

        //System.out.printf(" x =%s", FormatArray.toString(kde.hx, "%14.5e"));
        //System.out.printf("kde=%s", FormatArray.toString(kde.kde, "%14.5e"));

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f, cv=%.4f", h, cv));
        for (int i = 0; i < 10; ++i) {
            h /= 2.;
            kde = KernelDensityEstimator.viaFFTGaussKernel(hf.u, hf.hx, h);
            cv = KernelDensityEstimator.crossValidationScore(hf.u, hf.hx, h);
            plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f, cv=%.4f", h, cv));
        }
        plotter.writeFile("kde_2");
    }

    public void testCreateFineHistogram() throws NoSuchAlgorithmException, IOException {
        double[] data = getRandomGaussianData0();
        double h = calcH(data)[1];
        h = 0.2;
        double[][] hist = KernelDensityEstimator.createFineHistogram(data, h);
        double sum = KernelDensityEstimator.sumHistogram(hist);
        assertTrue(Math.abs(sum - 1.) < 1e-7);
        System.out.printf("dx=%.4e  sumHist=%.4e  data.length=%d  hist[0].length=%d\n",
                hist[0][1] - hist[0][0], sum, data.length, hist[0].length);

        KernelDensityEstimator.KDE kde = KernelDensityEstimator.viaFFTGaussKernel(data, h);
        double sumK = KernelDensityEstimator.sumKDE(kde.kde, (kde.hx[1] - kde.hx[0]));
        System.out.printf("dx=%.4e  sumKDE=%.4e  data.length=%d  hist[1].length=%d\n",
                hist[0][1] - hist[0][0], sumK, data.length, hist[1].length);

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(hist[0], hist[1], hist[0], hist[1], String.format("histogram, h=%.4e", h));

        plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("kde, h=%.4e", h));
        plotter.writeFile("kde_hist_kde");

    }

    public void estUnivariateKernelDensityEstimate() {
    }
}