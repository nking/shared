package algorithms.misc;

import algorithms.imageProcessing.FFTUtil;
import algorithms.matrix.MatrixUtil;
import algorithms.statistics.SMCFileReader;
import algorithms.statistics.Standardization;
import algorithms.statistics.UnivariateNormalDistribution;
import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.Random;

public class KernelDensityEstimatorTest extends TestCase {

    private Random rand = null;

    public KernelDensityEstimatorTest() {
        long seed = System.nanoTime();
        //seed = 188903032980675L;
        //byte[] bytes = ByteBuffer.allocate(Long.SIZE / Byte.SIZE).putLong(seed).array();
        System.out.println("SEED=" + seed);
        rand = new Random(seed);
    }

    /**
     * get a very small test data set. the values are already sorted, ascending.
     *
     @return
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
        double s = kMAD * mADMinMax[0];
        double r0 = mADMinMax[1] - 3 * s;
        double r1 = mADMinMax[1] + 3 * s;

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
        double s = kMAD * mADMinMax[0];
        double r0 = mADMinMax[1] - 3 * s;
        double r1 = mADMinMax[1] + 3 * s;

        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(data);

        double rt = KernelDensityEstimator.ruleOfThumbBandwidthGaussianKernel(s, medianAndIQR[1],
                data.length);
        double ob = KernelDensityEstimator.optimalBandwidthGaussianKernel(s, data.length);

        System.out.printf("stdev=%.4e,  s=%.4e,  ob=%.4e,  rt=%.4e\n",
                MiscMath0.getAvgAndStDev(data)[1], s, ob, rt);
        return new double[]{rt, ob};
    }

    public double[] getRandomGaussianData0() throws NoSuchAlgorithmException {

        double mean1 = 0.3;
        double mean2 = 0.8;
        double sigma1 = 0.13;
        double sigma2 = 0.13;
        int n1 = 100;
        int n2 = 3 * n1;

        return getRandomGaussianData(new double[]{mean1, mean2}, new double[]{sigma1, sigma2},
                new int[]{n1, n2});
    }

    public double[] getRandomGaussianData(double[] mean, double[] sigma, int[] n) throws NoSuchAlgorithmException {

        int nTotal = 0;
        int i;
        for (i = 0; i < n.length; ++i) {
            nTotal += n[i];
        }

        double[] data = new double[nTotal];
        double[] g;
        int d = 0;
        for (i = 0; i < n.length; ++i) {
            g = UnivariateNormalDistribution.randomSampleOf(mean[i], sigma[i], rand, n[i]);
            System.arraycopy(g, 0, data, d, g.length);
            d += g.length;
        }
        return data;
    }

    public void test0() throws IOException {

        double[] data = getData0();

        double h = calcH(data)[0];

        KernelDensityEstimator.KDE kde = null;

        kde = KernelDensityEstimator.viaFFTGaussKernel(data, h);

        //System.out.printf(" x =%s", FormatArray.toString(kde.hx, "%14.5e"));
        //System.out.printf("kde=%s", FormatArray.toString(kde.kde, "%14.5e"));

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));
        for (int i = 0; i < 10; ++i) {
            h /= 1.5;
            kde = KernelDensityEstimator.viaFFTGaussKernel(kde.u, kde.hx, h);
            plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f", h));
        }
        plotter.writeFile("kde_4");

    }

    private double[][] standardize(double[] data) {
        double[] mean = new double[1];
        double[] stDev = new double[1];
        double[] data2 = Standardization.standardUnitNormalization(data, 1, mean, stDev);
        double[][] s = new double[3][];
        s[0] = data2;
        s[1] = mean;
        s[2] = stDev;
        System.out.println("standardized mean, stdev=" + FormatArray.toString(
                MiscMath0.getAvgAndStDev(data2), "%.4e"));
        return s;
    }

    public void testViaFFTGaussKernel0() throws IOException, NoSuchAlgorithmException {

        double[] data;
        //data = getData0();
        //data = getCenteredData0();
        //data = getSMCData();
        //data = getRandomGaussianData0();
        data = getRandomGaussianData(new double[]{100, 200, 400}, new double[]{25, 10, 25},
                new int[]{100, 30, 100});

        // data[0] is unit standardized data, data[1] is the mean of data, data[2] is the standard deviation of the mean of data
        double[][] unitStandardized = standardize(data);
        data = unitStandardized[0];

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

        double minCV = cv;
        int minCVIdx = -1;
        for (int i = 0; i < 40; ++i) {
            //h /= 2.;
            h *= 0.85;
            kde = KernelDensityEstimator.viaFFTGaussKernel(data, h);
            System.out.println("i=" + i + ":");
            cv = KernelDensityEstimator.crossValidationScore(kde.u, kde.hx, h);
            plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f, cv=%.4f", h, cv));
            if (cv < minCV) {
                minCV = cv;
                minCVIdx = i;
            }
        }

        System.out.printf("min cv=%.4e , i=%d\n", minCV, minCVIdx);

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

    public void estViaFFTGaussKernelSMC() throws IOException {

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

    private double[] getDataEvery3rdBin() {

        // with a dataset that is populated in every 3rd bin,
        //    a bandwidth selection method should be able to tell
        //    that a kde formed with a bin size h should be better
        //    than a kde formed with a bin size >= 2*h

        int n = 100;
        double binSize0 = 0.2;
        double bw0 = 0.2;
        //[0]=0:0.2, [3]=0.6:0.8, [6]=1.2:1.4,[9]=1.8:2.0, ...[(10-1)*3=27]=5.4:5.6
        // range = 0:5.6
        int nBins0 = 28;// 5.6/0.2

        double[] expectedCount = new double[nBins0];
        double[] data = new double[n];
        int r;
        int bin;
        int i;

        double r2;
        for (i = 0; i < n; ++i) {
            //populate 10 bins which are spaced at intervals of every 3rd 0.2.
            r = rand.nextInt(10);
            bin = r * 3;
            expectedCount[bin]++;

            // generate a number between [bin*0.2, bin*0.2 + 0.2)

            /* this is code from open jdk for ThreadLocalRandom.current().nextDouble(bound):
            which states:
            Written by Doug Lea with assistance from members of JCP JSR-166
            Expert Group and released to the public domain, as explained at
            http://creativecommons.org/publicdomain/zero/1.0/
            */
            // rand.nextDouble(bw0):
            r2 = (rand.nextLong() >>> 11) * 0x1.0p-53 * bw0;
            data[i] = (r2 < bw0) ? r2 : // correct for rounding
                    Double.longBitsToDouble(Double.doubleToLongBits(bw0) - 1);
            data[i] += bw0 * bin;
        }

        Arrays.sort(data);

        return data;
    }

    public void testViaFFTGaussKernel_compare() throws IOException, NoSuchAlgorithmException {
        System.out.printf("testViaFFTGaussKernel_compare");
        double[] data = getDataEvery3rdBin();
        int nPeaksExpected = 10;
        double h0 = calcH(data)[0];

        System.out.printf("data.length=%d h0=%.4e\n", data.length, h0);

        KernelDensityEstimator.KDE kdeFFT0 = KernelDensityEstimator.viaFFTGaussKernel(data, h0);

        double[] xGrid = kdeFFT0.hx;
        //xGrid = data;
        KernelDensityEstimator.KDE kde0 = KernelDensityEstimator.viaGaussKernel(data, h0);
        KernelDensityEstimator.KDE kde00 = KernelDensityEstimator.viaGaussKernel(data, xGrid, h0);

        float minX = (float)kdeFFT0.hx[0];
        float maxX = (float)kdeFFT0.hx[kdeFFT0.hx.length - 1];
        float minY = 0.f;
        float maxY = 1.1f * (float)MiscMath0.getMinMax(kde00.kde)[1];

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(kdeFFT0.hx, kdeFFT0.kde, null, null, String.format("FFT: h=%.4f", h0));

        float[] x0 = MiscMath0.convertDoubleToFloat(kde0.hx);
        float[] y0 = MiscMath0.convertDoubleToFloat(kde0.kde);
        plotter.addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(x0), Misc0.convertToNumberArray(y0),
                null, null, null, null,
                String.format("0: h=%.4f", h0));
        plotter.addPlot(kde00.hx, kde00.kde, null, null, String.format("00: h=%.4f", h0));

        kde0 = KernelDensityEstimator.viaGaussKernel(data, h0/120.);
        plotter.addPlot(kde0.hx, kde0.kde, null, null, String.format("0: h=%.4f", h0/120.));

        plotter.writeFile("_kde_comp_");
    }

    public void testViaFFTGaussKernel3() throws IOException, NoSuchAlgorithmException {

        System.out.printf("testViaFFTGaussKernel3");

        double[] data = getDataEvery3rdBin();

        int nPeaksExpected = 10;

        // data[0] is unit standardized data, data[1] is the mean of data, data[2] is the standard deviation of the mean of data
        //double[][] unitStandardized = standardize(data);
        //data = unitStandardized[0];

        // the smoothest distribution that still has nPeaks=10 is
        //  h = approximately 1.3 so we'll use that as the expected answer.
        //  h0 here is about 0.33 which is 3.9 times smaller than what we might expect.
        //  note also that the bandwidth selection score (which is currently fudged) is a minimum for h=0.5.
        double h0 = calcH(data)[0];

        // start at wider bandwidth:
        double h = 10*h0;

        System.out.printf("data.length=%d h0=%.4e h=%.4e\n", data.length, h0, h);

        KernelDensityEstimator.KDE kde;

        kde = KernelDensityEstimator.viaFFTGaussKernel(data, h);

        System.out.printf("dx=%.4e  data.length=%d  hist[0].length=%d\n",
                kde.hx[1] - kde.hx[0], data.length, kde.hx.length);

        double cv = KernelDensityEstimator.crossValidationScore(kde.u, kde.hx, h);

        //System.out.printf(" x =%s", FormatArray.toString(kde.hx, "%14.5e"));
        //System.out.printf("kde=%s", FormatArray.toString(kde.kde, "%14.5e"));

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4f, cv=%.4f", h, cv));

        MinMaxPeakFinder peakFinder = new MinMaxPeakFinder();

        double minCV = cv;
        int minCVIdx = 0;
        int[] peaks = peakFinder.findPeaks(MatrixUtil.convertToFloat(kde.kde));
        int nPeaks = peaks.length;
        int minPeaks = nPeaks;
        for (int i = 1; i < 40; ++i) {
            h *= 0.9;
            kde = KernelDensityEstimator.viaFFTGaussKernel(kde.u, kde.hx, h);
            nPeaks = peakFinder.findPeaks(MatrixUtil.convertToFloat(kde.kde)).length;
            System.out.println("i=" + i + ":");
            cv = KernelDensityEstimator.crossValidationScore(kde.u, kde.hx, h);
            plotter.addPlot(kde.hx, kde.kde, kde.hx, kde.kde, String.format("h=%.4e, cv=%.4e", h, cv));
            System.out.printf("t3: nPeaks=%d, h=%.4e, cv=%.4e\n", nPeaks, h, cv);
            if (cv < minCV) {
                minCV = cv;
                minCVIdx = i;
                minPeaks = nPeaks;
                System.out.printf(" * min h=%.4e, cv=%.4e\n", h, cv);
            }
        }

        System.out.printf("t3: min cv=%.4e , i=%d\n", minCV, minCVIdx);
        plotter.writeFile("kde_3");
    }
}