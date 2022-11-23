package algorithms.statistics;

import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

public class GeneralizedExtremeValueTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public void estGenerateCurve() throws Exception {
        
        double[] xPoints;
        double[] yPoints;
        double[] dXPoints;
        double[] dYPoints;
        GeneralizedExtremeValue gev;
        double[] curve;
        double sigma, k, mu;
        yPoints = new double[0];
        dXPoints = new double[0];
        dYPoints = new double[0];
        
        // TypeII
        xPoints = new double[11];
        for (int i = 0; i < 11; i++) {
            xPoints[i] = -4+i;
        }
        k = 0.5;
        sigma = 1.0;
        mu = 0.0;
        gev = new GeneralizedExtremeValue(xPoints, yPoints, dXPoints, dYPoints);
        curve = gev.generateCurve(xPoints, mu, sigma, k);
        assertTrue(Math.abs(curve[3] - 0.15f) < 0.01f);
        assertTrue(Math.abs(curve[4] - 0.37f) < 0.05f);
        assertTrue(Math.abs(curve[5] - 0.19f) < 0.05f);
        assertTrue(Math.abs(curve[6] - 0.097f) < 0.01f);
        assertTrue(Math.abs(curve[7] - 0.05f) < 0.01f);
        assertTrue(Math.abs(curve[8] - 0.03f) < 0.01f);
        
        curve = GeneralizedExtremeValue.genCurve(xPoints, mu, sigma, k);
        assertTrue(Math.abs(curve[3] - 0.15f) < 0.01f);
        assertTrue(Math.abs(curve[4] - 0.37f) < 0.05f);
        assertTrue(Math.abs(curve[5] - 0.19f) < 0.05f);
        assertTrue(Math.abs(curve[6] - 0.097f) < 0.01f);
        assertTrue(Math.abs(curve[7] - 0.05f) < 0.01f);
        assertTrue(Math.abs(curve[8] - 0.03f) < 0.01f);
        
        for (int i = 0; i < xPoints.length; i++) {
            Double a = GeneralizedExtremeValue.generateYGEV(xPoints[i], mu, sigma, k);
            switch(i) {
                case 3:
                    assertTrue(Math.abs(a - 0.15f) < 0.01f);
                    break;
                case 4:
                    assertTrue(Math.abs(a - 0.37f) < 0.05f);
                    break;
                case 5:
                    assertTrue(Math.abs(a - 0.19f) < 0.05f);
                    break;
                case 6:
                    assertTrue(Math.abs(a - 0.097f) < 0.01f);
                    break;
                case 7:
                    assertTrue(Math.abs(a - 0.05f) < 0.01f);
                    break;
                case 8:
                    assertTrue(Math.abs(a - 0.03f) < 0.01f);
                    break;
                default:
                    break;
            }
        }
        
        // TypeI
        k = 0.0;
        sigma = 1.0;
        mu = 0.0;
        gev = new GeneralizedExtremeValue(xPoints, yPoints, dXPoints, dYPoints);
        curve = gev.generateCurve(xPoints, mu, sigma, k);
        assertNotNull(curve);        
        assertTrue(Math.abs(curve[3] - 0.18f) < 0.1f);
        assertTrue(Math.abs(curve[4] - 0.37f) < 0.01f);
        assertTrue(Math.abs(curve[5] - 0.25f) < 0.01f);
        assertTrue(Math.abs(curve[6] - 0.12f) < 0.05f);
        assertTrue(Math.abs(curve[7] - 0.05f) < 0.01f);
        
        curve = GeneralizedExtremeValue.genCurve(xPoints, mu, sigma, k);
        assertTrue(Math.abs(curve[3] - 0.18f) < 0.1f);
        assertTrue(Math.abs(curve[4] - 0.37f) < 0.01f);
        assertTrue(Math.abs(curve[5] - 0.25f) < 0.01f);
        assertTrue(Math.abs(curve[6] - 0.12f) < 0.05f);
        assertTrue(Math.abs(curve[7] - 0.05f) < 0.01f);
        
        for (int i = 0; i < xPoints.length; i++) {
            Double a = GeneralizedExtremeValue.generateYEVTypeI(xPoints[i], mu, sigma);
            switch(i) {
                case 3:
                    assertTrue(Math.abs(a - 0.18f) < 0.1f);
                    break;
                case 4:
                    assertTrue(Math.abs(a - 0.37f) < 0.01f);
                    break;
                case 5:
                    assertTrue(Math.abs(a - 0.25f) < 0.01f);
                    break;
                case 6:
                    assertTrue(Math.abs(a - 0.12f) < 0.05f);
                    break;
                case 7:
                    assertTrue(Math.abs(a - 0.05f) < 0.01f);
                    break;
                default:
                    break;
            }
        }
    }

    public void testFitParameters() throws NoSuchAlgorithmException, IOException {

        double[] yGEV;
        double[] X;
        double[] expectedParams = new double[]{};
        double[] params;
        HistogramHolder hist;

        int i;
        int j;
        int n = 1000;
        int nBins = 100;

        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);

        // Gumbel: shape = 0
        // Frechet: shape = 1/location and both are > 0.  and y=1+shape*(x-location/scale) > 0
        // Weibull: shape = 1/location and both are < 0.  and y=1+shape*(x-location/scale) > 0

        double[][] modelParams = new double[3][];
        modelParams[0] = new double[]{2., 0.25, 0}; // Gumbel
        modelParams[1] = new double[]{2., 0.25, 0.5}; // Frechet
        modelParams[2] = new double[]{-2., 0.25, -0.5}; // reverse Weibull

        float maxYHist;
        double binSize;

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        for (i = 0; i < modelParams.length; ++i) {
            expectedParams = modelParams[i];
            X = GeneralizedExtremeValue.sampleRandomlyFrom(expectedParams[0], expectedParams[1], expectedParams[2], n);
            Arrays.sort(X);

            params = GeneralizedExtremeValue.fitUsingMethodOfProbabilityWeightedMoments(X);

            binSize = (X[n-1] - X[0])/nBins;

            // plot the distribution and parameters
            hist = Histogram.createSimpleHistogram(X, binSize, X[0], X[X.length - 1]);
            maxYHist = MiscMath0.findMax(hist.getYHistFloat());
            for (j = 0; j < hist.getYHistFloat().length; ++j) {
                hist.getYHistFloat()[j] /= maxYHist;
            }

            yGEV = GeneralizedExtremeValue.genCurve(MiscMath0.convertFloatToDouble(hist.getXHist()),
                    expectedParams[0], expectedParams[1], expectedParams[2]);
            MatrixUtil.multiply(yGEV, 1./MiscMath0.findMax(yGEV));

            plotter.addPlot(hist.getXHist()[0], hist.getXHist()[hist.getXHist().length - 1],
                    0f, 1.1f,
                    hist.getXHist(), hist.getYHistFloat(),
                    null, null,
                    hist.getXHist(), MiscMath0.convertDoubleToFloat(yGEV),
                    "loc=" + expectedParams[0] + ", scale=" + expectedParams[1] +
                    " shape=" + expectedParams[2]);

            System.out.printf("given params=%s%n", FormatArray.toString(expectedParams, "%.3e"));
            System.out.printf("estimated params=%s%n", FormatArray.toString(params, "%.3e"));
        }
        String filePath = plotter.writeFile("gev_3_tests");
    }
}
