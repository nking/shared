package algorithms.statistics;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import junit.framework.TestCase;

import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;

public class GumbelTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public void testGenerateGumbelCurve() throws Exception {
        
        double[] xPoints;
        double[] curve;
        double sigma, k, mu;
        int n = 100;
        // TypeI
        double xMin = -100;
        double xMax = 100;
        xPoints = new double[n];
        double dx = (xMax - xMin)/n;
        xPoints[0] = xMin;
        for (int i = 1; i < n; i++) {
            xPoints[i] = xPoints[i - 1] + dx;
        }

        mu = 3.0;
        sigma = 7.0;
        k = 0.0;

        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        double[] randDist = Gumbel.sampleRandomlyFrom(mu, sigma, n, rand);
        Arrays.sort(randDist);
        double[] randDistParams = Gumbel.fitUsingMethodOfMoments(randDist);
        log.info(String.format("randDistParams=%s%n", FormatArray.toString(randDistParams, "%.3f")));

        curve = Gumbel.generateGumbelCurve(xPoints, mu, sigma);
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter((float)xMin, (float)xMax, 0f, 0.5f);
        plotter.addPlot(xPoints, curve, null, null, xPoints, curve, "gumbel: loc=" + mu +  " scale=" + sigma);

        // to generate the input for the estimators, we need to turn
        // the pair (xPoints, curve) into an ordered statistic for x,
        // that is, x points present in proportion to the curve (pdf).
        // TODO: add methods to do this properly, after take a quick rough look.
        //(Chap 19 of "Statistical Distributions").
                
        // where right tail begins to flatten
        double factor = 338.;
        //z=(x-mu)/sigma
        
        TDoubleList x = new TDoubleArrayList();
        int yC;
        for (int i = 0; i < curve.length; ++i) {
            yC = (int) Math.ceil(factor*curve[i]);
            for (int j = 0; j < yC; ++j) {
                x.add(xPoints[i]);
            }
        }
        double[] xOrd = x.toArray();

        double[] params = Gumbel.fitUsingMethodOfMoments(xOrd);
        log.info(String.format("params=%s%n", FormatArray.toString(params, "%.3f")));

        String filePath = plotter.writeFile("gumbel_mu1_sigma_1");

    }

    public void estSMC() throws Exception {

        double[] X = SMCFileReader.readDiffFile("smc118.1_diffs.txt");

        double sigma, k, mu;

        int n = X.length;

        mu = 3.0;
        sigma = 7.0;
        k = 0.0;

        int nBins = 100;

        int maxX = (int)Math.ceil(X[X.length - 1]);
        double binSize = maxX/nBins;
        HistogramHolder hist = Histogram.createSimpleHistogram(X, binSize, 0, maxX);
        float maxYHist = MiscMath0.findMax(hist.getYHistFloat());
        for (int i = 0; i < hist.getYHistFloat().length; ++i) {
            hist.getYHistFloat()[i] /= maxYHist;
        }
        double[] params = Gumbel.fitUsingMethodOfMoments(X);
        log.info(String.format("Gumbel MME params=%s%n", FormatArray.toString(params, "%.3f")));

        double tol = 1e-3;

        double[] params2 = Gumbel.fitUsingMaximumLikelihood(X);
        log.info(String.format("Gumbel MLE params=%s%n", FormatArray.toString(params2, "%.3f")));

        double[] params3 = GeneralizedExtremeValue.fitUsingMethodOfProbabilityWeightedMoments(X);
        log.info(String.format("GEV MPWME params=%s%n", FormatArray.toString(params3, "%.3f")));

        // matches results from scipy
        assertTrue(Math.abs(params2[0] - 1231.88073) < tol);
        assertTrue(Math.abs(params2[1] - 738.85117) < tol);

        // overplot a generated curve using params
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);

        double[] randX = Gumbel.sampleRandomlyFrom(params[0], params[1], 10*nBins, rand);
        HistogramHolder randHist = Histogram.createSimpleHistogram(randX, binSize, 0, maxX);
        maxYHist = MiscMath0.findMax(randHist.getYHistFloat());
        for (int i = 0; i < randHist.getYHistFloat().length; ++i) {
            randHist.getYHistFloat()[i] /= maxYHist;
        }

        double[] gev0x = MiscMath0.convertFloatToDouble(randHist.getXHist());
        double[] gev0y = GeneralizedExtremeValue.genCurve(gev0x, params3[0], params3[1], params3[2]);
        double yMax = MiscMath0.findMax(gev0y);
        MatrixUtil.multiply(gev0y, 1./yMax);

        // compare to scipy: -5.190956744208911, 3.207660475255879, 7.131822972996646
        double[] gev1y = GeneralizedExtremeValue._genCurve(gev0x, -5.190956744208911, 3.207660475255879, 7.131822972996646);
        yMax = MiscMath0.findMax(gev1y);
        MatrixUtil.multiply(gev1y, 1./yMax);

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter((float)0, (float)maxX, 0f, 1.f);

        plotter.addPlot((float)0, (float)maxX, 0f, 1.1f,
                hist.getXHist(), hist.getYHistFloat(),
                null, null,
                randHist.getXHist(), randHist.getYHistFloat(), "loc=" + params[0] + ", scale=" + params[1]);

        plotter.addPlot((float)0, (float)maxX, 0f, 1.1f,
                hist.getXHist(), hist.getYHistFloat(),
                null, null,
                randHist.getXHist(), MiscMath0.convertDoubleToFloat(gev0y),
                "gev loc=" + params3[0] + ", scale=" + params3[1] + ", shape=" + params3[2]);

        plotter.addPlot((float)0, (float)maxX, 0f, 1.1f,
                hist.getXHist(), hist.getYHistFloat(),
                null, null,
                randHist.getXHist(), MiscMath0.convertDoubleToFloat(gev1y),
                "gev loc=-5.19, scale=3.21, shape=7.13");

        String filePath = plotter.writeFile("gumbel_smc");

    }
}
