package algorithms.statistics;

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
    
    public void estGenerateGumbelCurve() throws Exception {
        
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
        double[] randDistParams = Gumbel.fitGumbelUsingMethodOfMoments(randDist);
        log.info(String.format("randDistParams=%s\n", FormatArray.toString(randDistParams, "%.3f")));

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

        double[] params = Gumbel.fitGumbelUsingMethodOfMoments(xOrd);
        log.info(String.format("params=%s\n", FormatArray.toString(params, "%.3f")));

        String filePath = plotter.writeFile("gumbel_mu1_sigma_1");

    }

    public void testSMC() throws Exception {

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
        double[] params = Gumbel.fitGumbelUsingMethodOfMoments(X);
        log.info(String.format("MME params=%s\n", FormatArray.toString(params, "%.3f")));

        double[] params2 = Gumbel.fitGumbelUsingML(X);
        log.info(String.format("MLE params=%s\n", FormatArray.toString(params2, "%.3f")));


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

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter((float)0, (float)maxX, 0f, 1.f);

        plotter.addPlot((float)0, (float)maxX, 0f, 1.1f,
                hist.getXHist(), hist.getYHistFloat(),
                null, null,
                randHist.getXHist(), randHist.getYHistFloat(), "loc=" + params[0] + ", scale=" + params[1]);

        String filePath = plotter.writeFile("gumbel_smc");

    }
}
