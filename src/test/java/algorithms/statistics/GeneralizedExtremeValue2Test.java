package algorithms.statistics;

import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import junit.framework.TestCase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;

public class GeneralizedExtremeValue2Test extends TestCase {

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
        double[] randDist = GumbelCDF.sampleRandomlyFrom(mu, sigma, n, rand);
        Arrays.sort(randDist);
        double[] randDistParams = GeneralizedExtremeValue.gumbelParamsViaMethodOfMoments(randDist);
        log.info(String.format("randDistParams=%s\n", FormatArray.toString(randDistParams, "%.3f")));

        curve = GeneralizedExtremeValue.generateGumbelCurve(xPoints, mu, sigma);
        
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

        double[] params = GeneralizedExtremeValue.gumbelParamsViaMethodOfMoments(xOrd);
        log.info(String.format("params=%s\n", FormatArray.toString(params, "%.3f")));

        String filePath = plotter.writeFile("gumbel_mu1_sigma_1");

    }

    public void estSMC() throws Exception {

        // read in the dataset
        String path = ResourceFinder.findFileInTestResources("smc118.1_diffs.txt");
        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        TDoubleList d = new TDoubleArrayList();
        try {
            in = new BufferedReader(new FileReader(new File(path)));
            String line = in.readLine();
            while (line != null) {
                d.add(Double.parseDouble(line));
                line = in.readLine();
            }
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }
        double[] X = d.toArray();
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
        double[] params = GeneralizedExtremeValue.gumbelParamsViaMethodOfMoments(X);
        log.info(String.format("params=%s\n", FormatArray.toString(params, "%.3f")));

        // overplot a generated curve using params
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);

        double[] randX = GumbelCDF.sampleRandomlyFrom(params[0], params[1], 10*nBins, rand);
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
