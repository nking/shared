package algorithms.statistics;

import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import junit.framework.TestCase;

import java.util.logging.Logger;

public class GeneralizedExtremeValue2Test extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public void testGenerateGumberCurve() throws Exception {
        
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
    
}
