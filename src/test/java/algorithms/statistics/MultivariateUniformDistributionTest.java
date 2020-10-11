package algorithms.statistics;

import algorithms.correlation.UnivariateDistance;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath0;
import algorithms.misc.MiscSorter;
import algorithms.optimization.GeometricMedian;
import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;
import junit.framework.TestCase;
import thirdparty.dlib.optimization.GeometricMedianWeightedFunction;

/**
 *
 * @author nichole
 */
public class MultivariateUniformDistributionTest extends TestCase {
    
    public MultivariateUniformDistributionTest(String testName) {
        super(testName);
    }
    
    private SecureRandom rand = null;

    @Override
    protected void setUp() throws Exception {
        super.setUp(); 
        
        rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
    }
    
    
    public void test0() throws NoSuchAlgorithmException, IOException {
        
        System.out.println("test0");
                
        int nTests = 1;
        int nDimensions = 2;
        int nSamples = 250;
        boolean onSurface = true;
        
        double[][] v = new double[nSamples][nDimensions];
        double[][] vT;
        double[] vl = new double[nSamples*nDimensions];
            
        for (int ii = 0; ii < nTests; ++ii) {

            for (int i = 0; i < nSamples; ++i) {
                // a quick look at some of the stats:
                if (onSurface) {
                    v[i] = MultivariateUniformDistribution.generateOnUnitStandardNSphere(
                        nDimensions, rand);
                } else {
                    v[i] = MultivariateUniformDistribution.generateInUnitStandardNBall(
                        nDimensions, rand);
                }
                double[] avgAndStDev = MiscMath0.getAvgAndStDev(v[i]);
                System.out.println("\nv[" + i + "]=" + FormatArray.toString(v[i], "%.3f"));
                System.out.println("   mean, stDev=" + FormatArray.toString(avgAndStDev, "%.3f"));
                
                System.arraycopy(v[i], 0, vl, i*nDimensions, v[i].length);
            }
            vT = MatrixUtil.transpose(v);
            
            float[] col1 = new float[nSamples];
            float[] col2 = new float[nSamples];

            // rewrite into col1 and col2
            for (int j = 0; j < nSamples; ++j) {
                col1[j] = (float) v[j][0];
                col2[j] = (float) v[j][1];
            }            
            
            for (int i = 0; i < nDimensions; ++i) {
                double[] avgAndStDev = MiscMath0.getAvgAndStDev(vT[i]);
                System.out.printf("dimension[%d]: mean,stDev=%s\n", i, FormatArray.toString(avgAndStDev, "%.3f"));
            }
            
            double[] w = new double[nSamples];
            Arrays.fill(w, 1.0);
            double[] init = new double[nDimensions];

            GeometricMedian gm = new GeometricMedian();
            GeometricMedianWeightedFunction f = new GeometricMedianWeightedFunction(
                vl, nDimensions, w);
            
            double minSum = gm.newtonsThenVardiZhang(f, init);

            double[] diff = f.calculateDifferences(init);

            int[] oIdx = MiscSorter.mergeSortIncreasing(diff);

            double[] avgAndStDevFromGM = MiscMath0.getAvgAndStDev(diff);
            
            HistogramHolder h = Histogram.createSimpleHistogram(vl, 0.1, -3.5, 3.5);
            String str = h.plotHistogram("[X] onSphere=" + onSurface, "hist");
            h = Histogram.createSimpleHistogram(diff, 0.1, -3.5, 3.5);
            str = h.plotHistogram("[X-g.m.] onSphere=" + onSurface, "diff_hist");

            
            System.out.println("g.m.=" + FormatArray.toString(init, "%.3f"));
            System.out.printf("mean diff = %.3f, stdev diff = %.3f, minSum=%.3f\n",
                avgAndStDevFromGM[0], avgAndStDevFromGM[1], minSum);
          
            double min, max, binSz;
            UnivariateDistance.DCov dcov;
            UnivariateDistance.DCor dcor;
            
            dcor = UnivariateDistance.fastDcor(vT[0], vT[1]);
            System.out.println("correlation within points=" + Math.sqrt(dcor.corSq));
            
            // sum_i(sum_j(xi - xj)) and sum_i(sum_j(yi - yj))
            //    fall radially in counts from min value of x-axis
            dcov = UnivariateDistance.fastDcov(vT[0], vT[1]);
            System.out.println("cov:" + dcov.toString2());
            min = MiscMath0.findMin(dcov.ai);
            max = MiscMath0.findMax(dcov.ai);
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.ai, binSz, min, max);
            str = h.plotHistogram("[sum_i(sum_j(xi - xj))],onSphere=" + onSurface, "dist_hist_1");
            
            min = MiscMath0.findMin(dcov.bi);
            max = MiscMath0.findMax(dcov.bi);
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.bi, binSz, min, max);
            str = h.plotHistogram("sum_i(sum_j(yi - yj))],onSphere=" + onSurface, "dist_hist_2");
            
            min = MiscMath0.findMin(dcov.iv1);
            max = MiscMath0.findMax(dcov.iv1);
            double[] avgAndStDev2 = MiscMath0.getAvgAndStDev(dcov.iv1);
            double[] qs = MiscMath0.calcMedianAndIQR(dcov.iv1);
            //min = qs[0] - 1.5*qs[1];
            //max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv1, binSz, min, max);
            str = h.plotHistogram("[iv1=cumulative 1],onSphere=" + onSurface, "dist_hist_iv1");
            
            min = MiscMath0.findMin(dcov.iv2);
            max = MiscMath0.findMax(dcov.iv2);
            qs = MiscMath0.calcMedianAndIQR(dcov.iv2);
            //min = qs[0] - 1.5*qs[1];
            //max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv2, binSz, min, max);
            str = h.plotHistogram("[iv2=cumulative X],onSphere=" + onSurface, "dist_hist_iv2");
            
            min = MiscMath0.findMin(dcov.iv3);
            max = MiscMath0.findMax(dcov.iv3);
            qs = MiscMath0.calcMedianAndIQR(dcov.iv3);
            //min = qs[0] - 1.5*qs[1];
            //max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv3, binSz, min, max);
            str = h.plotHistogram("[iv3=cumulative Y],onSphere=" + onSurface, "dist_hist_iv3");
            
            min = MiscMath0.findMin(dcov.iv4);
            max = MiscMath0.findMax(dcov.iv4);
            qs = MiscMath0.calcMedianAndIQR(dcov.iv4);
            //min = qs[0] - 1.5*qs[1];
            //max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv4, binSz, min, max);
            str = h.plotHistogram("[iv4=cumulative X*Y],onSphere=" + onSurface, "dist_hist_iv4");
            
            
            PolygonAndPointPlotter plotter;
            float yMin = -1.2f;
            float yMax = 1.2f;
            float xMin = -1.2f;
            float xMax = 1.2f;

            plotter = new PolygonAndPointPlotter();
            plotter.addPlot(
                xMin, xMax, yMin, yMax, col1, col2, null, null, "d=2, onSurface=" + onSurface
            );
            plotter.writeFile("onSurface_" + onSurface);
        }
    }
    
    public void testRejectionMethod() throws NoSuchAlgorithmException, IOException {
        
        System.out.println("testRejectionMethod");
        
        double tol = 0.001;
        
        int d = 2;
        int n = 250;
        boolean[] onSurface;
        onSurface = new boolean[]{true, false};
        onSurface = new boolean[]{true};
        double[][] x;
        
        PolygonAndPointPlotter plotter;
        float yMin = -1.2f;
        float yMax = 1.2f;
        float xMin = -1.2f;
        float xMax = 1.2f;
        
        float[] col1 = new float[n];
        float[] col2 = new float[n];
        double[] col1d = new double[n];
        double[] col2d = new double[n];
        double[] allX = new double[n*d];
        double[][] xT;
        
        int i, j, k;
        double dist;
        for (i = 0; i < onSurface.length; ++i) {
            x = MultivariateUniformDistribution.generateUnitStandardNSphereWithRejection(
                d, n, rand, onSurface[i]);
            
            assertEquals(n, x.length);
            assertEquals(d, x[0].length);
            for (j = 0; j < n; ++j) {
                for (k = 0; k < d; ++k) {
                    assertTrue(x[j][k] <= 1.0);
                }
                dist = MultivariateUniformDistribution.distSquaredFromOrigin(x[j]);
                dist = Math.sqrt(dist);
                if (onSurface[i]) {
                    assertTrue(Math.abs(dist - 1.) < tol);
                } else {
                    assertTrue(dist <= (1. + tol));
                }
            }
            
            xT = MatrixUtil.transpose(x);
            
            // rewrite into col1 and col2
            for (j = 0; j < n; ++j) {
                col1[j] = (float) x[j][0];
                col2[j] = (float) x[j][1];
                col1d[j] = x[j][0];
                col2d[j] = x[j][1];
                allX[j*2] = x[j][0];
                allX[j*2 + 1] = x[j][1];
            }
            
            plotter = new PolygonAndPointPlotter();
            plotter.addPlot(
                xMin, xMax, yMin, yMax, col1, col2, null, null, "d=2, onSurface=" + onSurface[i]
            );
            plotter.writeFile("onSurface_" + onSurface[i] + "_rej2");
            
            //====
            UnivariateDistance.DCor dcor = UnivariateDistance.fastDcor(xT[0], xT[1]);
            System.out.println("correlation within points=" + Math.sqrt(dcor.corSq));
            
            UnivariateDistance.DCov dcov;
            double min, max, binSz;
            HistogramHolder h;
            String str;
            
            dcov = UnivariateDistance.fastDcov(col1d, col2d);
            System.out.println("rej2: cov:" + dcov.toString2());
            min = MiscMath0.findMin(dcov.ai);
            max = MiscMath0.findMax(dcov.ai);
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.ai, binSz, min, max);
            str = h.plotHistogram("rej2: [sum_i(sum_j(xi - xj))],onSphere=" + onSurface[i], "dist_hist_1_rej2");
            
            min = MiscMath0.findMin(dcov.bi);
            max = MiscMath0.findMax(dcov.bi);
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.bi, binSz, min, max);
            str = h.plotHistogram("rej2: sum_i(sum_j(yi - yj))],onSphere=" + onSurface[i], "dist_hist_2_rej2");
            
            min = MiscMath0.findMin(dcov.iv1);
            max = MiscMath0.findMax(dcov.iv1);
            double[] avgAndStDev2 = MiscMath0.getAvgAndStDev(dcov.iv1);
            double[] qs = MiscMath0.calcMedianAndIQR(dcov.iv1);
            //min = qs[0] - 1.5*qs[1];
            //max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv1, binSz, min, max);
            str = h.plotHistogram("rej2: [iv1=cumulative 1],onSphere=" + onSurface[i], "dist_hist_iv1_rej2");
            
            min = MiscMath0.findMin(dcov.iv2);
            max = MiscMath0.findMax(dcov.iv2);
            qs = MiscMath0.calcMedianAndIQR(dcov.iv2);
            //min = qs[0] - 1.5*qs[1];
            //max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv2, binSz, min, max);
            str = h.plotHistogram("rej2: [iv2=cumulative X],onSphere=" + onSurface[i], "dist_hist_iv2_rej2");
            
            min = MiscMath0.findMin(dcov.iv3);
            max = MiscMath0.findMax(dcov.iv3);
            qs = MiscMath0.calcMedianAndIQR(dcov.iv3);
            //min = qs[0] - 1.5*qs[1];
            //max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv3, binSz, min, max);
            str = h.plotHistogram("rej2: [iv3=cumulative Y],onSphere=" + onSurface[i], "dist_hist_iv3_rej2");
            
            min = MiscMath0.findMin(dcov.iv4);
            max = MiscMath0.findMax(dcov.iv4);
            qs = MiscMath0.calcMedianAndIQR(dcov.iv4);
            //min = qs[0] - 1.5*qs[1];
            //max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv4, binSz, min, max);
            str = h.plotHistogram("rej2: [iv4=cumulative X*Y],onSphere=" + onSurface[i], "dist_hist_iv4_rej2");
            
            // ====
            double[] w = new double[n];
            Arrays.fill(w, 1.0);
            double[] init = new double[d];

            GeometricMedian gm = new GeometricMedian();
            GeometricMedianWeightedFunction f = new GeometricMedianWeightedFunction(
                allX, d, w);
            
            double minSum = gm.newtonsThenVardiZhang(f, init);

            double[] diff = f.calculateDifferences(init);

            int[] oIdx = MiscSorter.mergeSortIncreasing(diff);

            double[] avgAndStDevFromGM = MiscMath0.getAvgAndStDev(diff);
            
            h = Histogram.createSimpleHistogram(allX, 0.1, -3.5, 3.5);
            str = h.plotHistogram("rej2: [X] onSphere=" + onSurface[i], "hist_rej2");
            h = Histogram.createSimpleHistogram(diff, 0.1, -3.5, 3.5);
            str = h.plotHistogram("rej2: [X-g.m.] onSphere=" + onSurface[i], "diff_hist_rej2");

            
            System.out.println("rej2: g.m.=" + FormatArray.toString(init, "%.3f"));
            System.out.printf("rej2: mean diff = %.3f, stdev diff = %.3f, minSum=%.3f\n",
                avgAndStDevFromGM[0], avgAndStDevFromGM[1], minSum);
        }
         
    }
    
    /**
     * using the Rayleigh Test from the paper
     * "An Overview of Uniformity Tests on the Hypersphere"
     *  Garcia-Portugues and Verdebout, 2018, arXiv:1804.00286.
     * 
     * @throws NoSuchAlgorithmException
     * @throws IOException 
     */
    public void testRejectionMethodUnitCircleRayleighTest() throws NoSuchAlgorithmException, IOException {
          
        //The asymptotic distribution of Rn under H0 is a χ2, a chi-squared
        //distribution with 2 degrees of freedom.
        //
        //Rn = (2/n) * [
        //         (summation_{i=1_to_n}(cos(theta_i))^2 
        //         + (summation_{i=1_to_n}(sin(theta_i))^2 
        //     ];
        // --> looks like cos(theta_i) is x_i in the unit circle and
        //        sin(theta_i) is y_i in the unit circle
        //
        
        int d = 2;
        int n = 250;
        boolean[] onSurface = new boolean[] {false, true};
        
        for (boolean onS : onSurface) {
            
            double[][] x = MultivariateUniformDistribution
                .generateUnitStandardNSphereWithRejection(d, n, rand, onS);

            int i;
            double sum1 = 0;
            double sum2 = 0;
            for (i = 0; i < n; ++i) {
                sum1 += x[i][0];
                sum2 += x[i][1];
            }
            double rn = (2./n) * (sum1*sum1 + sum2*sum2);

            /*
            https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm

            chisq stat for deg of freedom=2

            For an upper-tail one-sided test, find the column corresponding to 
                1-α in the table containing upper-tail critical 
                 and reject the null hypothesis if the test statistic is 
                 greater than the tabled value.
            For a lower-tail one-sided test, find the column corresponding to α 
                in the lower-tail critical values table and reject the null 
                hypothesis if the computed test statistic is less than the 
                tabled value.

            Upper-tail critical values of chi-square distribution with ν degrees of freedom
                            Probability less than the critical value
               ν           0.90      0.95     0.975      0.99     0.999

               1          2.706     3.841     5.024     6.635    10.828
           ==> 2          4.605     5.991     7.378     9.210    13.816

            Lower-tail critical values of chi-square distribution with ν degrees of freedom
                    Probability less than the critical value
               ν           0.10     0.05     0.025      0.01     0.001

               1.          .016      .004      .001      .000      .000
           ==> 2.          .211      .103      .051      .020      .002

            */
            assertTrue(rn < 5.991);
        }
        
    }
}
