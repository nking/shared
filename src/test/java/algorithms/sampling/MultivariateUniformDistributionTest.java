package algorithms.sampling;

import algorithms.correlation.UnivariateDistance;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath0;
import algorithms.misc.MiscSorter;
import algorithms.optimization.GeometricMedian;
import algorithms.util.FormatArray;
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
    
    public void test0() throws NoSuchAlgorithmException, IOException {
        
        System.out.println("test0");
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        int nTests = 1;
        int nDimensions = 5;
        int nSamples = 100;
        boolean onSurface = true;
        
        double[][] v = new double[nSamples][nDimensions];
        double[] vl = new double[nSamples*nDimensions];
            
        for (int ii = 0; ii < nTests; ++ii) {

            for (int i = 0; i < nSamples; ++i) {
                // a quick look at some of the stats:
                if (onSurface) {
                    v[i] = MultivariateUniformDistribution
                        .generateUnitStandardOnNSphere(nDimensions, rand);
                } else {
                    v[i] = MultivariateUniformDistribution.generateUnitStandardInNBall(
                    nDimensions, rand);
                }
                double[] avgAndStDev = MiscMath0.getAvgAndStDev(v[i]);
                System.out.println("\nv[" + i + "]=" + FormatArray.toString(v[i], "%.3f"));
                System.out.println("   mean, stDev=" + FormatArray.toString(avgAndStDev, "%.3f"));
                System.arraycopy(v[i], 0, vl, i*nDimensions, 
                    v[i].length);
            }
            double[] avgAndStDev = MiscMath0.getAvgAndStDev(vl);
            System.out.println("all: mean, stDev=" + FormatArray.toString(avgAndStDev, "%.3f"));

            double[] w = new double[nDimensions];
            Arrays.fill(w, 1.0);
            double[] init = new double[nDimensions];

            GeometricMedian gm = new GeometricMedian();
            GeometricMedianWeightedFunction f = new GeometricMedianWeightedFunction(
                vl, nDimensions, w);
            
            double minSum = gm.newtonsThenVardiZhang(f, init);

            double[] diff = f.calculateDifferences(init);

            int[] oIdx = MiscSorter.mergeSortIncreasing(diff);

            double[] avgAndStDevFromGM = MiscMath0.getAvgAndStDev(diff);
            
            HistogramHolder h = Histogram.createSimpleHistogram(vl, 0.1, -1.1, 1.1);
            String str = h.plotHistogram("[X] onSphere=" + onSurface, "hist");
            h = Histogram.createSimpleHistogram(diff, 0.1, -1.1, 1.1);
            str = h.plotHistogram("[X-g.m.] onSphere=" + onSurface, "diff_hist");

            
            System.out.println("g.m.=" + FormatArray.toString(init, "%.3f"));
            System.out.printf("mean diff = %.3f, stdev diff = %.3f, minSum=%.3f\n",
                avgAndStDevFromGM[0], avgAndStDevFromGM[1], minSum);
          
            double min, max, binSz;
            double[] ww = new double[vl.length];
            Arrays.fill(ww, 0.5);
            UnivariateDistance.DCov dcov;
            UnivariateDistance.DCor dcor;
            
            double[] vl1 = Arrays.copyOf(vl, vl.length/2);
            double[] vl2 = Arrays.copyOfRange(vl, vl.length/2, vl.length);
            dcor = UnivariateDistance.fastDcor(vl1, vl2);
            System.out.println("correlation within points=" + Math.sqrt(dcor.corSq));
            
            dcov = UnivariateDistance.fastDcov(vl1, vl2);
            System.out.println("cov:" + dcov.toString2());
            min = MiscMath0.findMin(dcov.ai);
            max = MiscMath0.findMax(dcov.ai);
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.ai, binSz, min, max);
            str = h.plotHistogram("[x diff sums],onSphere=" + onSurface, "dist_hist_1");
            
            min = MiscMath0.findMin(dcov.bi);
            max = MiscMath0.findMax(dcov.bi);
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.bi, binSz, min, max);
            str = h.plotHistogram("[y diff sums],onSphere=" + onSurface, "dist_hist_2");
            
            min = MiscMath0.findMin(dcov.iv1);
            max = MiscMath0.findMax(dcov.iv1);
            double[] avgAndStDev2 = MiscMath0.getAvgAndStDev(dcov.iv1);
            double[] qs = MiscMath0.calcMedianAndIQR(dcov.iv1);
            min = qs[0] - 1.5*qs[1];
            max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv1, binSz, min, max);
            str = h.plotHistogram("[iv1=cumulative 1],onSphere=" + onSurface, "dist_hist_iv1");
            
            qs = MiscMath0.calcMedianAndIQR(dcov.iv2);
            min = qs[0] - 1.5*qs[1];
            max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv2, binSz, min, max);
            str = h.plotHistogram("[iv2=cumulative X],onSphere=" + onSurface, "dist_hist_iv2");
            
            qs = MiscMath0.calcMedianAndIQR(dcov.iv3);
            min = qs[0] - 1.5*qs[1];
            max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv3, binSz, min, max);
            str = h.plotHistogram("[iv3=cumulative Y],onSphere=" + onSurface, "dist_hist_iv3");
            
            qs = MiscMath0.calcMedianAndIQR(dcov.iv4);
            min = qs[0] - 1.5*qs[1];
            max = qs[0] + 1.5*qs[1];
            binSz = (max - min)/20.;
            h = Histogram.createSimpleHistogram(dcov.iv4, binSz, min, max);
            str = h.plotHistogram("[iv4=cumulative X*Y],onSphere=" + onSurface, "dist_hist_iv4");
            
            int z = 0;
        /*
        "Multivariate tests of uniformity"
        2015, Yang & Modarres, DOI 10.1007/s00362-015-0715-x
        
        They compare to unit hyper-cube S = [0, 1]^d support.
        */
        
        
        // see https://www.researchgate.net/publication/283875292_Multivariate_tests_of_uniformity
        
        //M^2, Q3, L_MST
        // (MS2) Multivariate Hegazy–Green tests using spherical (MS2)
        // (ML2) lens  
        // (MT2) random Tukey depth functions, 
        // (LMST) total MST weight , 
        // and interpoint distances of uniforms (Q1, Q2, and Q3),
        
        //  Ripley's K function
        
        // Kuiper test for the hypothesis of uniformity?
        
        }
    }
    
}
