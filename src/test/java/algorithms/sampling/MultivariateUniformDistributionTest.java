package algorithms.sampling;

import algorithms.correlation.DistanceTest;
import algorithms.misc.MiscMath0;
import algorithms.misc.MiscSorter;
import algorithms.misc.Standardization;
import algorithms.optimization.GeometricMedian;
import algorithms.util.FormatArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;
import junit.framework.TestCase;
import thirdparty.dlib.optimization.AbstractGeometricMedianFunction;
import thirdparty.dlib.optimization.GeometricMedianUnweightedFunction;
import thirdparty.dlib.optimization.GeometricMedianWeightedFunction;

/**
 *
 * @author nichole
 */
public class MultivariateUniformDistributionTest extends TestCase {
    
    public MultivariateUniformDistributionTest(String testName) {
        super(testName);
    }
    
    public void test0() throws NoSuchAlgorithmException {
        
        System.out.println("test0");
        
        //mean == 0
        // use UnivariateDistance w/ expected distr?  then chisqstat?
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        int nTests = 1;
        int nDimensions = 5;
        int nSamples = 10;
        boolean onSphere = true;
        
        double[][] v = new double[nSamples][nDimensions];
        double[] vl = new double[nSamples*nDimensions];
            
        for (int ii = 0; ii < nTests; ++ii) {

            for (int i = 0; i < nSamples; ++i) {
                // a quick look at some of the stats:
                v[i] = MultivariateUniformDistribution.generateUnitStandard(
                    nDimensions, rand, onSphere);
                double[] avgAndStDev = MiscMath0.getAvgAndStDev(v[i]);
                System.out.println("\nv[" + i + "]=" + FormatArray.toString(v[i], "%11.3f"));
                System.out.println("   mean, stDev=" + FormatArray.toString(avgAndStDev, "%11.3f"));
                System.arraycopy(v[i], 0, vl, i*nDimensions, nDimensions);
            }
            double[] avgAndStDev = MiscMath0.getAvgAndStDev(vl);
            System.out.println("all: mean, stDev=" + FormatArray.toString(avgAndStDev, "%11.3f"));

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

            System.out.println("g.m.=" + FormatArray.toString(init, "%11.3f"));
            System.out.printf("mean diff = %11.3f, stdev diff = %11.3f, minSum=%11.3f\n",
                    avgAndStDevFromGM[0], avgAndStDevFromGM[1], minSum);
          

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
