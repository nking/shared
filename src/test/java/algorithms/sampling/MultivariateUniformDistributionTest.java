package algorithms.sampling;

import algorithms.correlation.DistanceTest;
import algorithms.misc.MiscMath0;
import algorithms.misc.MiscSorter;
import algorithms.misc.Standardization;
import algorithms.optimization.GeometricMedian;
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
        
        int nDimensions = 5;
        boolean onSphere = true;
        
        // a quick look at some of the stats:
        
        double[] u = MultivariateUniformDistribution.generateUnitStandard(
            nDimensions, rand, onSphere);
        
        double[] w = new double[u.length];
        Arrays.fill(w, 1.0);
        
        double[] avgAndStDev = MiscMath0.getAvgAndStDev(u);
        
        GeometricMedian gm = new GeometricMedian();
        
        GeometricMedianWeightedFunction f = new GeometricMedianWeightedFunction(
            u, 1, w);
        double[] init = new double[]{0};
        double minSum = gm.newtonsThenVardiZhang(f, init);
        double gc = init[0];
        
        double[] diff = new double[u.length];
        for (int i = 0; i < u.length; ++i) {
            diff[i] = Math.abs(u[i] - gc);
        }
        
        int[] oIdx = MiscSorter.mergeSortIncreasing(diff);
        
        double[] avgAndStDevFromGM = MiscMath0.getAvgAndStDev(diff);
        
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
