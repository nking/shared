package algorithms.statistics;

import algorithms.misc.MiscMath0;
import static algorithms.statistics.HypersphereChordLength.chooseMCalcPairwiseDistances;
import algorithms.util.FormatArray;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.TDoubleSet;
import gnu.trove.set.hash.TDoubleHashSet;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class HypersphereChordLengthTest extends TestCase {
    
    public HypersphereChordLengthTest(String testName) {
        super(testName);
    }
    
    public void testMean() {
        double m, diff;
        double tol = 2e-3;
        
        int[] n;
        double[] expected;
        
        double r;
        double[] rs = new double[]{1, 3};
        n = new int[] {2, 3, 4, 5, 6};
        
        for (int j = 0; j < rs.length; ++j) {
            r = rs[j];
            expected = new double[] {r*4./Math.PI, r*4./3., r*1.358, r*1.371, r*1.38};
            for (int i = 0; i < n.length; ++i) {
                m = HypersphereChordLength.meanOfChordLengthDistribution(r, n[i]);
                diff = Math.abs(m - expected[i]);
                assertTrue(diff < tol);
            }
        }
    }
    
    public void testVariance() {
        double v, diff;
        double tol = 1e-2;
        
        int[] n;
        double[] expected;
        
        double r, rsq;
        double[] rs = new double[]{1, 3};
        n = new int[] {2, 3, 4, 5, 6};
        
        for (int j = 0; j < rs.length; ++j) {
            r = rs[j];
            rsq = r*r;
            expected = new double[] {rsq * 0.379, rsq * 0.222, rsq * 0.156, 
                rsq * 0.119, rsq * 0.0956};
            for (int i = 0; i < n.length; ++i) {
                v = HypersphereChordLength.varianceOfChordLengthDistribution(r, n[i]);
                diff = Math.abs(v - expected[i]);
                assertTrue(diff < tol);
            }
        }
        
    }
   
    public void testPDF() {
        
        double r, d, diff;
        int[] rs = new int[] {1};
        double[] ds, pdf;
        double tol = 1e-2;
        
        int[] ns = new int[]{2, 3, 4, 5, 6};
        ns = new int[]{3};
        
        for (int j = 0; j < rs.length; ++j) {
            r = rs[j];
            ds = new double[]{0.25*r, 0.5*r, 0.75*r, r, 
                1.25*r, 1.5*r, 1.75*r, 2.0*r};
            for (int i = 0; i < ns.length; ++i) {
                pdf = HypersphereChordLength.pdf(ds, r, ns[i]);
                System.out.println("pdf: " + FormatArray.toString(pdf, "%.3f"));
                
                if (i == 1) {
                    diff = Math.abs(pdf[3] - 0.5);
                    assertTrue(diff < tol);
                    
                    diff = Math.abs(pdf[7] - 1.0);
                    assertTrue(diff < tol);
                } else if (i == 2) {
                    int maxIdx = MiscMath0.findYMaxIndex(pdf);
                    assertEquals(5, maxIdx);
                    assertTrue(pdf[5] > 0.9 && pdf[5] < 1.);
                    assertTrue(pdf[3] > 0.5 && pdf[3] < 0.6);
                    assertEquals(0., pdf[7]);
                } else if (i == 3) {
                    int maxIdx = MiscMath0.findYMaxIndex(pdf);
                    assertEquals(5, maxIdx);
                    assertTrue(pdf[5] > 1.0 && pdf[5] < 1.25);
                    assertTrue(pdf[3] > 0.5 && pdf[3] < 0.6);
                    assertEquals(0., pdf[7]);
                } else if (i == 4) {
                    int maxIdx = MiscMath0.findYMaxIndex(pdf);
                    assertEquals(5, maxIdx);
                    assertTrue(pdf[5] > 1.0 && pdf[5] < 1.25);
                    assertTrue(pdf[3] > 0.5 && pdf[3] < 0.6);
                    assertEquals(0., pdf[7]);
                }
            }
        }
    }
    
    public void testCDF() {
        
        double r, d, diff;
        int[] rs = new int[] {1};
        double[] ds, cdf;
        double tol = 1e-2;
        
        int[] ns = new int[]{2, 3, 4, 5, 6};
        
        for (int j = 0; j < rs.length; ++j) {
            r = rs[j];
            ds = new double[]{0.25*r, 0.5*r, 0.75*r, r, 
                1.25*r, 1.5*r, 1.75*r, 2.0*r};
            for (int i = 0; i < ns.length; ++i) {
                cdf = HypersphereChordLength.cdf(ds, r, ns[i]);
                //System.out.println(FormatArray.toString(cdf, "%.9f"));
                
                diff = Math.abs(cdf[7] - 1.0);
                assertTrue(diff < tol);
                if (i == 0) {
                    diff = Math.abs(cdf[3] - (1./3));
                    assertTrue(diff < tol);
                    assertTrue(cdf[0] < 0.1 && cdf[0] > 0.);                    
                } else if (i == 1) {
                    diff = Math.abs(cdf[3] - (1./4));
                    assertTrue(diff < tol);
                    assertTrue(cdf[0] < 0.2 && cdf[0] > 0.);
                } else if (i == 2) {
                    diff = Math.abs(cdf[3] - (1./5));
                    assertTrue(diff < 0.01);
                    assertTrue(cdf[0] < 0.004 && cdf[0] > 0.);
                } else if (i == 3) {
                    diff = Math.abs(cdf[3] - (1./6));
                    assertTrue(diff < 0.1);
                    assertTrue(cdf[0] < 0.002 && cdf[0] > 0.);
                } else if (i == 4) {
                    diff = Math.abs(cdf[3] - (1./7));
                    assertTrue(diff < 0.1);
                    assertTrue(cdf[0] < 0.0002 && cdf[0] > 0.);
                }
            }
        }
    }
    
    public void testFindAlpha() {
        
        // quick look at finding critical values using the CDF.
        // no inverse function, so "trial-and-error".
        //    For r = 1.  n = [2:10:+1,15:50:+5,60:100:+10]
        //       find d's where alpha=0.95
        
        double r = 1;
        TIntList ns = new TIntArrayList();
        
        int i, j;
        double tol = 1e-2;
        
        for (i = 2; i < 10; ++i) {
            ns.add(i);
        }
        for (i = 10; i < 50; i+=5) {
            ns.add(i);
        }
        for (i = 50; i <= 100; i+=10) {
            ns.add(i);
        }
        
        // array to hold the chord lengths where alpha~0.95 in CDF
        TDoubleList da = new TDoubleArrayList();
        
        double[] ds = new double[]{ 1.5*r, 1.75*r,
            1.76*r, 1.77*r, 1.78*r, 1.79*r, 
            1.80*r, 1.81*r, 1.82*r, 1.83*r, 1.84*r, 1.85*r, 1.86*r, 1.87*r, 1.88*r, 1.89*r,
            1.90*r, 1.91*r, 1.92*r, 1.93*r, 1.94*r, 1.95*r, 1.96*r, 1.97*r, 1.98*r, 1.99*r,
            2.0*r};
        
        int n, idx;
        double[] cdf;
        
        for (j = 0; j < ns.size(); ++j) {
            
            n = ns.get(j);
            
            // binary search between d=1.5*r and d=2*R for n <= 100
            cdf = HypersphereChordLength.cdf(ds, r, n);
            System.out.println(FormatArray.toString(cdf, "%.9f"));
            
            // Type I error rejects a null hypothesis that is actually true.
            // Type II error accepts a null hypothesis that is actually false.
            
            // alpha=0.05 (probability of Type I error)
            // confidence level or a confidence coefficient, (1 - α)100% = 95%
            idx = CDFRandomSelect.binarySearchForNearest(cdf, 0.95, tol);
            // confidence interval is interval in x capturing 95% of area 
            //     under curve, e.g. mu +- 2*sigma/sqrt(n)

            da.add(ds[idx]);
        }
        
        for (j = 0; j < ns.size(); ++j) {
            System.out.printf("n=%d d(1-alpha=0.095)~%.2f\n", ns.get(j), da.get(j));
        }
    }
    
    public void testChooseMCalcPairwiseDistances() throws NoSuchAlgorithmException {
        int nPoints = 3;
        int nDimensions = 2;
        int m = 3;
        
        double[][] x = new double[nPoints][nDimensions];
        x[0] = new double[]{0.5, 0};
        x[1] = new double[]{1, 0.5};
        x[2] = new double[]{0.5, 1.0};
        
        double[] expectedD = new double[3];
        expectedD[0] = Math.sqrt(0.5);
        expectedD[1] = 1;
        expectedD[2] = Math.sqrt(0.5);
        
        TDoubleSet expectedDSet = new TDoubleHashSet();
        for (int i = 0; i < expectedD.length; ++i) {
            expectedDSet.add(expectedD[i]);
        }
        
        int[] xMIdx = HypersphereChordLength.chooseM(m, x.length);

        double[] d = HypersphereChordLength.chooseMCalcPairwiseDistances(x, xMIdx);
        
        assertEquals(3, d.length);
        double tol = 1.e-5;
        double diff;
        for (int i = 0; i < d.length; ++i) {
            for (int j = 0; j < expectedD.length; ++j) {
                diff = Math.abs(d[i] - expectedD[j]);
                if (diff < tol) {
                    expectedDSet.remove(expectedD[j]);
                    break;
                }
            }
        }
        assertEquals(0, expectedDSet.size());
        
    }
    
    public void testCalcL1UniformityStatistic() throws NoSuchAlgorithmException {
        
        double d, diff;
        int r = 1;
        double[] ds, pdf;
        double tol = 1e-2;
        
        int ns = 4;
        int nPoints = 250;
        boolean onSurfae = true;
        int m = 100; 
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        double[][] x = MultivariateUniformDistribution
            .generateUnitStandardNSphereWithRejection(ns, nPoints, rand, onSurfae);
        
        double l1Sum = HypersphereChordLength.calcL1UniformityStatistic(x, m);
        
        System.out.printf("l1Sum = %.4e\n", l1Sum);
    }
    
}
