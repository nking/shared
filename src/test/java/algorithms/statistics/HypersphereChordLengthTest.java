package algorithms.statistics;

import algorithms.misc.Distances;
import algorithms.misc.MiscMath0;
import algorithms.statistics.HypersphereChordLength.NonUniformityStats;
import algorithms.util.FormatArray;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.TDoubleSet;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TDoubleHashSet;
import gnu.trove.set.hash.TIntHashSet;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;
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
        
        TIntList ns = new TIntArrayList();
        
        int i, j;
        double tol = 1e-2;
        
        for (i = 2; i < 10; ++i) {
            ns.add(i);
        }
        for (i = 10; i < 50; i+=5) {
            ns.add(i);
        }
        for (i = 50; i < 100; i+=10) {
            ns.add(i);
        }
        for (i = 100; i <= 1000; i+=100) {
            ns.add(i);
        }
        for (i = 1000; i <= 5000; i+=500) {
            ns.add(i);
        }
        
        // array to hold the chord lengths where alpha~0.95 in CDF
        TDoubleList da = new TDoubleArrayList();
        
        int n;
        
        for (j = 0; j < ns.size(); ++j) {
            
            n = ns.get(j);

            da.add(HypersphereChordLength.findCVForAlpha95Percent(n));
        }
        
        for (j = 0; j < ns.size(); ++j) {
            System.out.printf("n=%d d(1-alpha=0.095)~%.2f%n", ns.get(j), da.get(j));
        }
    }
    
    public void testChooseMCalcPairwiseDistances() throws NoSuchAlgorithmException {
        int nPoints = 3;
        int nDimensions = 2;
        int m = 3;
        
         SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
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
        
        int[] xMIdx = HypersphereChordLength.chooseM(m, x.length, rand);

        double[] d = HypersphereChordLength.calcPairwiseDistances(x, xMIdx);
        
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
        
        int ns = 5;
        int nPoints = 250;
        boolean onSurface = true;
        int m = 100; 
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        double[][] x = MultivariateUniformDistribution
            .generateUnitStandardNSphereWithRejection(ns, nPoints, rand, onSurface);
        
        double l1Sum = HypersphereChordLength.calcL1UniformityStatistic(x, m, 
            HypersphereChordLength.POINT_DISTRIBUTION_TYPE.INTRA_DISTANCE_2, rand);
        
        System.out.printf("l1Sum = %.4e%n", l1Sum);
    }
    
    public void testCalcConfidenceOfNonUniformity_same() throws NoSuchAlgorithmException {
        
        System.out.println("testCalcConfidenceOfNonUniformity_same");

        int nTests = 10;
        
        int ns = 7;
        int nPoints = 250;
        boolean onSurface = true;
        int m = 50;//100; 
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        for (int ii = 0; ii < nTests; ++ii) {
        
            double[][] x = MultivariateUniformDistribution
                .generateUnitStandardNSphereWithRejection(ns, nPoints, rand, onSurface);

            NonUniformityStats stats = HypersphereChordLength.calcConfidenceOfNonUniformity(x, m, 
                HypersphereChordLength.POINT_DISTRIBUTION_TYPE.INTRA_DISTANCE_2, rand);

            System.out.println("stats=" + stats.toString() + "%n");

            assertFalse(stats.isConsistentWithNonUniform);
        }
    }
    
    public void testCalcConfidenceOfNonUniformity_different() throws NoSuchAlgorithmException {
        
        System.out.println("testCalcConfidenceOfNonUniformity_different");
                
        int nTests = 10;
        
        int ns = 7;
        int nPoints = 250;//400;//250;
        boolean onSurface = true;
        int m = 50;//100; 
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        for (int ii = 0; ii < nTests; ++ii) {
        
            double[][] x = makeNonUniform(ns, nPoints, rand, onSurface);

            NonUniformityStats stats = HypersphereChordLength.calcConfidenceOfNonUniformity(x, m, 
                HypersphereChordLength.POINT_DISTRIBUTION_TYPE.INTRA_DISTANCE_2, rand);

            //System.out.println("stats=" + stats.toString() + "%n");

            assertTrue(stats.isConsistentWithNonUniform);
        }
    }

    /*
    double[] u = new double[d];
                
        u = UnivariateNormalDistribution.randomSampleOfUnitStandard(rand, d);
        
        // then normalization of each is the sqrt(sum of squares of all)
        double norm = 0;
        for (double v : u) {
            norm += (v * v);
        }
        norm = 1./Math.sqrt(norm);
        
        if (!onSurface) {
            // to avoid bunching at the center of hypersphere
            double b = Math.pow(rand.nextDouble(), 1./d);
        
            norm *= b;
        }
        
        for (int i = 0; i < u.length; ++i) {
            u[i] *= norm;
        }
    */

    private double[][] makeNonUniform(int nDimensions, int nPoints, SecureRandom rand, boolean onSurface) {
        
        double[][] x = new double[nPoints][nDimensions];
        
        for (int i0 = 0; i0 < nPoints; ++i0) {
            x[i0] = new double[nDimensions];
            
            // distributed on a hypercube and projected to a hypersphere will have
            //   a density larger where corners were projected to the hypersphere
            double norm = 0;
            for (int j0 = 0; j0 < nDimensions; ++j0) {
                x[i0][j0] = rand.nextDouble();
                norm += (x[i0][j0] * x[i0][j0]);
            }
            norm = 1./Math.sqrt(norm);

            if (!onSurface) {
                // to avoid bunching at the center of hypersphere
                double b = Math.pow(rand.nextDouble(), 1./nDimensions);
                norm *= b;
            }

            for (int j0 = 0; j0 < nDimensions; ++j0) {
                x[i0][j0] *= norm;
            }
        }
        
        return x;
    }
    
    
}
