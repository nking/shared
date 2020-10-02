package algorithms.sampling;

import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
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
        
        for (int j = 0; j < rs.length; ++j) {
            r = rs[j];
            ds = new double[]{0.25*r, 0.5*r, 0.75*r, r, 
                1.25*r, 1.5*r, 1.75*r, 2.0*r};
            for (int i = 0; i < ns.length; ++i) {
                pdf = HypersphereChordLength.pdf(ds, r, ns[i]);
                //System.out.println(FormatArray.toString(pdf, "%.3f"));
                
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
}
