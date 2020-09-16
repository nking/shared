package algorithms.correlation;

import algorithms.correlation.BruteForceDistance.DCOV;
import algorithms.correlation.UnivariateDistance.DCor;
import algorithms.correlation.UnivariateDistance.DCov;
import algorithms.misc.Misc0;
import java.util.Arrays;
import java.util.Random;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;
    
/**
 *
 * @author nichole
 */
public class DistanceTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public DistanceTest() {
    }
    
    public void testSort() {
        
        //double[][] _sortCheck(double[] x, double[] y)
        
        double eps = 1.e-17;
        
        double[] x = new double[]{1., 1.1, 1.1, 2., 3.};
        double[] y = new double[]{1., 1.1, 2.1, 3., 4.};
        
        double[] eX = Arrays.copyOf(x, x.length);
        double[] eY = Arrays.copyOf(y, y.length);
        
        double[][] sorted = UnivariateDistance._sortCheck(x, y);
        double diff;
        for (int i = 0; i < x.length; ++i) {
            diff = Math.abs(sorted[0][i] - eX[i]);
            assertTrue(diff < eps);
            diff = Math.abs(sorted[1][i] - eY[i]);
            assertTrue(diff < eps);
        }
        
        x = new double[]{2, 1., 1.1, 3, 1.1};
        y = new double[]{3, 1., 2.1, 4, 1.1};
        sorted = UnivariateDistance._sortCheck(x, y);
        for (int i = 0; i < x.length; ++i) {
            diff = Math.abs(sorted[0][i] - eX[i]);
            assertTrue(diff < eps);
            diff = Math.abs(sorted[1][i] - eY[i]);
            assertTrue(diff < eps);
        }
    }
    
    public void testSort2() {
        
        //double[][] _sortCheck(double[] x, double[] y)
        
        double eps = 1.e-17;
        
        //                        0  1  2  3  4  5  6  7
        double[] x = new double[]{3, 5, 7, 3, 8, 4, 6, 7};
        double[] y = new double[]{1, 5, 3, 2, 4, 6, 7, 5};
        
        //                         0  3  5  1  6  2  7  4
        double[] eX = new double[]{3, 3, 4, 5, 6, 7, 7, 8};
        double[] eY = new double[]{1, 2, 6, 5, 7, 3, 5, 4};
        
        double[][] sorted = UnivariateDistance._sortCheck(x, y);
        //System.out.println("sortedX=" + Arrays.toString(sorted[0]));
        //System.out.println("sortedY=" + Arrays.toString(sorted[1]));
        
        double diff;
        for (int i = 0; i < x.length; ++i) {
            diff = Math.abs(sorted[0][i] - eX[i]);
            assertTrue(diff < eps);
            diff = Math.abs(sorted[1][i] - eY[i]);
            assertTrue(diff < eps);
        }
        
    }
    
    public void testSort3() {
                
        double eps = 1.e-17;
        
        int n = 1200;
        int nIter = 1000;
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        //long seed = 1595442111027L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        double[] x = new double[n];
        double[] y = new double[n];
        
        for (int i = 0; i < nIter; ++i) {
            for (int j = 0; j < n; ++j) {
                y[j] = rand.nextDouble();
                if (j > 0 && ((j&1)==0) && rand.nextBoolean()) {
                    // occassionally, duplicate x's with presumably different y's
                    x[j] = x[j - 1];
                } else {
                    x[j] = rand.nextDouble();
                }
            }
            
            double[][] sorted = UnivariateDistance._sortCheck(x, y);
            //System.out.println("sortedX=" + Arrays.toString(sorted[0]));
            //System.out.println("sortedY=" + Arrays.toString(sorted[1]));

            double diff;
            for (int ii = 1; ii < x.length; ++ii) {
                diff = Math.abs(sorted[0][ii] - sorted[0][ii-1]);
                if (diff < eps) {
                    assertTrue(sorted[1][ii] >= sorted[1][ii-1]);
                } else {
                    assertTrue(sorted[0][ii] >= sorted[0][ii-1]);
                }
            }
        }        
    }
    
    public void testCovariance() {
        
        //double[][] _sortCheck(double[] x, double[] y)
        
        double eps = 1.e-17;
        
        //                        0  1  2  3  4  5  6  7
        double[] x = new double[]{3, 5, 7, 3, 8, 4, 6, 7};
        double[] y = new double[]{1, 5, 3, 2, 4, 6, 7, 5};
        
        //                         0  3  5  1  6  2  7  4
        double[] eX = new double[]{3, 3, 4, 5, 6, 7, 7, 8};
        double[] eY = new double[]{1, 2, 6, 5, 7, 3, 5, 4};
        
        DCor dCor2 = UnivariateDistance.fastDcor(x, y);
        System.out.println("fastDcor(X,Y)\n: " + dCor2.toString());
        
        double diff = 0;
        /*for (int i = 0; i < x.length; ++i) {
            diff = Math.abs(sorted[0][i] - eX[i]);
            assertTrue(diff < eps);
            diff = Math.abs(sorted[1][i] - eY[i]);
            assertTrue(diff < eps);
        }
        
        x = new double[]{2, 1., 1.1, 3, 1.1};
        y = new double[]{3, 1., 2.1, 4, 1.1};
        sorted = UnivariateDistance._sortCheck(x, y);
        for (int i = 0; i < x.length; ++i) {
            diff = Math.abs(sorted[0][i] - eX[i]);
            assertTrue(diff < eps);
            diff = Math.abs(sorted[1][i] - eY[i]);
            assertTrue(diff < eps);
        }*/
    }
    
    public void testBruteForceCovariance() {
        
        //double[][] _sortCheck(double[] x, double[] y)
        
        double eps = 1.e-17;
        
        //                        0  1  2  3  4  5  6  7
        double[] x = new double[]{3, 5, 7, 3, 8, 4, 6, 7};
        double[] y = new double[]{1, 5, 3, 2, 4, 6, 7, 5};
        
        //                         0  3  5  1  6  2  7  4
        double[] eX = new double[]{8, 7, 7, 6, 5, 4, 3, 3};
        double[] eY = new double[]{4.0, 3.0, 5.0, 7.0, 5.0, 6.0, 1.0, 2.0};
        double[] eD = new double[]{11, 6, 21, 14, 1, 3, 0, 0};
        
        DCov dcov = UnivariateDistance._univariateCovariance(x, y);
        double diff = 0;
        for (int i = 0; i < x.length; ++i) {
            diff = Math.abs(dcov.sortedX[i] - eX[i]);
            assertTrue(diff < eps);
            diff = Math.abs(dcov.sortedY[i] - eY[i]);
            assertTrue(diff < eps);
            diff = Math.abs(dcov.dcov[i] - eD[i]);
            assertTrue(diff < eps);
        }
        
        System.out.printf("\n_univariateCovariance:\n");
        System.out.println("XY:" + dcov.toString());
        
        double[][] a = new double[x.length][2];
        for (int i = 0; i < x.length; ++i) {
            a[i] = new double[]{dcov.sortedX[i], dcov.sortedY[i]};
        }
        
        {
            DCov dcov2 = UnivariateDistance._univariateCovariance(x, x);
            System.out.println("XX:" + dcov2.toString());
            dcov2 = UnivariateDistance._univariateCovariance(y, y);
            System.out.println("YY:" + dcov2.toString());
        }
        
        a = new double[x.length][1];
        double[][] b = new double[x.length][1];
        for (int i = 0; i < x.length; ++i) {
            a[i] = new double[]{dcov.sortedX[i]};
            b[i] = new double[]{dcov.sortedY[i]};
        }
        
        DCOV dc2 = BruteForceDistance.correlation1(a, b);
        StringBuilder sb = new StringBuilder();
            sb.append("bf dist. correlation:\n ");
            sb.append("  dVarX=").append(dc2.dVarX).append("\n");
            sb.append("  dVarY=").append(dc2.dVarY).append("\n");
            sb.append("  dCovXY=").append(dc2.dCov).append("\n");
            sb.append("  dCor=").append(dc2.cor).append("\n");
        System.out.println(sb.toString());
        
        System.out.flush();
    }
}
