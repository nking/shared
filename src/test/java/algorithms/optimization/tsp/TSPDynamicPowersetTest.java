package algorithms.optimization.tsp;

import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import gnu.trove.list.TIntList;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TIntArrayList;
import java.math.BigInteger;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TSPDynamicPowersetTest extends TestCase {

    public TSPDynamicPowersetTest(String testName) {
        super(testName);
    }

    public void test0() throws InterruptedException {

        int n2 = 6;
        double[][] dist2 = new double[n2][n2];
        for (double[] row : dist2) {
            Arrays.fill(row, 10000);
        }
        for (int i = 0; i < n2; ++i) {
            dist2[i][i] = 0;
        }
        dist2[5][0] = 10;
        dist2[1][5] = 12;
        dist2[4][1] = 2;
        dist2[2][4] = 4;
        dist2[3][2] = 6;
        dist2[0][3] = 8;
        
        ////[3, 2, 4] =  100010011 = 275, sum = 10
        
        double expectedCost = 42; //8 + 6 + 4 + 2 + 12 + 10
        int[] expectedTour0 = new int[]{0, 3, 2, 4, 1, 5, 0};
        int[] expectedTour1 = new int[]{1, 5, 0, 3, 2, 4, 1};

        int[] tour = TSPDynamicPowerset.findMinCycle(dist2);

        assertTrue(Arrays.equals(expectedTour0, tour) || Arrays.equals(expectedTour1, tour));

    }

    public void test1() throws InterruptedException {

        int n2 = 4;
        double[][] x = new double[n2][];
        x[0] = new double[]{0, 0};
        x[1] = new double[]{0, 1};
        x[2] = new double[]{2, 0};
        x[3] = new double[]{3, 1};
        
        // optimal is 0, 1, 3, 2, 0 = 7.41
        // greedy is 0, 1, 2, 3, 0  = 7.81
        int i, j;
        double[] xi, xj;
        double xd, yd;
        double[][] dist = new double[n2][];
        for (i = 0; i < n2; ++i) {
            dist[i] = new double[n2];
            xi = x[i];
            for (j = 0; j < n2; ++j) {
                xj = x[j];
                xd = (xi[0] - xj[0]);
                yd = (xi[1] - xj[1]);
                dist[i][j] = Math.sqrt(xd*xd + yd*yd);
            }
        }

        System.out.printf("dist=\n%s\n", FormatArray.toString(dist, "%.4f"));
        
        int[] expectedGreedyFail = new int[]{0, 1, 2, 3, 0};
        double expectedOptimalCost = 7.41;
        double expectedGreedyCost = 7.81;
        
        double expectedCost = expectedOptimalCost;
        int[] expectedTour0 = new int[]{0, 2, 3, 1, 0};
        int[] expectedTour1 = new int[]{0, 1, 3, 2, 0};
        // or any cycling of that

        double sum0 = 0;
        double sum1 = 0;
        for (i = 1; i < expectedTour0.length; ++i) {
            sum0 += dist[expectedTour0[i-1]][expectedTour0[i]];
            sum1 += dist[expectedTour1[i-1]][expectedTour1[i]];
        }

        int[] tour = TSPDynamicPowerset.findMinCycle(dist);

        assertTrue(Arrays.equals(expectedTour0, tour) || Arrays.equals(expectedTour1, tour));

    }
  
    public void test3() throws InterruptedException {

        int n2 = 5;
        double[][] dist = new double[n2][n2];
        int zed = 0;
        dist[0] = new double[]{zed, 12, 10, 19, 8};
        dist[1] = new double[]{12, zed, 3, 7, 6};
        dist[2] = new double[]{10, 3, zed, 2, 20};
        dist[3] = new double[]{19, 7, 2, zed, 4};
        dist[4] = new double[]{8, 6, 20, 4, zed};
                
        //A B C D E
        //0 1 2 3 4
        //A, E, D, C, B -> 0, 4, 3, 2, 1, 0
        int[] expectedTour0 = new int[]{0, 4, 3, 2, 1, 0};
        int[] expectedTour1 = new int[]{1, 2, 3, 4, 0, 1};
        int[] expectedTour;
        double expectedCost = 29;//8+4+2+3+12

        int[] tour = TSPDynamicPowerset.findMinCycle(dist);

        assertTrue(Arrays.equals(expectedTour0, tour) || Arrays.equals(expectedTour1, tour));

    }
}
