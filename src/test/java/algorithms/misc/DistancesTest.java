package algorithms.misc;

import junit.framework.TestCase;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class DistancesTest extends TestCase {
    
    public DistancesTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
        double tol = 1e-5;
        
        double[] p1, p2;
        
        p1 = new double[]{Math.sqrt(2), Math.sqrt(2)};
        p2 = new double[]{Math.sqrt(2), 0};
        
        double dsq = Distances.calcEuclideanSquared(p1, p2);
        
        assertTrue(Math.abs(2. - dsq) < tol);
    }

    public void testMax() {
        int[][] points = new int[][]{
                {1,3}, {2,1}, {4,4}, {4,2}
        };
        long expMax = 5; // 1,2
        int[] expIdxs = new int[]{1,2};
        int[] idxs = Distances.maxManhattan(points);
        assertTrue(Arrays.equals(expIdxs, idxs));
    }

    public void testMax1() {
        int[][] points = new int[][]{
                {4,0}, {7,4}, {3,5}, {2,4}
        };
        long expMax = 7; // 0,1
        int[] expIdxs = new int[]{0,1};
        int[] idxs = Distances.maxManhattan(points);
        assertTrue(Arrays.equals(expIdxs, idxs));
    }
    public void testMax2() {
        int x0 = -4;
        int y0 = -4;
        int[][] points = new int[][]{
                {4+x0,0+y0}, {7+x0,4+y0}, {3+x0,5+y0}, {2+x0,4+y0}
        };
        long expMax = 7; // 0,1
        int[] expIdxs = new int[]{0,1};
        int[] idxs = Distances.maxManhattan(points);
        assertTrue(Arrays.equals(expIdxs, idxs));
    }
}
