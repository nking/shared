package algorithms.optimization.tsp;

import algorithms.misc.MiscMath0;
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
public class TSPHybridDynamicBruteForceTest extends TestCase {

    private int n = 6;
    private TSPHybridDynamicBruteForce tsp;

    public TSPHybridDynamicBruteForceTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        double[][] dist = new double[n][];
        for (int i = 0; i < n; ++i) {
            dist[i] = new double[n];
        }
        // set 0:4 and 3:0 in dist to 412 and 53:
        dist[0][4] = 412;
        dist[3][0] = 53;

        tsp = new TSPHybridDynamicBruteForce(dist);
    }

    public void testSetBits() {
        // for n=6, there are 5 nodes that aren't the startnode 0: 1-5
        // w = 3
        //4 = 100
        //1 = 001
        //2 = 010
        //5 = 101
        //3 = 011
        // path of 4,1,2 = 010 001 100

        //path of 4,1,2,5,3 =  011  101 010 001 100

        int nSet, nUnset;
        int base10Node = 4;
        long path = 0;
        int pathNodeNumber = 0;
        long path4 = tsp.setBits(base10Node, path, pathNodeNumber);
        assertEquals(path4, 4);

        int b10 = tsp.getBase10NodeIndex(pathNodeNumber, path4);
        assertEquals(base10Node, b10);

        assertEquals(n-1, tsp.numberOfUnsetNodes(path));
        assertEquals(0, tsp.numberOfSetNodes(path));

        assertEquals(n-2, tsp.numberOfUnsetNodes(path4));
        assertEquals(1, tsp.numberOfSetNodes(path4));


        base10Node = 1;
        pathNodeNumber++;
        long path41 = tsp.setBits(base10Node, path4, pathNodeNumber);
        assertEquals(path41, 12);

        b10 = tsp.getBase10NodeIndex(pathNodeNumber, path41);
        assertEquals(base10Node, b10);

        base10Node = 2;
        pathNodeNumber++;
        long path412 = tsp.setBits(base10Node, path41, pathNodeNumber);
        assertEquals(path412, 140);
        
        TIntList pathNodes = new TIntArrayList();
        tsp.readPathIntoBase10(path412, pathNodes);
        assertEquals(3, pathNodes.size());
        assertEquals(4, pathNodes.get(0));
        assertEquals(1, pathNodes.get(1));
        assertEquals(2, pathNodes.get(2));

        b10 = tsp.getBase10NodeIndex(pathNodeNumber, path412);
        assertEquals(base10Node, b10);

        int[] s412 = new int[]{4, 1, 2};
        long path412a = tsp.createAMemoNodeBitstring(s412);
        assertEquals(path412, path412a);//140

        long setNodes412Base10 =  tsp.findSetBitsBase10(path412);
        assertTrue((setNodes412Base10 & (1L << 1)) != 0);
        assertTrue((setNodes412Base10 & (1L << 2)) != 0);
        assertTrue((setNodes412Base10 & (1L << 4)) != 0);
        
        TIntList unsetNodes = new TIntArrayList();
        tsp.findUnsetBitsBase10(path412, unsetNodes);
        assertEquals(2, unsetNodes.size());
        assertEquals(3, unsetNodes.get(0));
        assertEquals(5, unsetNodes.get(1));

        long unsetNodes53 = tsp.findUnsetBitsBase10(path412);
        assertTrue( (unsetNodes53 & (1L << 3)) != 0);
        assertTrue( (unsetNodes53 & (1L << 5)) != 0);

        int nSet412 = tsp.numberOfSetNodes(path412);
        assertEquals(3, nSet412);

        long path412b = tsp.concatenate(0, 0, s412);
        assertEquals(path412, path412b);

        //path of 4,1,2,5,3 =  011 101 010 001 100
        int[] s53 = new int[]{5, 3};
        long path41253 = tsp.concatenate(path412, nSet412, s53);
        assertEquals(14988, path41253);

        // set 0:4 and 3:0 in dist to 412 and 53 in setUp:
        tsp.compareToMin(path41253, 0);
        double expected = 412 + 53;
        assertEquals(expected, tsp.getMinCost());

        //path of 4,1,2,5 =  101 010 001 100 = 2700
        long path4125 = tsp.createAMemoNodeBitstring(new int[]{4, 1, 2, 5});
        assertEquals(2700, path4125);
    }

    public void testCalculateAndStore3NodePaths() throws InterruptedException {

        int n2 = 6;
        double[][] dist2 = new double[n2][n2];
        for (double[] row : dist2) {
            java.util.Arrays.fill(row, 10000);
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
        
        tsp = new TSPHybridDynamicBruteForce(dist2);
        tsp.init3NodePaths();
        int sz = (int)MiscMath0.computeNDivNMinusK(dist2.length-1, 3);
        //assertEquals(sz, tsp.getMemoLength());
        
        //tsp.printMemo();
        tsp = new TSPHybridDynamicBruteForce(dist2);
        tsp.init4NodePaths();   
        sz = (int)MiscMath0.computeNDivNMinusK(dist2.length-1, 4);
        //assertEquals(sz, tsp.getMemoLength());
    }
    
    public void testInitNodePaths() throws InterruptedException {
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
        
        tsp = new TSPHybridDynamicBruteForce(dist);
        tsp.initNodePaths();
        
        int sze = (int)MiscMath0.factorial(dist.length-1);
        int sz = tsp.getMemoLength();
        assertEquals(sze, sze);
        
        tsp.printMemo();
    }
    
    public void testRecursive0() throws InterruptedException {

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
        
        tsp = new TSPHybridDynamicBruteForce(dist2);
        tsp.solveRecursively();
        double cost = tsp.getMinCost();
        TLongList pathsBitstrings = tsp.getMinPathBitstrings();
        
        System.out.printf("minCost=%.2f, expected=%.2f\n", cost, expectedCost);
                
        TIntList path0 = tsp.getMinPath(0);
        System.out.println(path0.toString());
        
        assertTrue(Math.abs(expectedCost - cost) < 0.01);
        assertTrue(Arrays.equals(expectedTour0, path0.toArray()));
    }

    public void testRecursive1() throws InterruptedException {

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
        
        int[] expectedGreedyFail = new int[]{0, 1, 2, 3, 0};
        double expectedOptimalCost = 7.41;
        double expectedGreedyCost = 7.81;
        
        double expectedCost = expectedOptimalCost;
        int[] expectedTour0 = new int[]{0, 2, 3, 1, 0};
        int[] expectedTour1 = new int[]{0, 1, 3, 2, 0};
        
        tsp = new TSPHybridDynamicBruteForce(dist);
        tsp.solveRecursively();
        double cost = tsp.getMinCost();
        TLongList pathsBitstrings = tsp.getMinPathBitstrings();
        
        System.out.printf("minCost=%.2f, expected=%.2f\n", cost, expectedCost);
                
        TIntList path0 = tsp.getMinPath(0);
        System.out.println(path0.toString());
        
        assertTrue(Math.abs(expectedCost - cost) < 0.01);
        assertTrue(Arrays.equals(expectedTour0, path0.toArray()));
    }
  
    public void testRecursive3() throws InterruptedException {

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
        
        tsp = new TSPHybridDynamicBruteForce(dist);
        tsp.solveRecursively();
        double cost = tsp.getMinCost();
        TLongList pathsBitstrings = tsp.getMinPathBitstrings();
        
        System.out.printf("minCost=%.2f, expected=%.2f\n", cost, expectedCost);
                
        TIntList path0 = tsp.getMinPath(0);
        
        System.out.println("tour=" + path0.toString());
        System.out.flush();
        
        assertTrue(Math.abs(expectedCost - cost) < 0.01);
        assertTrue(Arrays.equals(expectedTour0, path0.toArray()));
    }

    public void testIterative0() throws InterruptedException {

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
        
        tsp = new TSPHybridDynamicBruteForce(dist2);
        tsp.solveIteratively();
        double cost = tsp.getMinCost();
        TLongList pathsBitstrings = tsp.getMinPathBitstrings();
        
        System.out.printf("minCost=%.2f, expected=%.2f\n", cost, expectedCost);
                
        TIntList path0 = tsp.getMinPath(0);
        System.out.println(path0.toString());
        
        assertTrue(Math.abs(expectedCost - cost) < 0.01);
        assertTrue(Arrays.equals(expectedTour0, path0.toArray()));
    }

    public void testIterative1() throws InterruptedException {

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
        
        int[] expectedGreedyFail = new int[]{0, 1, 2, 3, 0};
        double expectedOptimalCost = 7.41;
        double expectedGreedyCost = 7.81;
        
        double expectedCost = expectedOptimalCost;
        int[] expectedTour0 = new int[]{0, 2, 3, 1, 0};
        int[] expectedTour1 = new int[]{0, 1, 3, 2, 0};
        
        tsp = new TSPHybridDynamicBruteForce(dist);
        tsp.solveIteratively();
        double cost = tsp.getMinCost();
        TLongList pathsBitstrings = tsp.getMinPathBitstrings();
        
        System.out.printf("minCost=%.2f, expected=%.2f\n", cost, expectedCost);
                
        TIntList path0 = tsp.getMinPath(0);
        System.out.println(path0.toString());
        
        assertTrue(Math.abs(expectedCost - cost) < 0.01);
        assertTrue(Arrays.equals(expectedTour0, path0.toArray()));
    }
  
    public void testIterative3() throws InterruptedException {

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
        
        tsp = new TSPHybridDynamicBruteForce(dist);
        tsp.solveIteratively();
        double cost = tsp.getMinCost();
        TLongList pathsBitstrings = tsp.getMinPathBitstrings();
        
        System.out.printf("minCost=%.2f, expected=%.2f\n", cost, expectedCost);
                
        TIntList path0 = tsp.getMinPath(0);
        
        System.out.println("tour=" + path0.toString());
        System.out.flush();
        
        assertTrue(Math.abs(expectedCost - cost) < 0.01);
        assertTrue(Arrays.equals(expectedTour0, path0.toArray()));
    }

    public void testCount() {
        
        int n = 8;
        BigInteger c;
        
        c = TSPHybridDynamicBruteForce.count2(n);
        System.out.printf("** dynamic: n=%d c=%s\n", n, c.toString());
        
        n=14;
        c = TSPHybridDynamicBruteForce.count2(n);
        System.out.printf("** dynamic: n=%d c=%s\n", n, c.toString());
        
        n=29;
        c = TSPHybridDynamicBruteForce.count2(n);
        System.out.printf("** dynamic: n=%d c=%s\n", n, c.toString());
        
        n=49;
        c = TSPHybridDynamicBruteForce.count2(n);
        System.out.printf("** dynamic: n=%d c=%s\n", n, c.toString());
        
        n=731;
        c = TSPHybridDynamicBruteForce.count2(n);
        System.out.printf("** dynamic: n=%d c=%s\n", n, c.toString());
        
    }
}
