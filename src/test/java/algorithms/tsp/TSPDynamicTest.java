package algorithms.tsp;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import gnu.trove.list.TIntList;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TSPDynamicTest extends TestCase {

    private int n = 6;
    private TSPDynamic tsp;

    public TSPDynamicTest(String testName) {
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

        tsp = new TSPDynamic(dist);
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
        
        tsp = new TSPDynamic(dist2);
        tsp.init3NodePaths();
        int sz = (int)MiscMath0.computeNDivNMinusK(dist2.length-1, 3);
        //assertEquals(sz, tsp.getMemoLength());
        
        //tsp.printMemo();
        tsp = new TSPDynamic(dist2);
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
        
        tsp = new TSPDynamic(dist);
        tsp.initNodePaths();
        
        int sze = (int)MiscMath0.factorial(dist.length-1);
        int sz = tsp.getMemoLength();
        assertEquals(sze, sze);
        
        tsp.printMemo();
    }
    
    public void test0() throws Exception {

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
        
        tsp = new TSPDynamic(dist2);
        tsp.solve();
        double cost = tsp.getMinCost();
        
        System.out.printf("minCost=%.2f, expected=%.2f\n", cost, expectedCost);
                
        TIntList path0 = tsp.getMinPath(0);
        System.out.println(path0.toString());
        
        assertTrue(Math.abs(expectedCost - cost) < 0.01);
        assertTrue(Arrays.equals(expectedTour0, path0.toArray()));
    }

    public void test() throws Exception {

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
        
        tsp = new TSPDynamic(dist);
        tsp.solve();
        double cost = tsp.getMinCost();
        
        System.out.printf("minCost=%.2f, expected=%.2f\n", cost, expectedCost);
                
        TIntList path0 = tsp.getMinPath(0);
        System.out.println(path0.toString());
        
        assertTrue(Math.abs(expectedCost - cost) < 0.01);
        assertTrue(Arrays.equals(expectedTour0, path0.toArray()));
    }
  
    public void test2() throws Exception {

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
        
        tsp = new TSPDynamic(dist);
        tsp.solve();
        double cost = tsp.getMinCost();
        
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
        
        c = TSPDynamic.count2(n);
        System.out.printf("** dynamic: n=%d c=%s\n", n, c.toString());
        
        n=14;
        c = TSPDynamic.count2(n);
        System.out.printf("** dynamic: n=%d c=%s\n", n, c.toString());
        
        n=29;
        c = TSPDynamic.count2(n);
        System.out.printf("** dynamic: n=%d c=%s\n", n, c.toString());
        
        n=49;
        c = TSPDynamic.count2(n);
        System.out.printf("** dynamic: n=%d c=%s\n", n, c.toString());
        
        n=731;
        c = TSPDynamic.count2(n);
        System.out.printf("** dynamic: n=%d c=%s\n", n, c.toString());
        
    }

    public void test10() throws Exception {

        /* Test of traveling salesman approximate tour. 
         * 
         * From Cormen et al. chap 35 Approximation Algorithms, 35.2 Traveling-Salesman problem, Fig 35.2
         * 
         * 
         *   6  |---|---|---|---|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   5  |---|--[A]--|--[D]--|---|---|----
         *      |   |   |   |   |   |   |   |
         *   4  |---|---|---|---|--[E]--|---|----
         *      |   |   |   |   |   |   |   |
         *   3  |---|--[B]--|--[F]--|--[G]--|----
         *      |   |   |   |   |   |   |   |
         *   2  |--[C]--|---|---|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   1  |---|---|--[H]--|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   0  |---|---|---|---|---|-------|----
         *      0   1   2   3   4   5   6   7

         a, b, c,  h, f, g, e, d
         0  1  2   7  5  6  4  3

         need 9 points to test part of the code, so will add
         a point called BC in between B and C that won't change
         the ideal path.

         a, b, bc, c,  h, f, g, e, d, a
         0, 1, 2,  3,  8, 6, 7, 5, 4, 0
        */
        int n2 = 9;
        double[][] x = new double[n2][];
        x[0] = new double[]{2, 5};
        x[1] = new double[]{2, 3};
        x[2] = new double[]{1.5, 2.5};
        x[3] = new double[]{1, 2};
        x[4] = new double[]{4, 5};
        x[5] = new double[]{5, 4};
        x[6] = new double[]{4, 3};
        x[7] = new double[]{6, 3};
        x[8] = new double[]{3, 1};

        double[][] dist = MatrixUtil.zeros(n2, n2);
        int i, j;
        double d, dx, dy, xi, yi, xj, yj;
        for (i = 0; i < n2; ++i) {
            xi = x[i][0];
            yi = x[i][1];
            for (j = 0; j < n2; ++j) {
                xj = x[j][0];
                yj = x[j][1];
                dx = xi - xj;
                dy = yi - yj;
                d = Math.sqrt(dx*dx + dy*dy);
                dist[i][j] = d;
            }
        }

        int[] expectedTour0 = new int[]{0, 1, 2,  3,  8, 6, 7, 5, 4, 0};
        int[] expectedTour1 = new int[]{0, 4, 5, 7, 6, 8, 3, 2, 1, 0};
        double expectedCost = 0;
        int idx0, idx1;
        for (i = 1; i < expectedTour1.length; ++i) {
            idx0 = expectedTour1[i-1];
            idx1 = expectedTour1[i];
            expectedCost += dist[idx0][idx1];
        }

        tsp = new TSPDynamic(dist);
        tsp.solve();
        double cost = tsp.getMinCost();

        System.out.printf("n=%d) minCost=%.2f, expected=%.2f\n", n2, cost, expectedCost);

        TIntList path0 = tsp.getMinPath(0);
        System.out.println(path0.toString());

        System.out.flush();

        assertTrue(Math.abs(expectedCost - cost) < 0.01);
        assertTrue(Arrays.equals(expectedTour1, path0.toArray()));
    }
    
    public void estATT48_15() throws Exception {
        
        TIntList expectedTour = new TIntArrayList();
        
        double[][] dist = MatrixUtil.zeros(15, 15);
        
        double expectedCost = populateATT48_15_cities(expectedTour, dist);
        
        tsp = new TSPDynamic(dist);
        tsp.solve();
        double cost = tsp.getMinCost();
        
        System.out.printf("minCost=%.2f, expected=%.2f\n", cost, expectedCost);
                
        TIntList path0 = tsp.getMinPath(0);
        
        System.out.println("tour=" + path0.toString());
        System.out.flush();
        
        assertTrue(Math.abs(expectedCost - cost) < 0.01);
        assertTrue(Arrays.equals(expectedTour.toArray(), path0.toArray()));
    }
    
    /**
     * 
     * @param expectedBestTour
     * @param dist
     * @return The minimal cost is 291
     * @throws Exception 
     */
    private double populateATT48_15_cities(TIntList expectedBestTour, 
        double dist[][]) throws Exception {

        String dir = 
            ResourceFinder.findTestResourcesDirectory()
            + "/att48";
            
        readBestTour(dir + "/p01_s.txt", expectedBestTour);
        
        readDistance(dir + "/p01_d.txt", dist);
        
        return 291.;
    }
    
    private void populateATT48_48_cities(List<PairInt> pointList, 
        TIntList expectedBestTour, 
        TIntObjectMap<TIntIntMap> adjCostMap) 
        throws Exception {

        String dir = 
            ResourceFinder.findTestResourcesDirectory()
            + "/att48";
            
        readCoords(dir + "/att48_xy.txt", pointList);
    
        readBestTour(dir + "/att48_s.txt", 
            expectedBestTour);
        
        readAdjacency(dir + "/att48_d.txt", adjCostMap);
    }

    private void readCoords(String filePath, 
        List<PairInt> pointList) throws Exception {
        
        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(new File(filePath)));
            
            String line = in.readLine();
            
            while (line != null) {
                
                String item1 = line.substring(0, 4).trim();
                
                String item2 = line.substring(5).trim();
                
                Integer x =
                    Integer.parseInt(item1);
                
                Integer y =
                    Integer.parseInt(item2);
                
                PairInt p = new PairInt(x.intValue(), y.intValue());
                
                pointList.add(p);
                
                line = in.readLine();
            }
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }   
    }

    /**
     * read expected best tour, transformed to starting index=0
     * (the file has city indexes that start with 1).
     * @param filePath
     * @param expectedBestTour
     * @throws Exception 
     */
    private void readBestTour(String filePath, 
        TIntList expectedBestTour) throws Exception {
        
        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(new File(filePath)));
            
            String line = in.readLine();
            
            while (line != null) {
                
                Integer index = Integer.parseInt(line.trim());
                
                // change to zero based indexes
                expectedBestTour.add(
                    index.intValue() - 1);
                
                line = in.readLine();
            }
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }
    }

    private void readAdjacency(String filePath, 
        TIntObjectMap<TIntIntMap> adjCostMap) 
        throws Exception {

        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(new File(filePath)));
            
            String line = in.readLine();
            
            int i = 0;
            
            while (line != null) {
            
                TIntIntMap idx2CostMap = new TIntIntHashMap();
                adjCostMap.put(i, idx2CostMap);
            
                String[] items = line.split("\\s+");
                if (items[0] == null || items[0].equals("")) {
                    items = Arrays.copyOfRange(items, 1, 
                        items.length);
                }
                
                for (int j = 0; j < items.length; ++j) {
                    Integer c = Integer.parseInt(items[j].trim());
                    if (!c.equals(Integer.valueOf(0))) {
                        idx2CostMap.put(j, c.intValue());
                    }
                }
                
                line = in.readLine();
                ++i;
            }
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }
    }
    
    private void readDistance(String filePath, double[][] dist) 
        throws Exception {

        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(new File(filePath)));
            
            String line = in.readLine();
            
            int i = 0;
            
            while (line != null) {
                        
                String[] items = line.split("\\s+");
                System.out.printf("items.length=%d", items.length);
                if (items[0] == null || items[0].equals("")) {
                    items = Arrays.copyOfRange(items, 1, 
                        items.length);
                }
                System.out.printf("->%d\n", items.length);
                
                for (int j = 0; j < items.length; ++j) {
                    Integer c = Integer.parseInt(items[j].trim());
                    dist[i][j] = c;
                }
                
                line = in.readLine();
                ++i;
            }
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }
    }
}
