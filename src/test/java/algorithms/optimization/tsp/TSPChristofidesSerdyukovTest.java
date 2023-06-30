package algorithms.optimization.tsp;

import algorithms.graphs.HierholzersEulerCircuit;
import algorithms.matrix.MatrixUtil;
import algorithms.optimization.tsp.TSPChristofidesSerdyukov;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TSPChristofidesSerdyukovTest extends TestCase {
    
    public TSPChristofidesSerdyukovTest(String testName) {
        super(testName);
    }

    public void testApproxTSPTour0() {
        System.out.println("testApproxTSPTour0");
        
        // using example from  
        //  https://en.m.wikipedia.org/wiki/Christofides_algorithm
        // edges obey triangle inequality
        int nVertexes = 5;
        int i, j;
        
        TIntObjectMap<TIntIntMap> adjCostMap = new TIntObjectHashMap<TIntIntMap>();        
        for (i = 0; i < nVertexes; ++i) {            
            adjCostMap.put(i, new TIntIntHashMap());
        }
        adjCostMap.get(0).put(1, 1);
        adjCostMap.get(0).put(2, 1);
        adjCostMap.get(0).put(3, 1);
        adjCostMap.get(0).put(4, 1);
        adjCostMap.get(1).put(0, 1);
        adjCostMap.get(1).put(2, 1);
        adjCostMap.get(1).put(3, 1);
        adjCostMap.get(1).put(4, 2);
        adjCostMap.get(2).put(0, 1);
        adjCostMap.get(2).put(1, 1);
        adjCostMap.get(2).put(3, 2);
        adjCostMap.get(2).put(4, 1);
        adjCostMap.get(3).put(0, 1);
        adjCostMap.get(3).put(1, 1);
        adjCostMap.get(3).put(2, 2);
        adjCostMap.get(3).put(4, 1);
        adjCostMap.get(4).put(0, 1);
        adjCostMap.get(4).put(1, 2);
        adjCostMap.get(4).put(2, 1);
        adjCostMap.get(4).put(3, 1);
        
        TSPChristofidesSerdyukov tsp = new TSPChristofidesSerdyukov();
        
        // MST 0:4, 0:3, 0:2, 0:1
        Map<Integer, LinkedList<Integer>> mstTree = tsp.buildMST(adjCostMap);
        assertEquals(1, mstTree.size());
        TIntSet expected = new TIntHashSet(new int[]{1,2,3,4});
        Iterator<Integer> intIter = mstTree.keySet().iterator();
        Iterator<Integer> intIter2;
        int u, v;
        while (intIter.hasNext()) {
            u = intIter.next();
            assertEquals(0, u);
            intIter2 = mstTree.get(u).iterator();
            while (intIter2.hasNext()) {
                v = intIter2.next();
                expected.remove(v);
            }
        }
        assertEquals(0, expected.size());
        
        //           0  1  2  3  4
        // degrees = 4  1  1  1  1
        int[] degrees = tsp.calculateDegrees(mstTree, nVertexes);
        int[] oddDVertexes = tsp.oddPassFilter(degrees);
        // there are an even number of odd vertexes
        assert((oddDVertexes.length & 1) != 1);
        expected = new TIntHashSet(new int[]{1,2,3,4});
        assertEquals(expected.size(), oddDVertexes.length);
        for (i = 0; i < oddDVertexes.length; ++i) {
            assertTrue(expected.contains(oddDVertexes[i]));
            expected.remove(oddDVertexes[i]);
        }
        assertEquals(0, expected.size());
        
        // subgraph of G from odd vertices:
        // 1:2, 1:3, 1:4(2)
        // 2:1, 2:3(2), 2:4
        // 3:1, 3:2(2), 3:4
        // 4:1(2), 4:2, 4:3
       
        // matching
        // 1:2, 3:4 (in context of oddDegress array, will see 0:1, 2:3)
        int[][] m = tsp.bipartiteMinCostMatchingFromSubgraph(
            oddDVertexes, adjCostMap);
        assertEquals(2, m.length);
        assertTrue((m[0][0]==1 && m[0][1]==2) || (m[1][0]==1 && m[1][1]==2));
        assertTrue((m[0][0]==3 && m[0][1]==4) || (m[1][0]==3 && m[1][1]==4));
        
        // union of MST and matching
        //MST 0:4, 0:3, 0:2, 0:1
        // 1:2, 3:4
        TIntObjectMap<TIntSet> h = tsp.unionMSTAndAssignments(mstTree, m);
        TIntObjectIterator<TIntSet> iterH = h.iterator();
        while (iterH.hasNext()) {
            iterH.advance();
            int node = iterH.key();
            TIntSet links = iterH.value();
           // System.out.printf("TunionM: node1=%d, links=%s\n", node, Arrays.toString(links.toArray()));
            assertTrue(node == 0 || node == 1 || node == 3);
            if (node == 0) {
                assertTrue(links.containsAll(new int[]{1,2,3,4}));
            } else if (node == 1) {
                assertTrue(links.containsAll(new int[]{2}));
            } else if (node == 3) {
                assertTrue(links.containsAll(new int[]{4}));
            }
        }
        assertEquals(3, h.size());
        
        
        // euler tour:
        //        0, 1, 2, 0, 3, 4, 0
        //?? ec: [0, 1, 2, 2, 3, 4, 4]
        HierholzersEulerCircuit hec = new HierholzersEulerCircuit();
        int[] ec = hec.createCircuit(h); // assuming start node = 0.
        //System.out.printf("ec: %s\n", Arrays.toString(ec));
        assertEquals(7, ec.length);
        assertEquals(0, ec[0]); assertEquals(1, ec[1]); assertEquals(2, ec[2]);
        assertEquals(3, ec[4]); assertEquals(4, ec[5]);
        
        // hamiltonian
        // 0, 1, 2, 3, 4, 0
        TIntSet visited = new TIntHashSet();
        TIntList t2 = new TIntArrayList();
        for (i = 0; i < ec.length; ++i) {
            if (!visited.contains(ec[i])) {
                t2.add(ec[i]);
                visited.add(ec[i]);
            }
        }
        t2.add(t2.get(0));
        int[] _t2 = t2.toArray();
        //System.out.printf("hamiltonian: %s\n", Arrays.toString(_t2));
        int[] expectedTSPApprox = new int[]{0, 1, 2, 3, 4, 0};
        assertTrue(Arrays.equals(expectedTSPApprox, _t2));
        
        int[] hTour = tsp.approxTSPTour(nVertexes, adjCostMap);
        assertTrue(Arrays.equals(expectedTSPApprox, hTour));
    }
    
    public void testApproxTSPTour1() {
        System.out.println("testApproxTSPTour1");
        
        /* Test of traveling salesman approximate tour. 
         * 
         * From Cormen, Leiserson, Rivest, and Stein chap 35 Approximation Algorithms, 35.2 Traveling-Salesman problem, Fig 35.2
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
        */
        
        
        int nVertexes = 8;
        int i, j;
        int n2 = nVertexes;
        
        double[][] x = new double[n2][];
        x[0] = new double[]{2, 5};
        x[1] = new double[]{2, 3};
        x[2] = new double[]{1, 2};
        x[3] = new double[]{4, 5};
        x[4] = new double[]{5, 4};
        x[5] = new double[]{4, 3};
        x[6] = new double[]{6, 3};
        x[7] = new double[]{3, 1};

        int[][] dist = new int[n2][];
        for (i = 0; i < n2; ++i) {
            dist[i] = new int[n2];
        }
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
                dist[i][j] = (int) Math.ceil(d);// rounding up to preserve triangle inequality
            }
        }
        
        TIntIntMap neighbors;
        TIntObjectMap<TIntIntMap> adjCostMap = new TIntObjectHashMap<TIntIntMap>();        
        for (i = 0; i < nVertexes; ++i) {
            neighbors = new TIntIntHashMap();
            for (j = 0; j < nVertexes; ++j) {
                neighbors.put(j, dist[i][j]);
            }
            adjCostMap.put(i, neighbors);
        }
        
        /*
        a, b, c,  h, d, e, f, g,
        0, 1, 2,  7, 3, 4, 5, 6

        In contrast, an optimal cost, non-crossing:
        a, b, c,  h, f, g, e, d
        0, 1, 2,  7, 5, 6, 4, 3
        */
        // 0  1 2 3 4 5 7 6 0
        // a  b c d e f h g a
        
        TSPChristofidesSerdyukov tsp = new TSPChristofidesSerdyukov();
        int[] hTour = tsp.approxTSPTour(nVertexes, adjCostMap);
        
        long optimalCost = TSPChristofidesSerdyukov.totalCost(
            new int[]{0, 1, 2,  7, 5, 6, 4, 3}, adjCostMap);
        
        long approxCost = TSPChristofidesSerdyukov.totalCost(
            hTour, adjCostMap);
                
        double ratio = (double)approxCost/(double)optimalCost;
        System.out.printf("approxCost=%d optCost=%d, approx/opt=%.2f, k-approx=%.1f\n",
            approxCost, optimalCost, ratio, 3./2.);
        assertTrue(ratio <= (1.5 + 0.2));
        
    }
    
    public void testApproxTSPTour2() {
        
        System.out.println("testApproxTSPTour2");
        
        int i, j, n2 = 6;
        int nVertexes = n2;
        int factor = 100;
        double[][] dist = new double[n2][n2];
        for (double[] row : dist) {
            Arrays.fill(row, 10000*factor);
        }
        for (i = 0; i < n2; ++i) {
            dist[i][i] = 0;
        }
        dist[5][0] = 10*factor;
        dist[1][5] = 12*factor;
        dist[4][1] = 2*factor;
        dist[2][4] = 4*factor;
        dist[3][2] = 6*factor;
        dist[0][3] = 8*factor;
        
        TIntIntMap neighbors;
        TIntObjectMap<TIntIntMap> adjCostMap = new TIntObjectHashMap<TIntIntMap>();        
        for (i = 0; i < nVertexes; ++i) {
            neighbors = new TIntIntHashMap();
            for (j = 0; j < nVertexes; ++j) {
                neighbors.put(j, (int)Math.ceil(dist[i][j]));
            }
            adjCostMap.put(i, neighbors);
        }
        
        double expectedCost = 42*factor; //8 + 6 + 4 + 2 + 12 + 10
        int[] optimalTour = new int[]{0, 3, 2, 4, 1, 5, 0};
        
        TSPChristofidesSerdyukov tsp = new TSPChristofidesSerdyukov();
        int[] hTour = tsp.approxTSPTour(nVertexes, adjCostMap);
        
        long optimalCost = TSPChristofidesSerdyukov.totalCost(optimalTour, adjCostMap);
        
        long approxCost = TSPChristofidesSerdyukov.totalCost(
            hTour, adjCostMap);
                
        double ratio = (double)approxCost/(double)optimalCost;
        System.out.printf("approxCost=%d optCost=%d, approx/opt=%.2f, k-approx=%.1f\n",
            approxCost, optimalCost, ratio, 3./2.);
        assertTrue(ratio <= 1.5);        
        
    }
    
    public void testApproxTSPTour3() throws Exception {
        
        System.out.println("testApproxTSPTour3");
        
        TIntList expectedTour = new TIntArrayList();
        
        double[][] dist = MatrixUtil.zeros(15, 15);
        
        double expectedOptimalCost = populateATT48_15_cities(expectedTour, dist);
        
        int i, j;
        int nVertexes = 15;
        
        TIntIntMap neighbors;
        TIntObjectMap<TIntIntMap> adjCostMap = new TIntObjectHashMap<TIntIntMap>();        
        for (i = 0; i < nVertexes; ++i) {
            neighbors = new TIntIntHashMap();
            for (j = 0; j < nVertexes; ++j) {
                neighbors.put(j, (int)Math.ceil(dist[i][j]));
            }
            adjCostMap.put(i, neighbors);
        }
        
        TSPChristofidesSerdyukov tsp = new TSPChristofidesSerdyukov();
        int[] hTour = tsp.approxTSPTour(nVertexes, adjCostMap);
        
        int[] optimalTour = expectedTour.toArray();
        long optimalCost = TSPChristofidesSerdyukov.totalCost(optimalTour, adjCostMap);
        
        long approxCost = TSPChristofidesSerdyukov.totalCost(
            hTour, adjCostMap);
                
        double ratio = (double)approxCost/(double)optimalCost;
        System.out.printf("approxCost=%d optCost=%d, approx/opt=%.2f, k-approx=%.1f\n",
            approxCost, optimalCost, ratio, 3./2.);
        
        assertTrue(ratio <= 1.5);        
        
    }
    
    /**
     * 
     @param expectedBestTour
     @param dist
     @return The minimal cost is 291
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
     @param filePath
     @param expectedBestTour
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
                //System.out.printf("items.length=%d", items.length);
                if (items[0] == null || items[0].equals("")) {
                    items = Arrays.copyOfRange(items, 1, 
                        items.length);
                }
                //System.out.printf("->%d\n", items.length);
                
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
