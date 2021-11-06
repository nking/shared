package algorithms.optimization;

import algorithms.bipartite.Graph;
import algorithms.bipartite.MinCostUnbalancedAssignment;
import algorithms.graphs.HierholzersEulerCircuit;
import algorithms.matrix.MatrixUtil;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;
import thirdparty.HungarianAlgorithm;

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

        // multiply distances by 100 and round them to integers
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
                dist[i][j] = (int) Math.round(100.*d);
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
        assertTrue(ratio <= 1.5);
        
        
        int z = 1;
    }
    
}
