package algorithms.msts;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import static org.junit.Assert.assertTrue;
import junit.framework.TestCase;
import java.util.*;

/**
 *
 * @author nichole
 */
public class PrimsSimpleTest extends TestCase {

    public void test0() {
        
        /*
        Following example in Cormen, Leiserson, Rivest, and Stein Chap 24, MST, Fig 23.4.
        
                 / [b] --- 8 ----/[c]-- 7  -- [d]\
               4    |          2     \         |    9
          [a]       11      [i]         4     14      [e]
               8    |     7     6              |   10
                 \ [h] / -- 1  --[g]-- 2 -- \ [f]/
        */
        
        int a = 0; int b = 1; int c = 2; int d = 3; int e = 4; int f = 5;
        int g = 6; int h = 7; int i = 8;
        
        int n = 9;
        int idx;

        //int[][] mst(int[][] edges, double[] weights, int src, double[] outSum)

        Map<Integer, Map<Integer, Double>> graph = new HashMap<>();
        for (idx = 0; idx < n; ++idx) {
            graph.put(idx, new HashMap<Integer, Double>());
        }

        /*
        add links for both directions since traversal order of vertexes
           matters in prims
                 / [b] --- 8 ----/[c]-- 7  -- [d]\
               4    |          2     \         |    9
          [a]       11      [i]         4     14      [e]
               8    |     7     6              |   10
                 \ [h] / -- 1  --[g]-- 2 -- \ [f]/
        */
        graph.get(a).put(b, 4.);
        graph.get(a).put(h, 8.);

        graph.get(b).put(c, 8.);
        graph.get(b).put(h, 11.);
        
        graph.get(c).put(d, 7.);
        graph.get(c).put(f, 4.);
        graph.get(c).put(i, 2.);

        graph.get(d).put(e, 9.);
        graph.get(d).put(f, 14.);
        
        graph.get(e).put(f, 10.);
        graph.get(f).put(g, 2.);
        graph.get(g).put(h, 1.);
        graph.get(g).put(i, 6.);
        graph.get(h).put(i, 7.);

        /*
                 / [b] --- 8 ----/[c]-- 7  -- [d]\
               4    |          2     \         |    9
          [a]       11      [i]         4     14      [e]
               8    |     7     6              |   10
                 \ [h] / -- 1  --[g]-- 2 -- \ [f]/
        
        */
        
        int r = a;
        
        double[] outSum = new double[1]; 
        List<int[]> mst = PrimsSimple.mst(graph, r, outSum);

        assertTrue(Math.abs(37. - outSum[0]) < 1E-3);
       
        Set<int[]> expA = new HashSet<>();
        expA.add(new int[]{a,b});
        expA.add(new int[]{b,c});
        expA.add(new int[]{c,d});
        expA.add(new int[]{c,i});
        expA.add(new int[]{c,f});
        expA.add(new int[]{d,e});
        expA.add(new int[]{f,g});
        expA.add(new int[]{g,h});

        assertEquals(expA.size(), mst.size());

        for (int[] mstEdge : mst) {
            int[] del = null;
            for (int[] ea : expA) {
               if (Arrays.equals(mstEdge, ea)) {
                   del = ea;
                   break;
               }
            }
            assert(del != null);
            expA.remove(del);
        }
        assertTrue(expA.isEmpty());

    }

    public void test1() {

        /*
        Following example in Cormen, Leiserson, Rivest, and Stein Chap 24, MST, Fig 23.4.

                 / [b] --- 8 ----/[c]-- 7  -- [d]\
               4    |          2     \         |    9
          [a]       11      [i]         4     14      [e]
               8    |     7     6              |   10
                 \ [h] / -- 1  --[g]-- 2 -- \ [f]/
        */

        int a = 0; int b = 1; int c = 2; int d = 3; int e = 4; int f = 5;
        int g = 6; int h = 7; int i = 8;

        int n = 9;
        int idx;

        //int[][] mst(int[][] edges, double[] weights, int src, double[] outSum)


        /*
        add links for both directions since traversal order of vertexes
           matters in prims
                 / [b] --- 8 ----/[c]-- 7  -- [d]\
               4    |          2     \         |    9
          [a]       11      [i]         4     14      [e]
               8    |     7     6              |   10
                 \ [h] / -- 1  --[g]-- 2 -- \ [f]/
        */
        int[][] edges = new int[14][0];
        double[] weights = new double[14];

        edges[0] = new int[]{a,b}; weights[0] = 4;
        edges[1] = new int[]{a, h}; weights[1] = 8.;
        edges[2] = new int[]{b, c}; weights[2] = 8.;
        edges[3] = new int[]{b, h}; weights[3] = 11.;

        edges[4] = new int[]{c, d}; weights[4] = 7.;
        edges[5] = new int[]{c, f}; weights[5] = 4.;
        edges[6] = new int[]{c, i}; weights[6] = 2.;

        edges[7] = new int[]{d, e}; weights[7] = 9.;
        edges[8] = new int[]{d, f}; weights[8] = 14.;

        edges[9] = new int[]{e, f}; weights[9] = 10.;
        edges[10] = new int[]{f, g}; weights[10] = 2.;
        edges[11] = new int[]{g, h}; weights[11] = 1.;
        edges[12] = new int[]{g, i}; weights[12] = 6.;
        edges[13] = new int[]{h, i}; weights[13] = 7.;

        /*
                 / [b] --- 8 ----/[c]-- 7  -- [d]\
               4    |          2     \         |    9
          [a]       11      [i]         4     14      [e]
               8    |     7     6              |   10
                 \ [h] / -- 1  --[g]-- 2 -- \ [f]/

        */

        int r = a;

        double[] outSum = new double[1];
        List<int[]> mst = PrimsSimple.mst(edges, weights, r, outSum);

        assertTrue(Math.abs(37. - outSum[0]) < 1E-3);

        Set<int[]> expA = new HashSet<>();
        expA.add(new int[]{a,b});
        expA.add(new int[]{b,c});
        expA.add(new int[]{c,d});
        expA.add(new int[]{c,i});
        expA.add(new int[]{c,f});
        expA.add(new int[]{d,e});
        expA.add(new int[]{f,g});
        expA.add(new int[]{g,h});

        assertEquals(expA.size(), mst.size());

        for (int[] mstEdge : mst) {
            int[] del = null;
            for (int[] ea : expA) {
                if (Arrays.equals(mstEdge, ea)) {
                    del = ea;
                    break;
                }
            }
            assert(del != null);
            expA.remove(del);
        }
        assertTrue(expA.isEmpty());
    }

}
