package algorithms.msts;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import junit.framework.TestCase;

import java.util.Arrays;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Set;
import java.util.Map;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 *
 * @author nichole
 */
public class KruskalsMinimumSpanningTreeTest extends TestCase {
    
    public KruskalsMinimumSpanningTreeTest(String testName) {
        super(testName);
    }

    public void testSortWeightsNonDecreasing() {
        
        TObjectDoubleMap<PairInt> edgeWeights = new TObjectDoubleHashMap<PairInt>();
        edgeWeights.put(new PairInt(1, 2), 10);
        edgeWeights.put(new PairInt(2, 5), 1);
        edgeWeights.put(new PairInt(5, 7), 5);
        edgeWeights.put(new PairInt(2, 4), 2);
        edgeWeights.put(new PairInt(4, 5), 9);
        
        PairInt[] expected = new PairInt[5];
        expected[0] = new PairInt(2, 5);
        expected[1] = new PairInt(2, 4);
        expected[2] = new PairInt(5, 7);
        expected[3] = new PairInt(4, 5);
        expected[4] = new PairInt(1, 2);
        
        PairInt[] sortedW = KruskalsMinimumSpanningTree.sortWeightsNonDecreasing(
            edgeWeights);
        
        assertEquals(expected.length, sortedW.length);
        
        int i;
        PairInt e, s;
        for (i = 0; i < expected.length; ++i) {
            e = expected[i];
            s = sortedW[i];
            assertEquals(e, s);
        }
    }
    public void test0() {
        /*
        Following example in Cormen, Leiserson, Rivest, and Stein Chap 23, MST, Fig 23.4.
        
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
        int expectedWeightSum = 37;
        
        SimpleLinkedListNode[] graph = new SimpleLinkedListNode[n]; 
        for (idx = 0; idx < n; ++idx) {
            graph[idx] = new SimpleLinkedListNode();
        }
        /*
                 / [b] --- 8 ----/[c]-- 7  -- [d]\
               4    |          2     \         |    9
          [a]       11      [i]         4     14      [e]
               8    |     7     6              |   10
                 \ [h] / -- 1  --[g]-- 2 -- \ [f]/
        */
        TObjectDoubleMap<PairInt> edgeWeights = new TObjectDoubleHashMap<PairInt>();

        int[][] weights = new int[14][];
        graph[a].insert(b); graph[a].insert(h);
        edgeWeights.put(new PairInt(a, b), 4);
        edgeWeights.put(new PairInt(a, h), 8);
        weights[0] = new int[]{a, b, 4};
        weights[1] = new int[]{a, h, 8};

        graph[b].insert(c); graph[b].insert(h);
        edgeWeights.put(new PairInt(b, c), 8);
        edgeWeights.put(new PairInt(b, h), 11);
        weights[2] = new int[]{b, c, 8};
        weights[3] = new int[]{b, h, 11};
        
        graph[c].insert(d); graph[c].insert(f); graph[c].insert(i);
        edgeWeights.put(new PairInt(c, d), 7);
        edgeWeights.put(new PairInt(c, f), 4);
        edgeWeights.put(new PairInt(c, i), 2);
        weights[4] = new int[]{c, d, 7};
        weights[5] = new int[]{c, f, 4};
        weights[6] = new int[]{c, i, 2};

        graph[d].insert(e); graph[d].insert(f);
        edgeWeights.put(new PairInt(d, e), 9);
        edgeWeights.put(new PairInt(d, f), 14);
        weights[7] = new int[]{d, e, 9};
        weights[8] = new int[]{d, f, 14};

        graph[e].insert(f);
        edgeWeights.put(new PairInt(e, f), 10);
        weights[9] = new int[]{e, f, 10};

        graph[f].insert(g);
        edgeWeights.put(new PairInt(f, g), 2);
        weights[10] = new int[]{f, g, 2};
        
        graph[g].insert(h); graph[g].insert(i);
        edgeWeights.put(new PairInt(g, h), 1);
        edgeWeights.put(new PairInt(g, i), 6);
        weights[11] = new int[]{g, h, 1};
        weights[12] = new int[]{g, i, 6};

        graph[h].insert(i);
        edgeWeights.put(new PairInt(h, i), 7);
        weights[13] = new int[]{h, i, 7};
        
        TIntObjectMap<SimpleLinkedListNode> mst = KruskalsMinimumSpanningTree.
           mst(graph, edgeWeights);

        // has 1 and 7 extra (1 has a,h) (7 has d,e)
        List<Integer> mst2 = KruskalsMinimumSpanningTree.mst(weights);
        Set<Integer> v2 = new HashSet<>();
        for (int _idx : mst2) {
            v2.add(weights[_idx][0]);
            v2.add(weights[_idx][2]);
        }

        Set<Integer> expected = new HashSet<Integer>(
                Arrays.stream(new int[]{0, 1, 4, 5, 6, 7, 10, 11}).boxed().collect(Collectors.toList()));
        assertTrue(mst2.removeAll(expected));
        assertTrue(mst2.isEmpty());

        /*
        0=a to b
        1=a to h
        4,5,6 = c to d,f,i
        7=d to e
        11=h to g
        10=g to f
         */

        assertTrue(mst.containsKey(a));
        assertTrue(mst.get(a).contains(b));
        assertTrue(mst.get(a).contains(h));
        
        assertTrue(mst.containsKey(c));
        assertTrue(mst.get(c).contains(d));
        assertTrue(mst.get(c).contains(f));
        assertTrue(mst.get(c).contains(i));

        assertTrue(mst.containsKey(d));
        assertTrue(mst.get(d).contains(e));
        
        assertTrue(mst.containsKey(f));
        assertTrue(mst.get(f).contains(g));
        
        assertTrue(mst.containsKey(g));
        assertTrue(mst.get(g).contains(h));
    }

    public void test1() {
        /*
        Following example in Cormen, Leiserson, Rivest, and Stein Chap 23, MST, Fig 23.4.
        
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
        int expectedWeightSum = 37;

        Map<Integer, Map<Integer, Double>> graph = new HashMap<>();
        for (idx = 0; idx < n; ++idx) {
	    graph.put(idx, new HashMap<Integer, Double>());
        }	    
        
        /*
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

        double[] outSum = new double[1];

        Map<Integer, Map<Integer, Double>> mst = KruskalsMinimumSpanningTree
	       .mst(graph, outSum);

	//TODO: finish here
	assertTrue(Math.abs(outSum[0] - expectedWeightSum) < 1E-7);

	assertTrue(mst.containsKey(a));
        assertTrue(mst.get(a).containsKey(b));
        assertTrue(mst.get(a).containsKey(h));

        assertTrue(mst.containsKey(c));
        assertTrue(mst.get(c).containsKey(d));
        assertTrue(mst.get(c).containsKey(f));
        assertTrue(mst.get(c).containsKey(i));

        assertTrue(mst.containsKey(d));
        assertTrue(mst.get(d).containsKey(e));

        assertTrue(mst.containsKey(f));
        assertTrue(mst.get(f).containsKey(g));

        assertTrue(mst.containsKey(g));
        assertTrue(mst.get(g).containsKey(h));

    }
}
    
