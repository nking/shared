package algorithms.msts;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import junit.framework.TestCase;

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
        Following example in Cormen et al. Chap 24, MST, Fig 23.4.
        
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
        
        graph[a].insert(b); graph[a].insert(h);
        edgeWeights.put(new PairInt(a, b), 4);
        edgeWeights.put(new PairInt(a, h), 8);
        
        graph[b].insert(c); graph[b].insert(h);
        edgeWeights.put(new PairInt(b, c), 8);
        edgeWeights.put(new PairInt(b, h), 11);
        
        graph[c].insert(d); graph[c].insert(f); graph[c].insert(i);
        edgeWeights.put(new PairInt(c, d), 7);
        edgeWeights.put(new PairInt(c, f), 4);
        edgeWeights.put(new PairInt(c, i), 2);
        
        graph[d].insert(e); graph[d].insert(f);
        edgeWeights.put(new PairInt(d, e), 9);
        edgeWeights.put(new PairInt(d, f), 14);
        
        graph[e].insert(f);
        edgeWeights.put(new PairInt(e, f), 10);
        
        graph[f].insert(g);
        edgeWeights.put(new PairInt(f, g), 2);
        
        graph[g].insert(h); graph[g].insert(i);
        edgeWeights.put(new PairInt(g, h), 1);
        edgeWeights.put(new PairInt(g, i), 6);
        
        graph[h].insert(i);
        edgeWeights.put(new PairInt(h, i), 7);
        
        TIntObjectMap<SimpleLinkedListNode> mst = KruskalsMinimumSpanningTree.
           mst(graph, edgeWeights);
        
        /*
        a to b
        c to d, c to f, c to i
        f to g
        g to h
        */
        assertTrue(mst.containsKey(a));
        assertTrue(mst.get(a).contains(b));
        
        assertTrue(mst.containsKey(c));
        assertTrue(mst.get(c).contains(d));
        assertTrue(mst.get(c).contains(f));
        assertTrue(mst.get(c).contains(i));
        
        assertTrue(mst.containsKey(f));
        assertTrue(mst.get(f).contains(g));
        
        assertTrue(mst.containsKey(g));
        assertTrue(mst.get(g).contains(h));
    }
}
