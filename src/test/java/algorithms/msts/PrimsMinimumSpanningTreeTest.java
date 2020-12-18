package algorithms.msts;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PrimsMinimumSpanningTreeTest extends TestCase {
    
    public PrimsMinimumSpanningTreeTest(String testName) {
        super(testName);
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
        add links for both directions since traversal order of vertexes
           matters in prims
                 / [b] --- 8 ----/[c]-- 7  -- [d]\
               4    |          2     \         |    9
          [a]       11      [i]         4     14      [e]
               8    |     7     6              |   10
                 \ [h] / -- 1  --[g]-- 2 -- \ [f]/
        */
        TObjectIntMap<PairInt> edgeWeights = new TObjectIntHashMap<PairInt>();
        
        graph[a].insert(b); 
        edgeWeights.put(new PairInt(a, b), 4);
        //graph[a].insert(h); edgeWeights.put(new PairInt(a, h), 8); <--- removing to let heap order find b to c before h to a
        
        graph[b].insert(h); graph[b].insert(c); graph[b].insert(a);
        edgeWeights.put(new PairInt(b, c), 8);
        edgeWeights.put(new PairInt(b, h), 11);
        edgeWeights.put(new PairInt(b, a), 4);
        
        graph[c].insert(i); graph[c].insert(f); graph[c].insert(d); 
        edgeWeights.put(new PairInt(c, d), 7);
        edgeWeights.put(new PairInt(c, f), 4);
        edgeWeights.put(new PairInt(c, i), 2);
        graph[c].insert(b);
        edgeWeights.put(new PairInt(c, b), 8);
        
        graph[d].insert(f); graph[d].insert(e); 
        edgeWeights.put(new PairInt(d, e), 9);
        edgeWeights.put(new PairInt(d, f), 14);
        graph[d].insert(c);
        edgeWeights.put(new PairInt(d, c), 7);
        
        graph[e].insert(f);
        edgeWeights.put(new PairInt(e, f), 10);
        graph[e].insert(d);
        edgeWeights.put(new PairInt(e, d), 9);
        
        graph[f].insert(g);
        edgeWeights.put(new PairInt(f, g), 2);
        graph[f].insert(e); graph[f].insert(d); graph[f].insert(c);
        edgeWeights.put(new PairInt(f, e), 10);
        edgeWeights.put(new PairInt(f, d), 14);
        edgeWeights.put(new PairInt(f, c), 4);
        
        graph[g].insert(i); graph[g].insert(h); 
        edgeWeights.put(new PairInt(g, h), 1);
        edgeWeights.put(new PairInt(g, i), 6);
        graph[g].insert(f);
        edgeWeights.put(new PairInt(g, f), 2);
        
        graph[h].insert(i);
        edgeWeights.put(new PairInt(h, i), 7);
        graph[h].insert(g); graph[h].insert(b); graph[h].insert(a);
        edgeWeights.put(new PairInt(h, g), 1);
        edgeWeights.put(new PairInt(h, b), 11);
        edgeWeights.put(new PairInt(h, a), 8);
        
        graph[i].insert(h); graph[i].insert(g); graph[i].insert(c);
        edgeWeights.put(new PairInt(i, h), 7);
        edgeWeights.put(new PairInt(i, g), 6);
        edgeWeights.put(new PairInt(i, c), 2);
        
        /*
                 / [b] --- 8 ----/[c]-- 7  -- [d]\
               4    |          2     \         |    9
          [a]       11      [i]         4     14      [e]
               8    |     7     6              |   10
                 \ [h] / -- 1  --[g]-- 2 -- \ [f]/
        */
        
        
        int r = a;
        
        TIntObjectMap<SimpleLinkedListNode> mst = PrimsMinimumSpanningTree.
           mst(graph, edgeWeights, r);
        
        /*
        int a = 0; int b = 1; int c = 2; int d = 3; int e = 4; int f = 5;
        int g = 6; int h = 7; int i = 8;
        a to b
        b to c
        c to d
        c to i
        c to f
        d to e
        f to g
        g to h
        
        another mst is
        0 to 1
        5 to 2
        2 to 3
        3 to 4
        6 to 5
        7 to 6
        0 to 7
        2 to 8
        */
        assertTrue(mst.containsKey(a));
        assertTrue(mst.get(a).contains(b));
        
        assertTrue(mst.containsKey(b));
        assertTrue(mst.get(b).contains(c));
        
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
    
}
