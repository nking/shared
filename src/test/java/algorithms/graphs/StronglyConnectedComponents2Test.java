package algorithms.graphs;

import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class StronglyConnectedComponents2Test extends TestCase {
    
    public StronglyConnectedComponents2Test() {
    }
    
    @Override
    protected void setUp() throws Exception {
       
    }

    @Override
    protected void tearDown() throws Exception {
        
    }
    
    public void test0() {
        // test from Cormen, Leiserson, Rivest, and Stein Fig 22.9
     
        SimpleLinkedListNode[] gr = new SimpleLinkedListNode[8];
        for (int i = 0; i < gr.length; ++i) {
            gr[i] = new SimpleLinkedListNode();
        }
        int a = 0;
        int b = 1;
        int c = 2;
        int d = 3;
        int e = 4;
        int f = 5;
        int g = 6;
        int h = 7;
        
        gr[a].insert(b);
        gr[b].insert(c); gr[b].insert(e); gr[b].insert(f);
        gr[c].insert(d); gr[c].insert(g);
        gr[d].insert(c); gr[d].insert(h);
        gr[e].insert(a); gr[e].insert(f);
        gr[f].insert(g); 
        gr[g].insert(f); gr[g].insert(h);
        gr[h].insert(h);
        
        StronglyConnectedComponents2 scc = new StronglyConnectedComponents2();
        int[] components = scc.findStronglyConnectedComponents(gr);
        
        TIntObjectMap<TIntSet> cSets = new TIntObjectHashMap<TIntSet>();
        int cI, idx;
        for (idx = 0; idx < components.length; ++idx) {
            cI = components[idx];
            if (!cSets.containsKey(cI)) {
                cSets.put(cI, new TIntHashSet());
            }
            cSets.get(cI).add(idx);
        }
        
        assertEquals(4, cSets.size());
        //[0, 0, 1, 1, 0, 2, 2, 3]
        TIntSet c0 = cSets.get(0);
        assertEquals(3, c0.size());
        assertTrue(c0.contains(a));
        assertTrue(c0.contains(b));
        assertTrue(c0.contains(e));
        
        TIntSet c1 = cSets.get(1);
        assertEquals(2, c1.size());
        assertTrue(c1.contains(c));
        assertTrue(c1.contains(d));
       
        TIntSet c2 = cSets.get(2);
        assertEquals(2, c2.size());
        assertTrue(c2.contains(f));
        assertTrue(c2.contains(g));
        
        TIntSet c3 = cSets.get(3);
        assertEquals(1, c3.size());
        assertTrue(c3.contains(h));
    }
    
}
