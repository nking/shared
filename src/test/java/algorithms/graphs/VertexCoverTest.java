package algorithms.graphs;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class VertexCoverTest extends TestCase {
    
    public VertexCoverTest(String testName) {
        super(testName);
    }
    
    /**
     * test is from Fig. 35.1 of Cormen et al. "Introduction to Algorithms".
     */
    public void testApprox2() {
        
        /*
        b=1----c=2---/d=3
         |     |   /  | \
         |     | /    |   \
        a=0    e=4----f=5   \g=6  
        */
        TIntObjectMap<TIntSet> g = new TIntObjectHashMap<TIntSet>();
        int i;
        for (i = 0; i < 7;++i) {
            g.put(i, new TIntHashSet());
        }
        g.get(0).add(1);
        g.get(1).add(2);
        g.get(2).add(3);
        g.get(2).add(4);
        g.get(3).add(4);
        g.get(3).add(5);
        g.get(3).add(6);
        g.get(4).add(2);
        g.get(4).add(3);
        g.get(4).add(5);
        g.get(5).add(3);
        g.get(5).add(4);
        g.get(6).add(3);
        
        VertexCover vc = new VertexCover();
        TIntSet cover = vc.approx2(g);
        //System.out.printf("2-appox set cover=%s\n", Arrays.toString(cover.toArray()));
        
        TIntSet expected = new TIntHashSet(new int[]{1, 2, 3, 4, 5, 6});
        
        assertEquals(expected.size(), cover.size());
        assertTrue(expected.containsAll(cover.toArray()));
    }
    
}
