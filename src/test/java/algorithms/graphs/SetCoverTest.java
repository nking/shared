package algorithms.graphs;

import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SetCoverTest extends TestCase {
    
    public SetCoverTest(String testName) {
        super(testName);
    }

    /**
     * Test of approx2 method, of class SetCover.
     
     0  1  2
     3  4  5
     6  7  8
     9  10 11
     */
    public void testApproxLgN() {
        /*
        from Gig 35.3 of Cormen at al. "Introduction to Comptuer Algorithms"
        */
        List<TIntSet> sets = new ArrayList<TIntSet>();
        int n = 6, i;
        for (i = 0; i < n; ++i) {
            sets.add(new TIntHashSet());
        }
        sets.get(0).addAll(new int[]{0, 1, 2, 3, 4, 5});
        sets.get(1).addAll(new int[]{4, 5, 7, 8});
        sets.get(2).addAll(new int[]{0, 3, 6, 9});
        sets.get(3).addAll(new int[]{1, 4, 6, 7, 10});
        sets.get(4).addAll(new int[]{2, 5, 8, 11});
        sets.get(5).addAll(new int[]{9, 10});
        
        SetCover sc = new SetCover();
    
        TIntSet setCover = sc.approxLgN(sets);
        
        // minimum set cover: s2, s3, s4
        TIntSet expected = new TIntHashSet();
        expected.addAll(sets.get(2));
        expected.addAll(sets.get(3));
        expected.addAll(sets.get(4));
        
        assertEquals(expected.size(), setCover.size());
        assertTrue(setCover.containsAll(expected));
    }

}
