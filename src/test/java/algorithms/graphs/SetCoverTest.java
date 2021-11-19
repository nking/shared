package algorithms.graphs;

import gnu.trove.list.TIntList;
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
        from Section 35.3 of Cormen at al. "Introduction to Comptuer Algorithms"
        */
        List<TIntSet> sets = new ArrayList<TIntSet>();
        int n = 6, i;
        for (i = 0; i < n; ++i) {
            sets.add(new TIntHashSet());
        }
        sets.get(0).addAll(new int[]{0, 1, 2, 3, 4, 5});//S1
        sets.get(1).addAll(new int[]{4, 5, 7, 8});      //S2
        sets.get(2).addAll(new int[]{0, 3, 6, 9});      //S3
        sets.get(3).addAll(new int[]{1, 4, 6, 7, 10});  //S4
        sets.get(4).addAll(new int[]{2, 5, 8, 11});     //S5
        sets.get(5).addAll(new int[]{9, 10});           //S6
        
        SetCover sc = new SetCover();
    
        TIntList setCoverIndexes = sc.approxLgN(sets);
        
        // optimal minimum set cover: s2, s3, s4
        TIntSet expected0 = new TIntHashSet(new int[]{0, 3, 4, 5});
        TIntSet expected1 = new TIntHashSet(new int[]{0, 3, 4, 2});
        
        assertEquals(expected0.size(), setCoverIndexes.size());
        assertTrue(setCoverIndexes.containsAll(expected0) || setCoverIndexes.containsAll(expected1));
    }

}
