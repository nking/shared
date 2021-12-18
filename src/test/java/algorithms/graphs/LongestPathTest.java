package algorithms.graphs;

import algorithms.util.PairInt;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LongestPathTest extends TestCase {
    
    public LongestPathTest(String testName) {
        super(testName);
    }

    public void testSolve() {
        //System.out.println("solve");
        
        // example test is from https://www.techiedelight.com/find-cost-longest-path-dag/
        
        TObjectDoubleMap<PairInt> g = new TObjectDoubleHashMap<PairInt>();
        g.put(new PairInt(0, 6), 2);
        g.put(new PairInt(1, 2), -4);
        g.put(new PairInt(1, 4), 1);
        g.put(new PairInt(1, 6), 8);
        g.put(new PairInt(3, 0), 3);
        g.put(new PairInt(3, 4), 5);
        g.put(new PairInt(5, 1), 2);
        g.put(new PairInt(7, 0), 6);
        g.put(new PairInt(7, 1), -1);
        g.put(new PairInt(7, 3), 4);
        g.put(new PairInt(7, 5), -4);
        
        // topological sort: srcIdx=7, sIdx=0, ts=[7, 5, 3, 1, 2, 4, 0, 6]
        
        int srcIdx = 7;
        
        int[] expResult0 = new int[]{7, 3, 4};
        int[] expResult1 = new int[]{7, 3, 0, 6};
        
        int[] result = LongestPath.solve(g, srcIdx);
        
        assertTrue(Arrays.equals(expResult0, result));
        
        result = LongestPath.solve(g);
        
        assertTrue(Arrays.equals(expResult0, result));
        
    }
    
}
