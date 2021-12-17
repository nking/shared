package algorithms.graphs;

import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class AllMaximalCliquesBronKerboschTest extends TestCase {
    
    public AllMaximalCliquesBronKerboschTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        // from https://en.m.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
        TIntObjectMap<TIntSet> g = new TIntObjectHashMap<TIntSet>();
        int n = 6, i, j;
        for (i = 0; i < n; ++i) {
            g.put(i, new TIntHashSet());
        }
        g.get(1-1).add(5-1);
        g.get(1-1).add(2-1);
        g.get(2-1).add(5-1);
        g.get(2-1).add(3-1);
        g.get(2-1).add(1-1);
        g.get(3-1).add(4-1);
        g.get(3-1).add(2-1);
        g.get(4-1).add(6-1);
        g.get(4-1).add(5-1);
        g.get(4-1).add(3-1);
        g.get(5-1).add(4-1);
        g.get(5-1).add(2-1);
        g.get(5-1).add(1-1);
        g.get(6-1).add(4-1);
        
        List<TIntSet> expected = new ArrayList<TIntSet>();
        expected.add(new TIntHashSet(new int[]{2-1, 3-1}));
        expected.add(new TIntHashSet(new int[]{1-1, 2-1, 5-1}));
        expected.add(new TIntHashSet(new int[]{3-1, 4-1}));
        expected.add(new TIntHashSet(new int[]{4-1, 5-1}));
        expected.add(new TIntHashSet(new int[]{4-1, 6-1}));
        
        AllMaximalCliquesBronKerbosch mc = new AllMaximalCliquesBronKerbosch(g);
        List<TIntSet> results = mc.bronKerBosch();
        
        for (j = 0; j < results.size(); ++j) {
            System.out.println(Arrays.toString(results.get(j).toArray()));
        }
        
        assertEquals(expected.size(), results.size());
        for (i = 0; i < expected.size(); ++i) {
            for (j = 0; j < results.size(); ++j) {
                if (expected.get(i).size() == results.get(j).size()) {
                    if (expected.get(i).containsAll(results.get(j))) {
                        results.remove(j);
                        break;
                    }
                }
            }
        }
        assertTrue(results.isEmpty());

    }
}
