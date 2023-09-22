package algorithms.graphs;

import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

import java.util.*;

/**
 *
 * @author nichole
 */
public class TrianglesTest extends TestCase {

    public TrianglesTest(String testName) {
        super(testName);
    }
    
    public void test0() {

        /*
            0----1
            |  /    \
            2 -- 3 --4
             \   |
                 5

          n = 6
          m = 8*2, but 8 unique
          sqrt(m) = 2

          hh:[1,2,3]
          nTriangles = 2
        */
        Map<Integer, Set<Integer>> adjMap = new HashMap<Integer, Set<Integer>>();
        for (int i = 0; i < 6; ++i) {
            adjMap.put(i, new HashSet<Integer>());
        }
        adjMap.get(0).add(1);adjMap.get(0).add(2);
        adjMap.get(1).add(0); adjMap.get(1).add(2); adjMap.get(1).add(4);
        adjMap.get(2).add(0); adjMap.get(2).add(1); adjMap.get(2).add(3); adjMap.get(2).add(5);
        adjMap.get(3).add(2); adjMap.get(3).add(4); adjMap.get(3).add(5);
        adjMap.get(4).add(1); adjMap.get(4).add(3);
        adjMap.get(5).add(2); adjMap.get(5).add(3);

        int nT = Triangles.count(adjMap);
        assertEquals(2, nT);

    }

}
