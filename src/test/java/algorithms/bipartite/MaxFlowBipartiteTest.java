package algorithms.bipartite;

import java.util.*;

import junit.framework.TestCase;

public class MaxFlowBipartiteTest extends TestCase {

    public void test0() {
        int nV = 8;
        int[][] edges = new int[][]{
                {1,5}, {2,7}, {3,5}, {3,6}, {3,8}, {4,7}
        };
        Set<Integer> outVertexCover = new HashSet<>();
        Map<Integer, Integer> matched = MaxFlowBipartite.pairsAndMinimumVertexCover(edges, nV, outVertexCover);
        assertEquals(3, matched.size());
        assertEquals(5, matched.get(1).intValue());
        assertEquals(7, matched.get(2).intValue());
        assertTrue(matched.get(3).intValue()==8 || matched.get(3).intValue()==6);

        assertEquals(3, outVertexCover.size());
        assertTrue(outVertexCover.contains(1));
        assertTrue(outVertexCover.contains(3));
        assertTrue(outVertexCover.contains(7) || outVertexCover.contains(2));
    }
}
