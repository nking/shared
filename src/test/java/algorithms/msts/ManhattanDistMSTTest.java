package algorithms.msts;

import junit.framework.TestCase;

import java.util.*;

public class ManhattanDistMSTTest extends TestCase {

    public void test0() {

        /*int[][] points = new int[][]{
                {1,3}, {2,1}, {4,4}, {4,2}
        };*/
        int[][] points = new int[][]{
                {4,0}, {7,4}, {3,5}, {2,4}
        };
        Set<ManhattanDistMST.Node> edges = ManhattanDistMST.manhattanMSTEdges(points);
        assertEquals(4, edges.size());

        int[][] expectedEdges = new int[][]{
                {2,3},  {3,1}, {3,0}, {0,1}
        };

        for (ManhattanDistMST.Node edge: edges) {
            boolean found = false;
            for (int[] expEdge : expectedEdges) {
                if ((expEdge[0] == edge.weightIJ[1] && expEdge[1] == edge.weightIJ[2])
                || (expEdge[0] == edge.weightIJ[2] && expEdge[1] == edge.weightIJ[1])){
                    found = true;
                    break;
                }
            }
            assertTrue(found);
        }
        int tt = 2;

        /*
            2-3 2
            3-1 5
            3-0 6 = 13
        */
        List<long[]> expected = new ArrayList<>();
        expected.add(new long[] {2, 3, 2});
        expected.add(new long[] {3, 1, 5});
        expected.add(new long[] {3, 0, 6});

        List<long[]> mst = ManhattanDistMST.manhattanMST(points);
        assertEquals(3, mst.size());

        for (int i = 0; i < expected.size(); ++i) {
            long[] e = expected.get(i);
            long[] r = mst.get(i);
            assertTrue(Arrays.equals(e, r));
        }

    }
}
