package algorithms.graphs;

import algorithms.util.PairInt;
import junit.framework.TestCase;

import java.util.*;

public class ColoringTest extends TestCase {

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testRLF() {

        rlf(GraphHelper.getGraph4(), 3);

        //rlf(getGraph3(), 3);

    }

    public void testDS() {

        ds(GraphHelper.getGraph4(), 3);

    }

    private void rlf(Map<Integer, Set<Integer>> g, int nExpected) {

        Map<Integer, Integer> cMap = new HashMap<Integer, Integer>();

        int n = Coloring.recursiveLargestFirst(g,cMap);

        assertEquals(nExpected, n);

        // assert that connected vertexes have different colors
        Iterator<Map.Entry<Integer, Set<Integer>>> iter = g.entrySet().iterator();
        int u;
        Map.Entry<Integer, Set<Integer>> entry;
        while (iter.hasNext()) {
            entry = iter.next();
            u = entry.getKey();
            for (int v : entry.getValue()) {
                assertTrue(cMap.get(u) != cMap.get(v));
            }
        }
    }

    private void ds(Map<Integer, Set<Integer>> g, int nExpected) {

        Map<Integer, Integer> cMap = new HashMap<Integer, Integer>();

        int n = Coloring.dSatur(g, cMap);

        assertEquals(nExpected, n);

        // assert that connected vertexes have different colors
        Iterator<Map.Entry<Integer, Set<Integer>>> iter = g.entrySet().iterator();
        int u;
        Map.Entry<Integer, Set<Integer>> entry;
        while (iter.hasNext()) {
            entry = iter.next();
            u = entry.getKey();
            for (int v : entry.getValue()) {
                assertTrue(cMap.get(u) != cMap.get(v));
            }
        }
    }

    public void testEdgeColoring() {

        Map<Integer, Set<Integer>> adjMap = GraphHelper.getGraph5();

        Map<PairInt, Integer> colorMap = new HashMap<PairInt, Integer>();

        int k = Coloring.edgeColoringMisraGies(adjMap, colorMap);

        assertEquals(4, k);

        //TODO: finish tests
    }
}
