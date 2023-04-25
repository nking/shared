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

    public void testEdgeColoring5() {

        Map<Integer, Set<Integer>> adjMap = GraphHelper.getGraph5();

        Map<PairInt, Integer> colorMap = new HashMap<PairInt, Integer>();
        int k = Coloring.edgeColoringMisraGies(adjMap, colorMap);

        assertEquals(4, k);
        assertEdgeColoring(adjMap, colorMap);

    }

    public void testEdgeColoring6() {

        Map<Integer, Set<Integer>> adjMap = GraphHelper.getGraph6();

        Map<PairInt, Integer> colorMap = new HashMap<PairInt, Integer>();
        int k = Coloring.edgeColoringMisraGies(adjMap, colorMap);

        assertEquals(5, k);
        assertEdgeColoring(adjMap, colorMap);
    }

    private void assertEdgeColoring(Map<Integer, Set<Integer>> adjMap, Map<PairInt, Integer> colorMap) {
        Iterator<Map.Entry<Integer, Set<Integer>>> iter = adjMap.entrySet().iterator();
        Set<Integer> chk = new HashSet<>();
        Map.Entry<Integer, Set<Integer>> entry;
        int u;
        PairInt edge;
        int c;
        while (iter.hasNext()) {
            entry = iter.next();
            u = entry.getKey();
            chk.clear();
            for (int v : entry.getValue()) {
                if (u <= v) {
                    edge = new PairInt(u, v);
                } else {
                    edge = new PairInt(v, u);
                }
                assertTrue(colorMap.containsKey(edge));
                c = colorMap.get(edge);
                assertTrue(chk.add(c));
            }
        }
    }
}
