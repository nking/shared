package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

import java.util.*;

public class GraphUtilTest extends TestCase {

    public void testConvertGraph0() {

        /*
        0  1  1
        0  0  1
        1  0  0
        */
        SimpleLinkedListNode[] g = new SimpleLinkedListNode[3];
        int i;
        for (i = 0; i < 3; ++i) {
            g[i] = new SimpleLinkedListNode();
        }
        g[0].insert(1); g[0].insert(2);
        g[1].insert(2);
        g[2].insert(0);

        TIntObjectMap<TIntSet> g2 = GraphUtil.convertGraph(g);
        assertEquals(3, g2.size());
        assertTrue(g2.containsKey(0));
        assertTrue(g2.containsKey(1));
        assertTrue(g2.containsKey(2));
        assertTrue(g2.get(0).contains(1)); assertTrue(g2.get(0).contains(2));
        assertTrue(g2.get(1).contains(2));
        assertTrue(g2.get(2).contains(0));
    }

    public void testConvertGraph1() {
        /*
        0  1  1
        0  0  1
        1  0  0
        */
        TIntObjectMap<TIntSet> g = new TIntObjectHashMap<>();
        g.put(0, new TIntHashSet());
        g.put(1, new TIntHashSet());
        g.put(2, new TIntHashSet());
        g.get(0).add(1); g.get(0).add(2);
        g.get(1).add(2);
        g.get(2).add(0);

        SimpleLinkedListNode[] g2 = GraphUtil.convertGraph(g);
        assertEquals(3, g2.length);
        SimpleLinkedListNode node0 = g2[0];
        SimpleLinkedListNode node1 = g2[1];
        SimpleLinkedListNode node2 = g2[2];
        TIntSet set = new TIntHashSet();
        while (node0 != null && node0.getNumberOfKeys() > 0) {
            set.add(node0.getKey());
            node0 = node0.getNext();
        }
        assertEquals(2, set.size());
        assertTrue(set.contains(1));
        assertTrue(set.contains(2));
        set.clear();
        while (node1 != null && node1.getNumberOfKeys() > 0) {
            set.add(node1.getKey());
            node1 = node1.getNext();
        }
        assertEquals(1, set.size());
        assertTrue(set.contains(2));

        set.clear();
        while (node2 != null && node2.getNumberOfKeys() > 0) {
            set.add(node2.getKey());
            node2 = node2.getNext();
        }
        assertEquals(1, set.size());
        assertTrue(set.contains(0));
    }

    public void testReverseMapping() {

        Map<Integer, Set<Integer>> g, revG;
        int u;
        int nV;
        g = GraphHelper.getGraph4();
        nV = g.size();
        revG = GraphUtil.createReverseMapping(g);

        assertEquals(nV, revG.size());

        Iterator<Map.Entry<Integer, Set<Integer>>> iter = g.entrySet().iterator();
        Map.Entry<Integer, Set<Integer>> entry;
        while (iter.hasNext()) {
            entry = iter.next();
            u = entry.getKey();
            for (int v : entry.getValue()) {
                assertTrue(revG.containsKey(v));
                assertTrue(revG.get(v).contains(u));
            }
        }
    }

    public void testSubtractVertex() {

        Map<Integer, Set<Integer>> g, revG;
        int i;
        int nV;
        g = GraphHelper.getGraph4();
        nV = g.size();
        revG = GraphUtil.createReverseMapping(g);

        for (i = 0; i < nV; ++i) {
            GraphUtil.subtractVertex(i, g, revG);
        }
        assertEquals(0, g.size());
        assertEquals(0, revG.size());
    }

    public void testDegrees() {

        Map<Integer, Set<Integer>> g;
        int i;
        int nV;
        g = GraphHelper.getGraph4();
        nV = g.size();

        Set<Integer> v = new HashSet<Integer>();
        for (i = 0; i < nV; ++i) {
            v.add(i);
        }
        Map<Integer, Integer> degreeMap = GraphUtil.createDegreeMapForVertices(v, g);
        int d;
        for (i = 0; i < nV; ++i) {
            d = degreeMap.get(i);
            if (i == 6) {
                assertEquals(6, d);
            } else {
                assertEquals(3, d);
            }
        }

        int vMax = GraphUtil.findMaxDegreeVertex(degreeMap);
        assertEquals(6, vMax);
    }

    public void testCopyToOrderedAdjMap() {
        Map<Integer, Set<Integer>> g = new HashMap<Integer, Set<Integer>>();
        int n = 5;
        int i;
        for (i = 0; i < n; ++i) {
            g.put(i, new HashSet<Integer>());
            for (int j = 0; j < i; ++j) {
                g.get(i).add(j);
            }
        }

        TIntObjectMap<TIntList> adjList = GraphUtil.copyToOrderedAdjMap(g, false);
        assertEquals(n, adjList.size());
        for (i = 0; i < n; ++i) {
            assertEquals(i, adjList.get(i).size());
        }

        Set<PairInt> edges = GraphUtil.extractEdges(g);
        PairInt edge;
        for (i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                edge = new PairInt(i, j);
                assertTrue(edges.contains(edge));
            }
        }
    }

    public void textExtractEdgesUsingOrder() {
        Map<Integer, Set<Integer>> adjMap = GraphHelper.getGraph6();
        Set<PairInt> edges = GraphUtil.extractEdgesUsingLexicographicOrder(adjMap);
        assertEquals(10, edges.size());
    }
}
