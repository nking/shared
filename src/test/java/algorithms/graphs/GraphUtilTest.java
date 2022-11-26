package algorithms.graphs;

import algorithms.matrix.MatrixUtil;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

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
}
