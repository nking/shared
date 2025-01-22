package algorithms.maxFlow;

import junit.framework.TestCase;

import java.util.*;

public class FordFulkersonTest extends TestCase {

    public void _test0() {
        //CLRS 4th ed. fig 24.6,
        int src = 0;
        int sink = 5;
        int nV = 6;
        Map<Integer, Map<Integer, Integer>> g = new HashMap<>();
        g.put(0, new HashMap<>());
        g.get(0).put(1, 16);
        g.get(0).put(2, 13);
        g.put(1, new HashMap<>());
        g.get(1).put(3, 12);
        g.put(2, new HashMap<>());
        g.get(2).put(1, 4);
        g.get(2).put(4, 14);
        g.put(3, new HashMap<>());
        g.get(3).put(2, 9);
        g.get(3).put(sink, 20);
        g.put(4, new HashMap<>());
        g.get(4).put(3, 7);
        g.get(4).put(sink, 4);

        FordFulkerson f = new FordFulkerson(g, nV, src, sink);
        long maxFlow = f.maxFlow();
        assertEquals(23, maxFlow);
    }

    public void test1() {
        //competitive programmer's handbook,  chap 20.1
        int src = 1-1;
        int sink = 6-1;
        int nV = 6;
        Map<Integer, Map<Integer, Integer>> g = new HashMap<>();
        g.put(1-1, new HashMap<>());
        g.get(1-1).put(2-1, 5);
        g.get(1-1).put(4-1, 4);
        g.put(2-1, new HashMap<>());
        g.get(2-1).put(3-1, 6);
        g.put(3-1, new HashMap<>());
        g.get(3-1).put(6-1, 5);
        g.get(3-1).put(5-1, 8);
        g.put(4-1, new HashMap<>());
        g.get(4-1).put(2-1, 3);
        g.get(4-1).put(5-1, 1);
        g.put(5-1, new HashMap<>());
        g.get(5-1).put(6-1, 2);

        FordFulkerson f = new FordFulkerson(g, nV, src, sink);
        long maxFlow = f.maxFlow();
        assertEquals(7, maxFlow);
    }
}
