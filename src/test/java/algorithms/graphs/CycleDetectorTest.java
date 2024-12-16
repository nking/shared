package algorithms.graphs;

import junit.framework.TestCase;

import java.util.*;

public class CycleDetectorTest extends TestCase {

    public void test0() {
        Map<Integer, Map<Integer, Integer>> graph = new HashMap<>();
        /*
        a->b ->c -> d->b
         */
        graph.putIfAbsent(0, new HashMap<>());
        graph.get(0).put(1, 1);
        graph.putIfAbsent(1, new HashMap<>());
        graph.get(1).put(2, 1);
        graph.putIfAbsent(2, new HashMap<>());
        graph.get(2).put(3, 1);
        graph.get(2).put(4, 1);
        graph.putIfAbsent(3, new HashMap<>());
        graph.get(3).put(2, 1);

        CycleDetector det = new CycleDetector();

        assertTrue(det.hasCycle(graph));

        graph.remove(3);
        assertFalse(det.hasCycle(graph));

        graph.putIfAbsent(3, new HashMap<>());
        graph.get(3).put(3, 1);
        assertTrue(det.hasCycle(graph));
    }
}
