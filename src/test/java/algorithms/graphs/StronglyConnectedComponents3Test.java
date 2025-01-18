package algorithms.graphs;

import junit.framework.TestCase;

import java.util.*;
import java.util.stream.Collectors;

/**
 *
 * @author nichole
 */
public class StronglyConnectedComponents3Test extends TestCase {

    public void test0() {
        // test from Cormen, Leiserson, Rivest, and Stein Fig 22.9
        // fig 20.9 in 4th ed
        int a = 0;
        int b = 1;
        int c = 2;
        int d = 3;
        int e = 4;
        int f = 5;
        int g = 6;
        int h = 7;
        Map<Integer, Collection<Integer>> gr = new HashMap<>();
        for (int i = 0; i < 8; ++i) {
            gr.put(i, new HashSet<>());
        }
        gr.get(a).add(b);
        gr.get(b).add(c); gr.get(b).add(e); gr.get(b).add(f);
        gr.get(c).add(d); gr.get(c).add(g);
        gr.get(d).add(c); gr.get(d).add(h);
        gr.get(e).add(a); gr.get(e).add(f);
        gr.get(f).add(g); 
        gr.get(g).add(f); gr.get(g).add(h);
        gr.get(h).add(h);

        List<Set<Integer>> expComp = new ArrayList<>();
        Map<Integer, Set<Integer>> expOutMap = new HashMap<>();

        expComp.add(Arrays.stream(new int[]{a, b, e}).boxed()
                .collect(Collectors.toSet()));
        expComp.add(Arrays.stream(new int[]{c, d}).boxed()
                .collect(Collectors.toSet()));
        expComp.add(Arrays.stream(new int[]{f, g}).boxed()
                .collect(Collectors.toSet()));
        expComp.add(Arrays.stream(new int[]{h}).boxed()
                .collect(Collectors.toSet()));

        Map<Integer, Set<Integer>> expAdjMap = new HashMap<>();
        expAdjMap.putIfAbsent(a, new HashSet<>());
        expAdjMap.get(a).add(c);  expAdjMap.get(a).add(f);
        expAdjMap.putIfAbsent(c, new HashSet<>());
        expAdjMap.get(c).add(f);  expAdjMap.get(c).add(h);
        expAdjMap.put(f, new HashSet<>());
        expAdjMap.get(f).add(h);
        // abe (gr0) -> cd (gr1)  repr 0->2
        // abe (gr0) -> fg (gr2)  repr:0->5
        // cd (gr1) -> fg (gr2)   repr:2->5
        // cd (gr1) -> h (gr3)    repr:2->7
        // fg (gr2) -> h (gr3)    repr:5->7

        Map<Integer, Set<Integer>> outputComp = new HashMap<>();
        Map<Integer, Set<Integer>> outAdjMap = new HashMap<>();
        StronglyConnectedComponents3 scc = new StronglyConnectedComponents3();
        scc.find(gr, outputComp, outAdjMap);

        assertEquals(expComp.size(), outputComp.size());

        assertEquals(expAdjMap.size(), outAdjMap.size());

        for (Set<Integer> comp : outputComp.values()) {
            Set<Integer> found = null;
            for (Set<Integer> exp : expComp) {
                if (exp.containsAll(comp) && exp.size() == comp.size()) {
                    found = exp;
                    break;
                }
            }
            assertNotNull(found);
            expComp.remove(found);
        }
        assertEquals(0, expComp.size());

        for (Map.Entry<Integer, Set<Integer>> entry : outAdjMap.entrySet()) {
            assertTrue(expAdjMap.containsKey(entry.getKey()));
            assertTrue(expAdjMap.get(entry.getKey()).containsAll(entry.getValue()));
        }

    }
    
}
