package algorithms.graphs;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class GraphHelper {

    /**
     * cube with 2 maximal independent sets
     * @return graph w/ 2 colors
     */
    public static Map<Integer, Set<Integer>> getGraph0() {
        //https://en.wikipedia.org/wiki/Maximal_independent_set#Listing_all_maximal_independent_sets
        Map<Integer, Set<Integer>> adj = new HashMap<Integer, Set<Integer>>();

        int n = 8;
        int i;
        for (i = 0; i < n; ++i) {
            adj.put(i, new HashSet<Integer>());
        }

        adj.get(0).add(1); adj.get(0).add(2); adj.get(0).add(6);
        adj.get(1).add(0); adj.get(1).add(3); adj.get(1).add(7);
        adj.get(2).add(0); adj.get(2).add(3); adj.get(2).add(4);
        adj.get(3).add(1); adj.get(3).add(2); adj.get(3).add(5);
        adj.get(4).add(2); adj.get(4).add(5); adj.get(4).add(6);
        adj.get(5).add(3); adj.get(5).add(4); adj.get(5).add(7);
        adj.get(6).add(0); adj.get(6).add(4); adj.get(6).add(7);
        adj.get(7).add(1); adj.get(7).add(5); adj.get(7).add(6);

        return adj;
    }

    /**
     * star graph with 1 maximal independent set
     * @return graph w/ 2-colors
     */
    public static Map<Integer, Set<Integer>> getGraph1() {
        //https://en.wikipedia.org/wiki/Maximal_independent_set#Listing_all_maximal_independent_sets
        Map<Integer, Set<Integer>> adj = new HashMap<Integer, Set<Integer>>();

        int n = 9;
        int i;
        for (i = 0; i < n; ++i) {
            adj.put(i, new HashSet<Integer>());
        }
        for (i = 1; i < n; ++i) {
            adj.get(0).add(i);
        }

        return adj;
    }

    /**
     *
     * @return graph with 2 -colors
     */
    public static Map<Integer, Set<Integer>> getGraph2() {
        //https://en.wikipedia.org/wiki/Maximal_independent_set#Listing_all_maximal_independent_sets
        Map<Integer, Set<Integer>> adj = new HashMap<Integer, Set<Integer>>();

        int n = 6;
        int i;
        for (i = 0; i < n; ++i) {
            adj.put(i, new HashSet<Integer>());
        }

        adj.get(0).add(1);
        adj.get(1).add(0); adj.get(1).add(3); adj.get(1).add(4); adj.get(1).add(5);
        adj.get(2).add(3); adj.get(2).add(5);
        adj.get(3).add(1); adj.get(3).add(2); adj.get(3).add(5);
        adj.get(4).add(1);
        adj.get(5).add(1); adj.get(5).add(2); adj.get(5).add(3);

        return adj;
    }

    /**
     * Fig 1 from Brelaz 1979
     * @return graph with 3-colors
     */
    public static Map<Integer, Set<Integer>> getGraph3() {
        Map<Integer, Set<Integer>> adj = new HashMap<Integer, Set<Integer>>();

        int n = 10;
        int i;
        for (i = 0; i < n; ++i) {
            adj.put(i, new HashSet<Integer>());
        }

        adj.get(0).add(1);//1
        adj.get(1).add(0); adj.get(1).add(2); adj.get(1).add(3); adj.get(1).add(6); adj.get(1).add(8); adj.get(1).add(9);//6
        adj.get(2).add(1); adj.get(2).add(4); //2
        adj.get(3).add(1); adj.get(3).add(2); adj.get(3).add(4); adj.get(3).add(5); adj.get(3).add(8); adj.get(3).add(9); //6
        adj.get(4).add(3); adj.get(4).add(6); adj.get(4).add(7); //3
        adj.get(5).add(3); //1
        adj.get(6).add(1); adj.get(6).add(4); adj.get(6).add(7); adj.get(6).add(8);//4
        adj.get(7).add(4); adj.get(7).add(6); adj.get(7).add(8);//3
        adj.get(8).add(1); adj.get(8).add(3); adj.get(8).add(6); adj.get(8).add(7);//4
        adj.get(9).add(1); adj.get(9).add(3); //2

        return adj;
    }

    /**
     * wheel graph
     * @return graph with 3-colors
     */
    public static Map<Integer, Set<Integer>> getGraph4() {
        // https://en.m.wikipedia.org/wiki/Recursive_largest_first_algorithm
        Map<Integer, Set<Integer>> adj = new HashMap<Integer, Set<Integer>>();

        int n = 7;
        int i;
        for (i = 0; i < n; ++i) {
            adj.put(i, new HashSet<Integer>());
        }

        adj.get(0).add(1); adj.get(0).add(5); adj.get(0).add(6);
        adj.get(1).add(0); adj.get(1).add(2); adj.get(1).add(6);
        adj.get(2).add(1); adj.get(2).add(3); adj.get(2).add(6);
        adj.get(3).add(2); adj.get(3).add(4); adj.get(3).add(6);
        adj.get(4).add(3); adj.get(4).add(5); adj.get(4).add(6);
        adj.get(5).add(0); adj.get(5).add(4); adj.get(5).add(6);
        adj.get(6).add(0); adj.get(6).add(1); adj.get(6).add(2);
            adj.get(6).add(3); adj.get(6).add(4); adj.get(6).add(5);

        return adj;
    }
}
