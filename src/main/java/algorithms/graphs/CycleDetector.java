package algorithms.graphs;

import java.util.*;

public class CycleDetector {

    public <T> boolean hasCycle(Map<Integer, Map<Integer, T>> graph) {

        Map<Integer, Set<Integer>> g = new HashMap<>();
        for (int key : graph.keySet()) {
            for (int val : graph.get(key).keySet()) {
                g.putIfAbsent(key, new HashSet<Integer>());
                g.get(key).add(val);
            }
        }

        return hasCycle2(g);
    }

    private <T> boolean hasCycle(int u, int[] visited, int[] cycle, Map<Integer, Set<Integer>> graph) {
        visited[u] = 1;
        cycle[u] = 1;
        if (graph.containsKey(u)) {
            for (int v : graph.get(u)) {
                if (cycle[v] == 1) {
                    return true;
                }
                if (visited[v] != 0) continue;
                if (hasCycle(v, visited, Arrays.copyOf(cycle, cycle.length), graph)) {
                    return true;
                }
            }
        }
        cycle[u] = 2;
        return false;
    }

    public boolean hasCycle2(Map<Integer, Set<Integer>> graph) {
        // uses DFS to look for cycles in the graph
        GraphUtil util = new GraphUtil();
        int nV = util.countNodes2(graph);

        int[] visited = new int[nV];
        for (int i = 0; i < nV; ++i) {
            if (visited[i] == 0) {
                if (hasCycle(i, visited, new int[nV], graph)) {
                    return true;
                }
            }
        }
        return false;
    }
}
