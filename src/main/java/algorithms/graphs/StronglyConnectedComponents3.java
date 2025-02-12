package algorithms.graphs;

import java.util.*;
import java.util.stream.Collectors;

public class StronglyConnectedComponents3 {

    /*
    a more readable version of strongly connected components.
    adapted from
    https://cp-algorithms.com/graph/strongly-connected-components.html

    This is Kosaraju’s algorithm
     */

    public StronglyConnectedComponents3() {}

    protected Set<Integer> visited = null;

    /**
     * given an DAG as an adjacency map, find strongly connected components and return
     * them in outputComponents, and return the condensed DAG of those components in
     * outputAdjMap.  Note that the condensed graph node numbers are the smallest vertexes
     * in their component.
     * r.t.c. is O(|V| + |E|).  s.c. is O(|V| + |E|)
     * @param adjMap map w/ key = vertex, value = set of adjacent vetices.  Note that the vertices can
     *               have negative numbers.
     * @param outputComponents
     * @param outputAdjMap
     */
    public void find(Map<Integer, Collection<Integer>> adjMap,
                     Map<Integer, Set<Integer>> outputComponents,
                     Map<Integer, Set<Integer>> outputAdjMap) {

        visited = new HashSet<>();
        outputComponents.clear();
        outputAdjMap.clear();

        List<Integer> order = new ArrayList<>();

        // first series of depth first searches
        for (int u : adjMap.keySet()) {
            if (!visited.contains(u)) {
                // order is by finish times for a node (after has visited all its descendants)
                dfs(u, adjMap, order);
            }
        }

        // create adjacency list of complementary graph G^T
        Map<Integer, Collection<Integer>> adjMapRev = new HashMap<>();
        for (Map.Entry<Integer, Collection<Integer>> entry : adjMap.entrySet()) {
            for (int v : entry.getValue()) {
                adjMapRev.putIfAbsent(v, new HashSet<>());
                adjMapRev.get(v).add(entry.getKey());
            }
        }

        visited.clear();

        // root vertex of the SCC that the vertex is in
        Map<Integer, Integer> roots = new HashMap<>();

        // reverse the DFS order list to use for traversal from top
        Collections.reverse(order);

        // second series of depth first searches
        for (int u : order) {
            if (visited.contains(u)) continue;

            List<Integer> component = new ArrayList<>();

            dfs(u, adjMapRev, component);

            // TODO: edit the DFS to return the min component to remove O(n) here
            int root = Collections.min(component); //O(n) at worst.
            //int root = component.get(0);

            outputComponents.put(root, component.stream().collect(Collectors.toSet()));

            for (int v : component) {
                roots.put(v, root);
            }
        }

        // add edges to condensation graph
        for (Map.Entry<Integer, Collection<Integer>>  entry : adjMap.entrySet()) {
            for (int v : entry.getValue()) {
                int pU = roots.get(entry.getKey());
                int pV = roots.get(v);
                if (pU != pV) {
                    outputAdjMap.putIfAbsent(pU, new HashSet<>());
                    outputAdjMap.get(pU).add(pV);
                }
            }
        }
    }

    protected void dfs(int u, Map<Integer, Collection<Integer>> adjMap, List<Integer> output) {
        visited.add(u);
        if (adjMap.containsKey(u)) {
            for (int v : adjMap.get(u)) {
                if (!visited.contains(v)) {
                    dfs(v, adjMap, output);
                }
            }
        }
        output.add(u);
    }

}
