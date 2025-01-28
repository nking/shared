package algorithms.bipartite;

import algorithms.maxFlow.FordFulkerson;

import java.util.*;

/**
 *  this class determines whether a perfect matching is possible,
 * given graph data, and if so, calculates the perfect matching.
 */
public class MaxFlowBipartite {

    /**
     * find a perfect matching of the graph given by edges, and if outVertexCover is not null,
     * calculate the minimum vertex cover.
     * @param edges array of edge endpoints.  e.g. edges[0] = {1, 4} for edge between node 1 and 4.
     * @param outVertexCover if not null, the minimum vertex cover is also calculated and returned
     *                       in this set.  WARNING: the minimum vertex cover is in a work in progress.
     *                       A quick greedy method is implemented temporarily.
     *                       NOTE that independent set can be constructed from the complement of the minimum
     *                       vertex cover.
     * @return the maximum size perfect matching of left and right nodes of the bipartite graph.
     * If the return is null, a perfect matching was not possible.
     */
    public static Map<Integer, Integer> pairsAndMinimumVertexCover(int[][] edges, int nVertices,
                                                      Set<Integer> outVertexCover) {
        // add a src and sink node
        // renumber the vertices to allow space for src and sink

        // build graph without src and sink first to use Hall's theorem
        Map<Integer, Integer> idxMap = new HashMap<>();
        Map<Integer, Integer> idxMapRev = new HashMap<>(); // to write to out

        //key = uIdx, value = map of vIdx, val=flow
        Map<Integer, Map<Integer, Integer>> graph = new HashMap<>();

        int u, v, uIdx, vIdx;

        Map<Integer, Integer> degreesMap = null;
        if (outVertexCover != null) {
            degreesMap = new HashMap<>();
        }

        for (int[] edge : edges) {
            u = edge[0];
            v = edge[1];
            if (idxMap.containsKey(u)) {
                uIdx = idxMap.get(u);
            } else {
                uIdx = idxMap.size() + 1;
                idxMap.put(u, uIdx);
            }
            if (idxMap.containsKey(v)) {
                vIdx = idxMap.get(v);
            } else {
                vIdx = idxMap.size() + 1;
                idxMap.put(v, vIdx);
            }

            graph.putIfAbsent(uIdx, new HashMap<>());
            graph.get(uIdx).put(vIdx, 1);

            idxMapRev.put(vIdx, v);
            idxMapRev.put(uIdx, u);

            if (degreesMap != null) {
                degreesMap.put(uIdx, degreesMap.getOrDefault(uIdx, 0) + 1);
                degreesMap.put(vIdx, degreesMap.getOrDefault(vIdx, 0) + 1);
            }
        }

        if (!isPossible(graph)) {
            return null;
        }

        // add src and sink nodes to use FordFulkerson
        int src = 0;
        int sink = nVertices + 1;
        graph.put(src, new HashMap<>());
        for (int[] edge : edges) {
            uIdx = idxMap.get(edge[0]);
            vIdx = idxMap.get(edge[1]);
            // add an edge from src to u
            graph.get(src).put(uIdx, 1);
            // add an edge from v to sink
            graph.putIfAbsent(vIdx, new HashMap<>());
            graph.get(vIdx).put(sink, 1);
        }

        FordFulkerson f = new FordFulkerson(graph, nVertices + 2, src, sink);
        long nMatchings = f.maxFlow();

        Map<Integer, Integer> matched = new HashMap<>();
        for (int[] edge : edges) {
            u = edge[0];
            if (matched.containsKey(u)) continue;
            uIdx = idxMap.get(u);
            for (Map.Entry<Integer, Integer> vF : f.getRemG().get(uIdx).entrySet()) {
                if (vF.getValue() == 0) {
                    vIdx = vF.getKey();
                    v = idxMapRev.get(vIdx);
                    matched.put(u, v);
                    break;
                }
            }
        }
        assert(matched.size() == nMatchings);
        if (outVertexCover == null) {
            return matched;
        }
        outVertexCover.clear();

        //we can remove half of the nodes.
        /*
        TODO: replace this with formal methods.  see VertexCover.java, though those are for more complex graphs.

        for each node in matched, calc the degree and store it as node, degree in an array.
        sort the array in desc order by degree.
        for each node in the list, remove all adjacent from the outVertexCover
        if there are not nMatchings left when done, this greedy approach is not a good solution.
        If the matching is small enough, exhaustive methods are possible.
         */
        int[][] degrees = new int[matched.size()*2][];
        int i = 0;
        for (Map.Entry<Integer, Integer> entry : matched.entrySet()) {
            u = entry.getKey();
            v = entry.getValue();
            uIdx = idxMap.get(u);
            vIdx = idxMap.get(v);
            degrees[i++] = new int[]{u, degreesMap.get(uIdx)};
            degrees[i++] = new int[]{v, degreesMap.get(vIdx)};
            outVertexCover.add(u);
            outVertexCover.add(v);
        }

        Arrays.sort(degrees, new Comparator<int[]>() {
            @Override
            public int compare(int[] o1, int[] o2) {
                // descending sort by degree
                return Integer.compare(o2[1], o1[1]);
            }
        });

        for (int[] uD : degrees) {
            u = uD[0];
            if (!outVertexCover.contains(u)) continue;
            uIdx = idxMap.get(u);
            // remove adjacent
            for (Map.Entry<Integer, Integer> entry : graph.get(uIdx).entrySet()) {
                vIdx = entry.getKey();
                if (vIdx == src || vIdx == sink) continue;
                v = idxMapRev.get(vIdx);
                outVertexCover.remove(v);
            }
        }

        return matched;
    }

    /**
     * using Hall's theorem, check whether a bipartite perfect matching is possible
     * @param g
     * @return
     */
    protected static boolean isPossible(Map<Integer, Map<Integer, Integer>> g) {
        return isPossible(g, g.keySet());
    }
    protected static boolean isPossible(Map<Integer, Map<Integer, Integer>> g, Set<Integer> leftNodes) {
        int nLeft = leftNodes.size();
        int nRight = 0;
        for (int u : leftNodes) {
            if (g.containsKey(u)) {
                nRight += g.get(u).size();
            }
        }
        return (nLeft <= nRight);
    }
}
