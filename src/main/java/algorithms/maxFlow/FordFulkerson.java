package algorithms.maxFlow;

import java.util.*;

/**
 * a class holding an implementation of the Ford-Fulkerson method with
 * ability to provide other implmentations.
 *
 * The Ford-Fulkerson method is a greedy solution to the maximum flow problem
 * of a flow network.
 *
 * A flow network is a graph whose edges have a capacity
 * and each edge receives a flow.  The sum of flow into a node must equal the sum
 * of flow out of a node unless the node is a source or sink.
 *
 * The max flow of a flow network is a state of the network that obtains the max
 * possible flow rate.
 */
public class FordFulkerson {

    /**
     * the graph with edges holding the remaining amount of flow that the edge can hold
     * (the reduced capacity due to flow having been push through the edge.  the flow that
     * was push through the edge is stored in revG).
     * The amount of flow the edge can hold is decreased from capacity as flow moves
     * out of it.
     * the map key is the start node of an edge, the value is a map with key = stop node of an
     * edge, value = amount of flow the edge can handle (reduced capacity).
     */
    private Map<Integer, Map<Integer, Integer>> remG;

    /**
     * reverse g tracks the edge flow, but in opposite direction for the edge.
     * the map key is the stop node of an edge, the value is a map with key = start node of an
     * edge, value = amount of flow through the edge.
     */
    protected Map<Integer, Map<Integer, Integer>> revG;

    /**
     * the the source node of the flow network
     */
    protected final int src;

    /**
     * the sink node of the flow network.
     */
    protected final int sink;

    /**
     * number of vertices in the flow networkd
     */
    protected final int nVertices;

    /**
     * number of edges in the flow network.
     */
    protected final int nEdges;

    /**
     * flag indicating whether the max flow has been solved.
     */
    protected boolean finished = false;

    /**
     * the amount of flow pushed through the system, maximized.
     */
    protected long maxFlow = -1;

    /**
     * constructor for the max flow algorithm.
     * @param g graph with nodes numbered from 0 to nVertices - 1.  the network
     *          should not have any anti-parallel edges, that is for an edge u to v
     *          there should be no edge v to u
     *          for a pair(u,v) of vertices.
     * @param nVertices the number of vertices in the graph
     * @param src
     * @param sink
     */
    public FordFulkerson(Map<Integer, Map<Integer, Integer>> g, int nVertices, int src,
                       int sink) {
        this.src = src;
        this.sink = sink;
        this.nVertices = nVertices;

        this.remG = new HashMap<>();
        this.revG = new HashMap<>();
        int nEdges = 0;
        for (Map.Entry<Integer, Map<Integer, Integer>> entry : g.entrySet()) {
            int u = entry.getKey();
            getRemG().putIfAbsent(u, new HashMap<>());
            for (Map.Entry<Integer, Integer> entry1 : entry.getValue().entrySet()) {
                int v = entry1.getKey();
                int cap = entry1.getValue();

                getRemG().get(u).put(v, cap);

                revG.putIfAbsent(v, new HashMap<>());
                revG.get(v).put(u, 0);
                ++nEdges;
            }
        }
        this.nEdges = nEdges;
    }

    /**
     * calculate the maximum flow that can be push through the system.
     * @return the maximum flow that can be push through the system.
     */
    public long maxFlow() {
        if (finished) return maxFlow;

        int[] visited = new int[nVertices];
        int[] prev = new int[nVertices];

        long flowSum = 0;

        long nIterMax = nVertices * nEdges;
        int nIter = 0;
        while (nIter < nIterMax) {
            // find an augmenting flow
            int flow = findAugPathEK(visited, prev);

            //printDebug(String.format("\nfound flow=%d, prev=%s", flow,
            //        Arrays.toString(prev)));

            if (flow == 0) {
                break;
            }

            flowSum += flow;

            // augment the path in prev from sink to src
            int v = sink;
            int u = prev[v];
            while (v != src) {

                // subtract flow from u,v
                getRemG().get(u).put(v, getRemG().get(u).get(v) - flow);

                //and store in v,u
                revG.get(v).put(u, revG.get(v).get(u) + flow);

                v = u;
                u = prev[v];
            }

            //printDebug("augmented path:");

            ++nIter;
        };
        this.finished = true;
        this.maxFlow = flowSum;
        return maxFlow;
    }

    private void printDebug(String label) {
        System.out.println(label);
        for (Map.Entry<Integer, Map<Integer, Integer>> entry: getRemG().entrySet()) {
            int u = entry.getKey();
            for (Map.Entry<Integer, Integer> entryVW : entry.getValue().entrySet()) {
                int v = entryVW.getKey();
                int rem = entryVW.getValue();
                int flow = revG.get(v).get(u);
                System.out.printf("(%2d, %2d) : %3d, %3d\n", u, v, rem, flow);
            }
        }
    }

    /**
     * find an augmenting path using the Edmonds-Karp algorithm.
     * @param visited an array that will be used internally to track visited nodes.  must be of
     *                length nVertices.
     * @param prev an array that will be used internally to track the parents of a node.
     *             must be of length nVertices.
     * @return the flow found for the path composed in prev array.  the value will be 0 when
     * no augmenting path is found in the graph remG.
     */
    protected int findAugPathEK(int[] visited, int[] prev) {
        Arrays.fill(visited, 0);
        Arrays.fill(prev, -1);

        // use BFS to find first path to reach destination sink.
        // for each edge, store the minimum remaining edge capacity
        // of the edge itself and the it's on path.

        //queue key = idx, minflow
        Queue<int[]> q = new ArrayDeque<>();
        q.offer(new int[]{src, Integer.MAX_VALUE});

        while (!q.isEmpty()) {
            int[] idxMF = q.poll();
            if (idxMF[0] == sink) {
                return idxMF[1];
            }

            if (!getRemG().containsKey(idxMF[0])) continue;

            for (Map.Entry<Integer, Integer> entry : getRemG().get(idxMF[0]).entrySet()) {
                int w = entry.getValue();
                if (w < 0) {
                    throw new IllegalStateException("error in alg.  remainder graph has neg weight");
                }
                if (w == 0) continue;

                int v = entry.getKey();
                if (visited[v] != 0) continue;

                prev[v] = idxMF[0];
                visited[v] = 1;

                q.offer(new int[]{v, Math.min(idxMF[1], w)});
            }
        }
        // didn't find an augmenting path
        return 0;
    }

    /*
    protected int findAugPathDFS(int[] visited, int[] prev, Map<Integer, Map<Integer, Integer>> g) {
        throw new UnsupportedOperationException("not implemented");
    }

    protected int findAugPathDFSScaling(int[] visited, int[] prev, Map<Integer, Map<Integer, Integer>> g) {
        throw new UnsupportedOperationException("not implemented");
    }*/

    /**
     * get the graph  with edges that hold remaining capacity for flow.  when there is no path without
     * a 0 weight in it, the graph flow is maximum.
     * @return the graph  with edges that hold remaining capacity for flow
     */
    public Map<Integer, Map<Integer, Integer>> getRemG() {
        return remG;
    }
}
