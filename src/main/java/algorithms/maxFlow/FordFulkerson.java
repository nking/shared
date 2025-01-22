package algorithms.maxFlow;

import java.util.*;

public class FordFulkerson {

    private Map<Integer, Map<Integer, Integer>> remG;

    /**
     * reverse g tracks the edge flow, but in opposite direction.
     */
    protected Map<Integer, Map<Integer, Integer>> revG;

    protected final int src;

    protected final int sink;

    protected final int nVertices;
    protected final int nEdges;

    protected boolean finished = false;
    protected long maxFlow = -1;

    /**
     * constructor for the max flow algorithm.
     * @param g graph with nodes numbered from 0 to nVertices - 1.  the network
     *          should not have any anti-parallel edges, that is for an edge u->v
     *          there should be no edge v-> u
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

        // init maps
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

    protected int findAugPathEK(int[] visited, int[] prev) {
        Arrays.fill(visited, 0);
        Arrays.fill(prev, -1);

        // use BFS to find first path to reach dest sink

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

    protected int findAugPathDFS(int[] visited, int[] prev, Map<Integer, Map<Integer, Integer>> g) {
        throw new UnsupportedOperationException("not implemented");
    }

    protected int findAugPathDFSScaling(int[] visited, int[] prev, Map<Integer, Map<Integer, Integer>> g) {
        throw new UnsupportedOperationException("not implemented");
    }

    /**
     * the graph weighted by the remaining capacity for flow.  when there is no path without
     * a 0 weight in it, the graph flow is maximum.
     */
    public Map<Integer, Map<Integer, Integer>> getRemG() {
        return remG;
    }
}
