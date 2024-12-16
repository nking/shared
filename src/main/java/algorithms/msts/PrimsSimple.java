package algorithms.msts;

import algorithms.graphs.GraphUtil;

import java.util.*;

public class PrimsSimple {

    /**
     given an undirected weighted graph adjMap, determine the minimum weight tree that 
     connects all vertices by a single edge between each pair and at the minimum total
     weight of edges among all spanning trees.  a minimum spanning tree is returned,
     though there may be more than one possible for the given graph with the same
     mst weight.

     runtime complexity using a heap
     is O(log(|V|) * (|V| + |E|))

     @param graph weighted undirected graph
     @src node to start the search from.  since the resulting minimum spanning tree contains
     all nodes, src can be any node. 
     */
    public static List<int[]> mst(Map<Integer, Map<Integer, Double>> graph,
        int src, double[] outSum) {

        GraphUtil util = new GraphUtil();
        int nV = util.countNodes(graph);

        reset(outSum);

        Map<Integer, Map<Integer, Double>> dirGraph = addBiDirection(graph);

        double[] dist = new double[nV];
        int[] prev = new int[nV];
        Arrays.fill(prev, -1);

        double sentinel = Double.POSITIVE_INFINITY;
        Arrays.fill(dist, sentinel);
        dist[src] = 0;

        TreeSet<double[]> q = new TreeSet<>((o1, o2)-> {
            int c = Double.compare(o1[1], o2[1]);
            if (c != 0) return c;
            return (int)o1[0] - (int)o2[0];
            });

        double[][] nodes = new double[nV][];

        for (int u = 0; u < nV; ++u) {
            nodes[u] = new double[]{u, dist[u]};
            q.add(nodes[u]);
        }

        int u, v;
        double w;
        double[] uw;
        while (!q.isEmpty()) {
            uw = q.pollFirst();
            u = (int)uw[0];
            nodes[u] = null;
            if (!dirGraph.containsKey(u)) continue;
            for (Map.Entry<Integer, Double> entry: dirGraph.get(u).entrySet()) {
                v = entry.getKey();
                if (nodes[v] == null || u==v) continue;
                w = entry.getValue();
                if (dist[v] > w) {
                    dist[v] = w;
                    prev[v] = u;
                    q.remove(nodes[v]);
                    nodes[v] = new double[]{v, w};
                    q.add(nodes[v]);
                }
            }
        }

        List<int[]> tree = new ArrayList<>();
        for (v = 0; v < nV; ++v) {
            if (prev[v] != -1) {
                tree.add(new int[]{prev[v], v});
            }
            addTo(outSum, dist[v]);
        }

        assert(tree.size() == nV-1);

        return tree;
    }

    private static Map<Integer, Map<Integer, Double>> addBiDirection(Map<Integer, Map<Integer, Double>> graph) {
        Map<Integer, Map<Integer, Double>> out = new HashMap<>();
        for (int u : graph.keySet()) {
            for (Map.Entry<Integer, Double> entry : graph.get(u).entrySet()) {
                int v = entry.getKey();
                double w = entry.getValue();
                out.putIfAbsent(u, new HashMap<>());
                out.putIfAbsent(v, new HashMap<>());
                out.get(u).put(v, w);
                out.get(v).put(u, w);
            }
        }
        return out;
    }

    /**
     given an undirected weighted graph adjMap, determine the minimum weight tree that 
     connects all vertices by a single edge between each pair and at the minimum total
     weight of edges among all spanning trees.  a minimum spanning tree is returned,
     though there may be more than one possible for the given graph with the same
     mst weight.
     @param adjMap weighted undirected graph
     @src node to start the search from.  since the resulting minimum spanning tree contains
     all nodes, src can be any node. 
     */
    public static List<int[]> mst(int[][] edges, double[] weights, int src, double[] outSum) {

        Map<Integer, Map<Integer, Double>> adjMap = new HashMap<>();
        for (int i = 0; i < edges.length; ++i) {
            adjMap.putIfAbsent(edges[i][0], new HashMap<Integer, Double>());
            adjMap.get(edges[i][0]).put(edges[i][1], weights[i]);
        }

        return mst(adjMap, src, outSum);
    }

    protected static void reset(double[] outSum) {
        if (outSum != null && outSum.length > 0) {
            outSum[0] = 0;
        }
    }
    protected static void addTo(double[] outSum, double value) {
        if (outSum != null && outSum.length > 0) {
            outSum[0] += value;
        }
    }
}
