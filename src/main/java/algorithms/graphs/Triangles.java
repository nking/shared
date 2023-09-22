package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.PolyInt;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class Triangles {

    /**
     * count the edges which form triangles in the given undirected graph.
     * runtime complexity is O(m^(3/2)) where m is the number of edges in the graph.
     *
     * Note that for very large undirected graphs, one can create a Count Triangles algorithm to run
     * on a massive parallel software computing architecture
     * and distributed file system like MapReduce for the same runtime complexity, O(m^(3/2)),
     * including the computation cost.
     * A single MapReduce job can be used to make multiway joins for edges using order constraints
     * and the relation E in a natural join on E(X, Y) &#8882; &#x22B3;  E(X, Z) &#8882; &#x22B3; E(Y, Z)
     <pre>
     reference: Leskovec, Rajaraman, and Ullman, "Mining of Massive Datasets" (a.k.a. MMDS), chap 10.7 and chap 2.
     </pre>
     * @param adjMap undirected graph
     * @return
     */
    public static int count(Map<Integer, Set<Integer>> adjMap) {

        //O(|V|*k) where k is the avg degree of a node
        Map<Integer, Integer> nodeDegreeMap = GraphUtil.createDegreeMapForVertices(adjMap.keySet(), adjMap);

        //O(|V|*k) where k is the avg degree of a node
        Set<PairInt> edges = GraphUtil.extractEdgesUsingLexicographicOrder(adjMap);

        int n = adjMap.size();
        int m = edges.size();

        int mHH = (int)Math.sqrt(m);

        Set<PairInt> hh = new HashSet<>();
        Set<PolyInt> nonHHTriangles = new HashSet<PolyInt>();

        int u, v;
        Set<Integer> uAdj;
        PairInt p2;
        PolyInt p3;
        //O(|V|*k) where k is the avg degree of a node
        for (PairInt p : edges) {
            u = p.getX();
            v = p.getY();
            if (nodeDegreeMap.get(u) >= mHH && nodeDegreeMap.get(v) >= mHH) {
                hh.add(p);
                continue;
            }
            if (nodeDegreeMap.get(u) >= mHH) {
                continue;
            }
            // count non-hh triangles
            uAdj = adjMap.get(u);
            if (uAdj == null || uAdj.isEmpty()) {
                continue;
            }
            for (int uA : uAdj) {
                //count the triangle {u, v, uA} if and only if the edge (uA, v) exists, and u ≺ uA.
                if (u < uA) {
                    // edge is undirected in the original graph so e(v, uA) == e(uA, v),
                    // but in the ordered edges set, it is present only as (uA, v) since uA < v.
                    // the triangle (u, v, uA) then has order: u < uA, uA < v, u < v
                    p2 = new PairInt(uA, v);
                    if (edges.contains(p2)) {
                        p3 = new PolyInt(new int[]{u, v, uA});
                        //if u==v and in this undirected adjMap, v adj is ame as u adj, this will already be present
                        nonHHTriangles.add(p3);
                    }
                }
            }
        }

        Set<PolyInt> hhTriangles = new HashSet<PolyInt>();

        // count triangles in hh
        // number of hh is <= mHH and degree of each hh is >= mHH
        // consider sets of 3 of these nodes: C(sqrt(m), 3) ~ O(m^(3/2))
        for (PairInt p : hh) {
            u = p.getX();
            v = p.getY();
            uAdj = adjMap.get(u);
            if (uAdj == null || uAdj.isEmpty()) {
                continue;
            }
            for (int uA : uAdj) {
                if (nodeDegreeMap.get(uA) >= mHH && u < uA) {
                    p2 = new PairInt(v, uA);
                    if (hh.contains(p2)) {
                        p3 = new PolyInt(new int[]{u, v, uA});
                        hhTriangles.add(p3);
                    }
                }
            }
        }
        int nTriangles = nonHHTriangles.size() + hhTriangles.size();
        return nTriangles;
    }
}
