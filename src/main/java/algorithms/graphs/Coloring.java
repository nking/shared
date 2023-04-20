package algorithms.graphs;

import java.util.*;

/**
 The graph coloring problem (GCP) asks to color the n vertices of a given undirected graph G = (V, E), that is,
 to map each vertex in V to a color in a given set of available colors K, in such a way that adjacent vertices
 take different colors. The set of available colors is commonly mapped to the set of integers {1, . . . , |K|}.

 In the decision version of the problem, also called the (vertex) k-coloring problem,
 we are asked whether for some given k = |K| a k-coloring exists.

 In the optimization version, we are asked for the smallest number k, called the chromatic number χ(G),
 for which a k-coloring exists.

 In general, the k-coloring problem is NP-complete [9]. For the optimization version, bad results exist also
 in terms of approximation, for example, a ratios of n1−ε cannot be achieved in polynomial time unless ZPP=NP [8]
<pre>
 Reference:
paper "Efficiency issues in the RLF heuristic for graph coloring"
by Marco Chiarandini, Giulia Galbiati, and Stefano Gualandi,
 MIC 2011: The IX Metaheuristics International Conference S1-47–1
 </pre>
 * see also the interval partitioning methods in algorithms.scheduling.Misc.java
 */
public class Coloring {

    /**
     * A polynomial time constructive algorithm to solve heuristically the graph coloring problem
     * It doesn't exhibit guaranteed approximation ratios but is very fast and produces good solutions in practice.
     * These features make them very appealing in practical applications.
     * Recursive Largest First (RLF) algorithm has a strategy to sequentially color stable sets,
     * that is, it sequentially builds sets of vertices that can take the same color.
     * Computational studies on these algorithms show that RLF clearly outperforms the other two in terms of
     * quality on a wide range of graph classes [5]. However, RLF comes with a higher computational cost,
     * having a O(n^3) worst-case complexity,
     <pre>
     Reference:
     paper "Efficiency issues in the RLF heuristic for graph coloring"
     by Marco Chiarandini, Giulia Galbiati, and Stefano Gualandi,
     MIC 2011: The IX Metaheuristics International Conference S1-47–1
     </pre>
     * @param adjMap the graph represented as an adjacency map
     * @param colorMap the output map having key=vertex, value= color where the range of color is [0, k)
     * @return the number of colors k.
     */
    public static int recursiveLargestFirst(Map<Integer, Set<Integer>> adjMap, Map<Integer, Integer> colorMap) {

        if (colorMap == null) {
            throw new IllegalArgumentException("colorMap must be instantiated by caller");
        }
        if (adjMap == null) {
            return 0;
        }

        adjMap = GraphUtil.copy(adjMap);

        Map<Integer, Set<Integer>> revAdjMap = GraphUtil.createReverseMapping(adjMap);

        int k = 0;
        while (!adjMap.isEmpty()) {
            findStableSet(adjMap, revAdjMap, k, colorMap);
            ++k;
        }
        return k;

        //TODO: consider generating graphs in tests using these models:
        //Erdo ̈s-Re ́nyi model (see, e.g., [7]).
        //Culberson’s generator for weight based models to produce graphs specifying the parameters α and γ to be 0 and 1 respectively.

        /*
        NOTE:  A stable set is also known as an independent set, coclique or anticliqu.
        It is a set of vertices in a graph, no two of which are adjacent.

        from Chiarandini et al. 2011:

        Let G[X] be the subgraph of G induced by the set of vertices X, i.e.,
            G[X] = (X, {uv ∈ E(G)|u, v ∈ X}).
        Denote by δX(v) the set of vertices adjacent to v in G[X ∪ v].
        Let dX(v) = |δX(v)| be the degree of v induced by X.
        Further, let P be the set of uncolored vertices
        and U the set of vertices that cannot be selected for becoming part of the current stable set.

        Initially, P is set equal to V, the set of vertices in the reduced graph passed to the
        FINDSTABLESET procedure, U is empty and the degree induced by U for all vertices is equal to 0 (Line 2).

        The procedure successively selects a vertex v to be colored (Lines 4 and 5),
        moves its neighbors δP(v) from P to U (Line 8), reduces the input graph G (Lines 6 and 9),
        and updates the degree induced by U of every vertex adjacent to a vertex that has moved from P to U .

        Algorithm 1 is an iterative extraction of stable sets from the (reduced) graph G = (V,E).
        The core ideas are given in Algorithm 2 that details the FINDSTABLESET procedure.

        Algorithm 1 RECURSIVE LARGEST FIRST(G)
           In G = (V, E) : input graph
           Out k : upper bound on χ(G)
           Out c : a coloring c:V 􏰄→K of G

        1. k←0
        2. while |V|>0 do
        3.    k ← k + 1 // increment the color
        4.    FINDSTABLESET(V, E, k) // G=(V,E) is reduced
        5. end while
        6. return k

        Algorithm 2 - procedure FINDSTABLESET(G, k)
            In G = (V, E) : input graph (in output G is the reduced graph)
            In k : color for current stable set
            Var P : set of potential vertices for the stable set
            Var U : set of vertices not in the current stable set

        1:  P ← V, U ← ∅
        2:  for all v∈P do d_U(v)←0
        3:  while |P|>0 do
        4:      v ← argmax_{w∈P} d_U(w)  // vertex w/ max degree induced by U
        5:      c(v) ← k // vertex v takes color k
        6:      V←V\{v}
        7:      for all w∈δ(v) do
        8:          if w∈P then P←P\{w}, U←U∪{w} end if // move w from P to U
        9:          E ← E\{v, w} // remove v from δ(w)
        10:         for all u ∈ δ(w) do
        11:             if u∈P then d_U(u)←d_U(u) + 1
        12:         end for
        13:     end for
        14: end while

         */
    }

    /**
     * A polynomial time constructive algorithm to solve heuristically the graph coloring problem
     * It doesn't exhibit guaranteed approximation ratios but is very fast and produces good solutions in practice.
     * These features make them very appealing in practical applications.
     * DSATUR [2] uses a dynamic order of the vertices,
     * instead of a static precomputed order. The idea is to sequentially color the vertices with the
     * smallest color, but the order is based on a more elaborated idea:
     * the next vertex is the one with the highest saturation degree, that is,
     * it has the highest degree induced by the colored vertices (ties are broken by the original degree,
     * in non increasing order).
     * The worst-case complexity is O(n^2).
     <pre>
     Reference:
     paper "Efficiency issues in the RLF heuristic for graph coloring"
     by Marco Chiarandini, Giulia Galbiati, and Stefano Gualandi,
     MIC 2011: The IX Metaheuristics International Conference S1-47–1
     </pre>
     NOTE: dSatur is a good algorithm to learn whether a graph is bipartite in O(n^2).
     * @param adjMap the graph represented as an adjacency map
     * @param colorMap the output map having key=vertex, value= color where the range of color is [0, k)
     * @return the number of colors k.
     */
    public static int dSatur(Map<Integer, Set<Integer>> adjMap, Map<Integer, Integer> colorMap) {

        //https://en.m.wikipedia.org/wiki/DSatur
        // Let the "degree of saturation" of a vertex be the number of different colours being used by its neighbors

        throw new UnsupportedOperationException("not yet implemented");

    }

    // NOTE:  A stable set is also known as an independent set, coclique or anticliqu.
    //        It is a set of vertices in a graph, no two of which are adjacent.
    private static void findStableSet(Map<Integer, Set<Integer>> adjMap,
                                      Map<Integer, Set<Integer>> revAdjMap, final int k, Map<Integer, Integer> c) {

        Set<Integer> p = new HashSet<Integer>(adjMap.keySet());
        Set<Integer> u = new HashSet<Integer>();

        // δX(v) = the set of vertices adjacent to v in G[X ∪ v].
        // dX(v) = |δX(v)| be the degree of v induced by X.

        //for all v∈P do d_U(v)←0
        Map<Integer, Integer> degreeUMap = createDegreeMapForVertices(u, p, adjMap);

        Set<Integer> vAdj, wAdj;

        int v;
        while (!p.isEmpty()) {
            //vertex with max degree induced by U
            //v ← argmax_{w∈P} d_U(w)
            v = GraphUtil.findMaxDegreeVertex(p, degreeUMap);
            /*if (v == -1) {
                //find the vertex in G w/ the largest degree
                v = GraphUtil.findMaxDegreeVertex(GraphUtil.createDegreeMapForVertices(p, adjMap));
            }*/
            c.put(v, k);

            vAdj = adjMap.get(v);

            p.remove(v);// p is only the uncolored vertexes, so remove v now
            GraphUtil.subtractVertex(v, adjMap, revAdjMap); // this removes v from V and E

            //δX(v) is the set of vertices adjacent to v in G[X ∪ v].
            // δX(v) = adjMap.get(v)
            // Let d_X(v) = |δX(v)| be the degree of v induced by X = adjMap.get(v).size()
            if (vAdj == null) {
                continue;
            }
            for (int w : vAdj) {
                //if w∈P then P←P\{w}, U←U∪{w}endif
                if (p.contains(w)) {
                    p.remove(w);
                    u.add(w);
                }

                // remove v from δ(w).  This was already done above with GraphUtil.subtractVertex(v...
                /*if (GraphUtil.removeEdge(v, w, adjMap)) {
                    GraphUtil.removeEdge(w, v, revAdjMap);
                }*/

                wAdj = adjMap.get(w);
                // update degreeUMap
                if (wAdj == null) {
                    continue;
                }
                for (int uu : wAdj) {
                    if (p.contains(uu)) {
                        assert(degreeUMap.containsKey(uu));
                        degreeUMap.put(uu, degreeUMap.get(uu) + 1);
                    }
                }
            }
        }
    }

    private static Map<Integer, Integer> createDegreeMapForVertices(Set<Integer> subsetOfG, Set<Integer> targetVertices,
                                                                   Map<Integer, Set<Integer>> g) {
        Map<Integer, Integer> degreeMap = new HashMap<Integer, Integer>();

        int nA;
        for (int v : targetVertices) {
            if (!subsetOfG.contains(v) || !g.containsKey(v)) {
                nA = 0;
            } else {
                nA = g.get(v).size();
            }
            degreeMap.put(v, nA);
        }
        return degreeMap;
    }
}
