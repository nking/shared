package algorithms.graphs;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.linked.TIntLinkedList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

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
     * A polynomial time constructive algorithm to solve heuristically the graph coloring problem for vertexes.
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
     * A polynomial time constructive algorithm to solve heuristically the graph coloring problem for vertexes.
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

     Other references used in implementing the algorithm:
         https://en.m.wikipedia.org/wiki/DSatur
         and the DSatur implementation by R. M. R. Lewis, School of Mathematics, Cardiff University, Wales
         http://rhydlewis.eu/gcol/
     </pre>
     NOTE: dSatur is a good algorithm to learn whether a graph is bipartite in O(n^2).
     * @param adjMap the graph represented as an adjacency map
     * @param colorMap the output map having key=vertex, value= color where the range of color is [0, k)
     * @return the number of colors k.
     */
    public static int dSatur(final Map<Integer, Set<Integer>> adjMap, final Map<Integer, Integer> colorMap) {

        if (colorMap == null) {
            throw new IllegalArgumentException("colorMap cannot be null");
        }

        int n = adjMap.size();

        // The edge chromatic number is never more than △(G)+1 where is the maximum degree of the graph.
        Map<Integer, Integer> degreeMap = GraphUtil.createDegreeMapForVertices(adjMap.keySet(), adjMap);
        int maxDegree = degreeMap.get(GraphUtil.findMaxDegreeVertex(degreeMap));

        Map<Integer, Set<Integer>> adjColors = new HashMap<>();

        SaturatedItem[] sItems = new SaturatedItem[n];

        // current largest color used
        int k = 0;

        PriorityQueue<SaturatedItem> pQ = new PriorityQueue<>();

        Iterator<Map.Entry<Integer, Set<Integer>>> iter = adjMap.entrySet().iterator();
        Map.Entry<Integer, Set<Integer>> entry;
        int i;
        while (iter.hasNext()) {
            entry = iter.next();
            i = entry.getKey();

            sItems[i] = new SaturatedItem(0, degreeMap.get(i), i);
            pQ.offer(sItems[i]);

            adjColors.put(i, new HashSet<>());
        }

        boolean[] used = new boolean[maxDegree + 1];

        SaturatedItem s;
        int u;
        int c;
        int d;
        Set<Integer> uAdj;
        while (!pQ.isEmpty()) {
            s = pQ.poll();
            u = s.vertex;

            // get smallest feasible color, not used by adjacent vertexes
            uAdj = adjMap.get(u);
            c = 0;
            if (uAdj != null) {
                Arrays.fill(used, false);
                for (int v : uAdj) {
                    if (colorMap.containsKey(v)) {
                        used[colorMap.get(v)] = true;
                    }
                }
                for (i = 0; i < used.length; ++i) {
                    if (!used[i]) {
                        c = i;
                        break;
                    }
                }
            }

            colorMap.put(u, c);
            if (c > k) {
                k = c;
            }
            if (uAdj != null) {
                for (int v : uAdj) {
                    if (!colorMap.containsKey(v)) {
                        adjColors.get(v).add(c);
                        d = degreeMap.get(v) - 1;
                        degreeMap.put(v, d);
                        // for priority queue to change order need to remove and insert again.
                        //TODO: modify fibonacci heap to use a comparator including the decreaseKey method
                        sItems[v].deg = d;
                        pQ.remove(sItems[v]);
                        pQ.add(sItems[v]);
                    }
                }
            }
        }

        return k + 1;

        /*
        adapted from pseudocode in
        https://en.m.wikipedia.org/wiki/DSatur
        and the DSatur implementation by R. M. R. Lewis, School of Mathematics, Cardiff University, Wales
        http://rhydlewis.eu/gcol/
        http://rhydlewis.eu/resources/gCol.zip
        Lewis, R. (2021) A Guide to Graph Colouring: Algorithms and Applications (second ed.).
        Springer, isbn: 978-3-030-81053-5, doi: 10.1007/978-3-030-81054-2
        The copyright on the web page http://rhydlewis.eu/gcol/   follows:
        <pre>
        8.   Copyright Notice
         Redistribution and use in source and binary forms, with or without modification, of the code associated with
         this document are permitted provided that a citation is made to the publication given at the start of this
         document. Neither the name of the University nor the names of its contributors may be used to endorse or
         promote products derived from this software without specific prior written permission. This software is
         provided by the contributors “as is”' and any express or implied warranties, including, but not limited to,
         the implied warranties of merchantability and fitness for a particular purpose are disclaimed. In no event
         shall the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential
         damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or
         profits; or business interruption) however caused and on any theory of liability, whether in contract, strict
         liability, or tort (including negligence or otherwise) arising in any way out of the use of this software,
         even if advised of the possibility of such damage. This software is supplied without any support services.
        </pre>

         */
    }

    /**
     Find a nearly optimal edge coloring using the Misra & Gries 1992 algorithm.
     It uses at most one color more than the optimal and has a runtime complexity O(|V|*|E|).
     In contrast, the optimal edge coloring, in general, is NP-complete.
     <pre>
     references used in implementing the algorithm:
         https://en.m.wikipedia.org/wiki/Misra_%26_Gries_edge_coloring_algorithm
         https://en.wikipedia.org/wiki/Vizing%27s_theorem
     and
         Liang, Shen and Hu 1996
         "Parallel Algorithms for the Edge-Coloring and Edge-Coloring Update Problems"
         JOURNAL OF PARALLEL AND DISTRIBUTED COMPUTING 32, 66–73 (1996) ARTICLE NO. 0005
     and
         https://www.boost.org/doc/libs/1_67_0/boost/graph/edge_coloring.hpp
         which uses Boost Software License - Version 1.0 - August 17th, 2003
     </pre>
     * @param adjMap adjacency map for a simple undirected graph (no self-loops and no more than 1 edge between same vertices).
     *               It is expected that the invoker has already handled the bi-directional mappings in adjMap for edges (u,v) and (v,u).
     * @param colorMap output map of edge colorings where key=edge vertexes and value is the color assigned to the edge
     * @return the number of colors assigned to edges
     */
    public static int edgeColoringMisraGries(final Map<Integer, Set<Integer>> adjMap, final Map<PairInt, Integer> colorMap) {

        //For building the "fans" we need a fixed ordering of the edges, presumably by label = vertex number which are integers here.

        // key = vertex
        // value = adjacent vertexes
        TIntObjectMap<TIntList> adjMapOrdered = GraphUtil.copyToOrderedAdjMap(adjMap, true);

        int maxDegreeP1 = findMaxDegree(adjMapOrdered) + 1;

        Set<PairInt> edges = GraphUtil.extractEdgesUsingLexicographicOrder(adjMap);

        /*
        https://en.m.wikipedia.org/wiki/Misra_%26_Gries_edge_coloring_algorithm
        U = edges
        while U ≠ ∅ do
            Let (u, v) be any edge in U.
            Let F[1:k] be a maximal fan of u starting at F[1] = v.
            Let c be a color that is free on u and d be a color that is free on F[k].
            Invert the cdu path
            Let w ∈ V(G) be such that w ∈ F, F' = [F[1]...w] is a fan and d is free on w.
            Rotate F' and set c(u, w) = d.
            U:= U − {(u, v)}
        end while
         */
        int j;
        int maxColors = -1;
        PairInt edge;
        int u;
        int k;
        int colorC;
        int colorD;
        while (!edges.isEmpty()) {

            edge = edges.iterator().next();

            TIntList f = findMaximalFan(edge, adjMapOrdered, colorMap);
            int w;
            assert(!f.isEmpty());
            u = edge.getX();
            k = f.get(f.size() - 1);

            // for x.adj, find first missing colors in 0:max
            colorC = findFreeColor(u, adjMapOrdered, colorMap, maxDegreeP1);
            // for k.adj, find first missing colors in 0:max
            colorD = findFreeColor(k, adjMapOrdered, colorMap, maxDegreeP1);

            // invert cd path.  does nothing if no colors are set yet
            invertCDPath(u, colorC, colorD, adjMapOrdered, colorMap, maxDegreeP1);

            // find the first member of f which has no adjacent vertexes of color colorD
            w = -1;
            for (j = 0; j < f.size(); ++j) {
                w = f.get(j);
                if (isFree(w, colorD, adjMapOrdered, colorMap)) {
                    // w.adjacent colors are not colorD
                    break;
                }
            }
            if (w > -1) {
                // does nothing if f.get(0) == w
                rotateFan(u, f.get(0), w, adjMapOrdered, colorMap);

                // assign color d to edge (u, w)
                colorMap.put(new PairInt(u, w), colorD);
            }

            maxColors = Math.max(maxColors, Math.max(colorC, colorD));

            edges.remove(edge);
        }

        return maxColors + 1;
    }

    private static PairInt makeOrderedNode(int u, int v) {
        if (u <= v) {
            return new PairInt(u, v);
        }
        return new PairInt(v, u);
    }

    private static int findFreeColor(int vertex, TIntObjectMap<TIntList> adjMapOrdered, Map<PairInt, Integer> colorMap,
                                     int maxC) {
        int c = 0;
        while (c < maxC && !isFree(vertex, c, adjMapOrdered, colorMap)) {
            c++;
        }
        return c;
    }

    private static int findMaxDegree(TIntObjectMap<TIntList> adjMapOrdered) {
        if (adjMapOrdered == null) {
            throw new IllegalArgumentException("adjMapOrdered cannot be null");
        }
        int max = 0;
        int s;
        TIntObjectIterator<TIntList> iter = adjMapOrdered.iterator();
        while (iter.hasNext()) {
            iter.advance();
            s = iter.value().size();
            if (s > max) {
                max = s;
            }
        }
        return max;
    }

    private static void rotateFan(int vertex, int f0, int f1, TIntObjectMap<TIntList> adjMapOrdered,
      Map<PairInt, Integer> colorMap) {
        if (f0 == f1) {
            return;
        }
        // f0 and f1 are adjacent to vertex
        PairInt previous = makeOrderedNode(vertex, f0);
        PairInt current;
        int color;
        PairInt chk;
        for (int b = f0 + 1; b < f1; ++b) {
            current = new PairInt(vertex, b);

            chk = makeOrderedNode(vertex, b);
            assert(colorMap.containsKey(chk));
            color = colorMap.get(chk);

            // assign color to previous
            colorMap.put(previous, color);
        }
    }

    private static void invertCDPath(int vertex, int colorC, int colorD,
         TIntObjectMap<TIntList> adjMapOrdered, Map<PairInt, Integer> colorMap, int maxColors) {

        TIntList adj = adjMapOrdered.get(vertex);
        if (adj == null) {
            return;
        }
        int v;
        PairInt edge;
        PairInt chk;
        for (int i = 0; i < adj.size(); ++i) {
            v = adj.get(i);
            edge = new PairInt(vertex, v);
            chk = makeOrderedNode(vertex, v);

            if (colorMap.containsKey(chk) && colorMap.get(chk) == colorD) {
                invertCDPathRecursive(edge, colorD, colorC, adjMapOrdered, colorMap);
                return;
            }
        }
    }

    /**
     * @param edge
     * @param colorC
     * @param colorD
     * @param adjMapOrdered
     * @param colorMap
     */
    private static void invertCDPathRecursive(PairInt edge, int colorC, int colorD,
         TIntObjectMap<TIntList> adjMapOrdered, Map<PairInt, Integer> colorMap) {

        int x = edge.getX();
        int v = edge.getY();

        // assign color d to edge
        PairInt chk = makeOrderedNode(x, v);
        colorMap.put(chk, colorD);

        TIntList adj = adjMapOrdered.get(v);
        if (adj == null) {
            return;
        }
        int v2;
        PairInt edge2;
        PairInt chk2;
        for (int i = 0; i < adj.size(); ++i) {
            v2 = adj.get(i);
            edge2 = new PairInt(v, v2);
            chk2 = makeOrderedNode(v, v2);
            if (!chk.equals(chk2) && colorMap.containsKey(chk2) && colorMap.get(chk2) == colorD) {
                invertCDPathRecursive(edge2, colorD, colorC, adjMapOrdered, colorMap);
                return;
            }
        }
    }

    /**
     Given an uncolored edge (u, v), the following procedure constructs a fan f of u
     that is maximal in that it cannot be extended.
     * @param edge
     * @param adjMapOrdered map w/ key = vertex, value= adjacent vertexes
     * @param cMap
     * @return
     */
    private static TIntList findMaximalFan(PairInt edge, TIntObjectMap<TIntList> adjMapOrdered, Map<PairInt, Integer> cMap) {

        if (cMap.containsKey(edge)) {
            throw new IllegalArgumentException("edge has already been colored");
        }

        int v = edge.getX();

        TIntSet inF = new TIntHashSet();
        TIntList f = new TIntArrayList();

        int xi = edge.getY();
        inF.add(xi);
        f.add(xi);

        TIntList vAdj = adjMapOrdered.get(v);
        if (vAdj == null || vAdj.isEmpty()) {
            return f;
        }
        int xii;
        int i;
        int cVXii;
        TIntList xiAdj;
        for (i = 0; i < vAdj.size(); ++i) {
            xii = vAdj.get(i);
            if (inF.contains(xii)) {
                continue;
            }
            // add xii if the color (v, xii) is not present in edges of xi
            // set xi = xii

            xiAdj = adjMapOrdered.get(xi);
            if (xiAdj == null) {
                f.add(xii);
                xi = xii;
                continue;
            }

            PairInt chk = makeOrderedNode(v, xii);
            if (!cMap.containsKey(chk)) {
                break;
            }

            cVXii = cMap.get(chk);

            //if xi edges do not have color of edge (v, xii), add xii
            if (isFree(xi, cVXii, adjMapOrdered, cMap)) {
                f.add(xii);
                xi = xii;
            }
        }
        return f;
    }

    private static boolean isFree(int u, int color, TIntObjectMap<TIntList> adjMapOrdered, Map<PairInt, Integer> cMap) {
        TIntList adj = adjMapOrdered.get(u);
        if (adj == null) {
            return true;
        }
        PairInt e;
        int v;
        for (int i = 0; i < adj.size(); ++i) {
            v = adj.get(i);
            e = makeOrderedNode(u, v);
            if (cMap.containsKey(e) && (cMap.get(e) == color)) {
                return false;
            }
        }
        return true;
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

    private static class SaturatedItem implements Comparable<SaturatedItem> {
        int sat;
        int deg;
        int vertex;
        public SaturatedItem(int saturation, int degree, int vertex) {
            this.sat = saturation;
            this.deg = degree;
            this.vertex = vertex;
        }
        @Override
        public int compareTo(SaturatedItem other) {
            if (sat > other.sat) {
                return 1;
            }
            if (sat < other.sat) {
                return -1;
            }
            if (deg > other.deg) {
                return 1;
            }
            if (deg < other.deg) {
                return -1;
            }
            // they're equal, so pick one
            if (vertex > other.vertex) {
                return 1;
            }
            return 0;
        }
    }
}
