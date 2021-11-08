package algorithms.optimization;

import algorithms.graphs.HierholzersEulerCircuit;
import algorithms.msts.KruskalsMinimumSpanningTree;
import algorithms.msts.PrimsMST;
import algorithms.tsp.TSPPrimsMST;
import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import thirdparty.HungarianAlgorithm;

/**
 * An approximate solution to the Traveling Salesman Problem.
 * It is an approximation algorithm that guarantees that its solutions will be 
 * within a factor of 3/2 of the optimal solution length, and is named after 
 * Nicos Christofides and Anatoliy I. Serdyukov, who discovered it 
 * independently in 1976.  The runtime compexity is O(n^3)
 * 
 * following the pseudocode in
 * https://en.m.wikipedia.org/wiki/Christofides_algorithm
 * 
 * <pre>
 * (1) T=MST(G) where G is a complete graph with vertices v and non-negative edge weights w.
        (a complete graph is an undirected graph w/ edge between every pair of nodes)
        ==> can use PrimsMST
    (2) O = the vertices in T w/ odd degree.  in the subgraph O, connect all 
            vertices to one another.
    (3) M = min weight perfect matching in O
        ==> can use Hungarian algorithm or MinCostUnbalancedAssignment.java
    (4) H = connected multigraph from combining the edges of M and T, such that
        each vertex has even degree.
        (a multigraph may have more than 1 edge between same 2 end nodes).
    (5) EC = form a eulerian circuit
        ==> can use HierholzersEulerCircuit
    (6) Make the EC circuit into a Hamiltonian circuit by skipping repeated 
          vertices (shortcutting).
 * </pre>
 * @author nichole
 */
public class TSPChristofidesSerdyukov {
        
    /**
     * find a Hamiltonian tour of the given graph (simple cycle including all vertexes) 
     * that is 3/2 - approximate for minimum total cost.
     * @param nVertexes
     * @param adjCostMap
     * @return the Hamiltonian cycle within a factor of no more than 1.5 of the
     * optimal tour's minimum cost.  the array returned contains the vertex numbers
     * from the adjacency cost map.
     */
    public int[] approxTSPTour(final int nVertexes, final TIntObjectMap<TIntIntMap> adjCostMap) {
        
        //(1) T = mst(G)
        Map<Integer, LinkedList<Integer>> mstTree = buildMST(adjCostMap);
        
        //(2) W = odd degree vertices in T
        int[] degrees = calculateDegrees(mstTree, nVertexes);
        int[] oddDVertexes = oddPassFilter(degrees);
        // there are an even number of odd vertexes
        assert((oddDVertexes.length & 1) != 1);
        
        int i;
        
        //(3) O is subgraph of G induced by the vertices in oddDVertexes.
        //    make perfect min-cost matching from it
        // format: [nMatcings][2]
        int[][] m = bipartiteMinCostMatchingFromSubgraph(oddDVertexes, adjCostMap);
                
        //(4) H is the union of graphs T and M
        //    where each edges present in both T and M are present twice in H
        TIntObjectMap<TIntSet> h = unionMSTAndAssignments(mstTree, m);
        
        //(5) EC is the eulerian circuit in H. each edge is visited exactly once.
        HierholzersEulerCircuit hec = new HierholzersEulerCircuit();
        int[] ec = hec.createCircuit(h); // assuming start node = 0.

        //(6) T2 is the hamiltonian circuit of EC made by skipping over previously visited vertices.
        TIntSet visited = new TIntHashSet();
        TIntList t2 = new TIntArrayList();
        for (i = 0; i < ec.length; ++i) {
            if (!visited.contains(ec[i])) {
                t2.add(ec[i]);
                visited.add(ec[i]);
            }
        }
        t2.add(t2.get(0));
        int[] _t2 = t2.toArray();
        return t2.toArray();
    }
    
    /**
     * return an array of the indexes which have odd degrees.
     * @param degrees array where indexes are the vertex number and values are
     * the degree for the vertex.
     * @return array of indexes in degrees which have odd values stored in the degrees array.
     */
    protected int[] oddPassFilter(int[] degrees) {
        int[] odd = new int[degrees.length];
        int i, c = 0;
        for (i = 0; i < degrees.length; ++i) {
            if ((degrees[i] & 1) == 1) {
                odd[c] = i;
                c++;
            }
        }
        odd = Arrays.copyOf(odd, c);
        return odd;
    }

    protected int[] calculateDegrees(Map<Integer, LinkedList<Integer>> mstTree, 
        int nVertexes) {
        
        int[] d = new int[nVertexes];
        
        int u, v;
        LinkedList<Integer> neighbors;
        Iterator<Integer> iter = mstTree.keySet().iterator();
        Iterator<Integer> iterV;
        
        while (iter.hasNext()) {
            u = iter.next();
            neighbors = mstTree.get(u);
            
            iterV = neighbors.iterator();
            while (iterV.hasNext()) {
                v = iterV.next();
                d[u]++;
                d[v]++;
            }
        }
        
        return d;
    }

    /**
     * create a cost matrix from the vertexes listed in oddDVertexes where the
     * adjacency and costs are within adjCostMap.
     * @param oddDVertexes values are the vertex numbers with odd degrees
     * @param adjCostMap the cost map within the original adjacency map.
     * @return a cost matrix whose indexes are relative to oddDVertexes.
     * Note that non-existing connections have a cost of Float.MAX_VALUE.
     */
    protected float[][] buildCostMatrix(int[] oddDVertexes, 
        TIntObjectMap<TIntIntMap> adjCostMap) {

        int i, u, v, uvCost;
        float[][] out = new float[oddDVertexes.length][];
        for (i = 0; i < oddDVertexes.length; ++i) {
            out[i] = new float[oddDVertexes.length];
            Arrays.fill(out[i], Float.MAX_VALUE);
        }
  
        // a reverse index map to find where v is in oddDVertexes
        TIntIntMap rIdxs = new TIntIntHashMap();
        for (i = 0; i < oddDVertexes.length; ++i) {
            u = oddDVertexes[i];
            rIdxs.put(u, i);
        }
                
        int idxV;
        TIntIntMap neighborCost;
        TIntIntIterator iter;
        for (i = 0; i < oddDVertexes.length; ++i) {            
            u = oddDVertexes[i];
            neighborCost = adjCostMap.get(u);
            if (neighborCost != null) {
                iter = neighborCost.iterator();
                while (iter.hasNext()) {
                    iter.advance();
                    v = iter.key();
                    if (!rIdxs.containsKey(v)) {
                        continue;
                    }
                    uvCost = iter.value();
                    idxV = rIdxs.get(v);
                    out[i][idxV] = uvCost;
                }
            }
        }
        
        return out;
    }

    protected Map<Integer, LinkedList<Integer>> buildMST(TIntObjectMap<TIntIntMap> adjCostMap) {
                
        // finding the max cost in the graph G needed for a trie used in Prim's MST
        int maxCost = PrimsMST.maxEdgeCost(adjCostMap);
        
        PrimsMST prims = new PrimsMST();
        
        prims.calculateMinimumSpanningTree(adjCostMap, maxCost);
        
        //(1) T = mst(G)
        Map<Integer, LinkedList<Integer>> mstTree = prims.makeTreeFromPrev();
        //print(mstTree);
        
        //TIntList treeWalk = prims.getPreorderIndexes();
        //System.out.printf("treeWalk=%s\n", Arrays.toString(treeWalk.toArray()));
        
        return mstTree;
    }
    
    /*
    protected Map<Integer, LinkedList<Double>> buildMST2(TIntObjectMap<TIntDoubleMap> adjCostMap) {

        SimpleLinkedListNode[] graph, TObjectDoubleMap<PairInt> edgeWeights;
        
        TIntObjectMap<SimpleLinkedListNode> mst = KruskalsMinimumSpanningTree.mst(
            graph, edgeWeights);
        
        return mstTree;
    }*/

    /**
     * perfect min-cost bipartite matchings of the subgraph of G induced by the
     * odd vertexes.  The results are in a double array where each row
     * is a pair of matching vertexes in context of graph G.
     * @param oddDVertexes
     * @param adjCostMap
     * @return 
     */
    protected int[][] bipartiteMinCostMatchingFromSubgraph(
        int[] oddDVertexes, TIntObjectMap<TIntIntMap> adjCostMap) {
        
        // building O as a cost matrix for input to the Hungarian algorithm
        float[][] oCostMatrix = buildCostMatrix(oddDVertexes, adjCostMap);
        // indexes of the oCostMatrix are the values of oddDVertexes
        //   so should be transformed back to vertex numbers after
        //   mincost matching.
                
        // for min-cost perfect matching can use Hungarian Algorithm.
        // alternatively, can use MinCostUnbalancedAssignment.
        // hungarian algorithm accepts argument: computeAssignments(float[][] matrix).
        // MinCostUnbalancedAssignment needs a Graph g.
        int[][] assignmentsM = new HungarianAlgorithm().computeAssignments(oCostMatrix);
        
        int i, n = 0, k, v;
        TIntIntMap keep = new TIntIntHashMap();
        TIntSet included = new TIntHashSet();
        for (i = 0; i < assignmentsM.length; ++i) {
            k = assignmentsM[i][0];
            v = assignmentsM[i][1];
            //System.out.printf("%d, %d\n", k, v);
            if (included.contains(k) || included.contains(v)) {
                continue;
            }
            keep.put(k, v);
            included.add(k);
            included.add(v);
            n++;
        }
        
        int[][] m = new int[n][];
        TIntIntIterator iterK = keep.iterator();
        for (i = 0; i < n; ++i) {
            iterK.advance();
            k = iterK.key();
            v = iterK.value();
            m[i] = new int[]{oddDVertexes[k], oddDVertexes[v]};
        }
        
        return m;
    }

    protected TIntObjectMap<TIntSet> unionMSTAndAssignments(
        Map<Integer, LinkedList<Integer>> mstTree, int[][] m) {
        
        int i, u, v;
        
        TIntObjectMap<TIntSet> h = new TIntObjectHashMap<TIntSet>();
        
        Iterator<Integer> iterMST = mstTree.keySet().iterator();
        LinkedList<Integer> lList;
        Iterator<Integer> iterList;
        TIntSet s;
        while (iterMST.hasNext()) {
            u = iterMST.next();
            lList = mstTree.get(u);
            iterList = lList.iterator();
            while (iterList.hasNext()) {
                v = iterList.next();
                s = h.get(u);
                if (s == null) {
                    s = new TIntHashSet();
                    h.put(u, s);
                }
                s.add(v);
            }
        }

        for (i = 0; i < m.length; ++i) {
            s = h.get(m[i][0]);
            if (s == null) {
                s = new TIntHashSet();
                h.put(m[i][0], s);
            }
            s.add(m[i][1]);
        }        
        
        return h;
    }
    
    public static long totalCost(int[] hamiltonian, TIntObjectMap<TIntIntMap> adjCostMap) {
        long sum = 0;
        int i, u, v, cost;
        TIntIntMap assoc;
        for (i = 1; i < hamiltonian.length;++i) {
            u = hamiltonian[i - 1];
            v = hamiltonian[i];
            assoc = adjCostMap.get(u);
            if (assoc == null) {
                throw new IllegalStateException("node " + u + " is not a key in the adjCostMap");
            }
            if (!assoc.containsKey(v)) {
                throw new IllegalStateException("node " + v + " is not a value in the adjCostMapfor key=" + u);
            }
            sum += assoc.get(v);
        }
        return sum;
    }

    private void print(Map<Integer, LinkedList<Integer>> mstTree) {
        Iterator<Integer> intIter = mstTree.keySet().iterator();
        Iterator<Integer> intIter2;
        int u, v;
        System.out.println("MST:");
        while (intIter.hasNext()) {
            u = intIter.next();
            intIter2 = mstTree.get(u).iterator();
            while (intIter2.hasNext()) {
                v = intIter2.next();
                System.out.printf("%d:%d\n", u, v);
            }
        }
    }
}
