package algorithms.optimization;

import algorithms.msts.PrimsMST;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
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
        ==> can use Hungarian algorithm
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
    
    public int[] approxTSPTour(final int nVertexes, final TIntObjectMap<TIntIntMap> adjCostMap) {
        
        //(1)
        int maxCost = PrimsMST.maxEdgeCost(adjCostMap);
        
        PrimsMST prims = new PrimsMST();
        
        prims.calculateMinimumSpanningTree(adjCostMap, maxCost);
        
        Map<Integer, LinkedList<Integer>> mstTree = prims.makeTreeFromPrev();
        
        //(2)
        int[] degrees = calculateDegrees(mstTree, nVertexes);
        int[] oddDVertexes = oddPassFilter(degrees);
        
        // building O as a cost matrix for input to the Hungarian algorithm
        float[][] oCostMatrix = buildCostMatrix(oddDVertexes, adjCostMap);
        
        // for min-cost perfect matching can use Hungarian Algorithm.
        // alternatively, can use MinCostUnbalancedAssignment.
        // hungarian algorithm accepts argument: computeAssignments(float[][] matrix).
        // MinCostUnbalancedAssignment needs a Graph g.
        HungarianAlgorithm ha = new HungarianAlgorithm();
        int[][] assignments = ha.computeAssignments(oCostMatrix);

        throw new UnsupportedOperationException("not yet finished");
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

        int i, j, u, v, uvCost;
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
                for (j = 0; j < neighborCost.size(); ++j) {
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
}
