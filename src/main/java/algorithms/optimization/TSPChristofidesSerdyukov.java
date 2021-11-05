package algorithms.optimization;

import algorithms.msts.PrimsMST;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;

/**
 * An approximate solution to the Traveling Salesman Problem.
 * It is an approximation algorithm that guarantees that its solutions will be 
 * within a factor of 3/2 of the optimal solution length, and is named after 
 * Nicos Christofides and Anatoliy I. Serdyukov, who discovered it 
 * independently in 1976.
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
        
        int[] degrees = calculateDegrees(mstTree, nVertexes);
        
        //(2)

        throw new UnsupportedOperationException("not yet finished");
    }

    private int[] calculateDegrees(Map<Integer, LinkedList<Integer>> mstTree, int nVertexes) {
        
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
}
