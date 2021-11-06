package algorithms.graphs;

import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Stack;

/**
 Hierholzer’s Algorithm to create a euler circuit from a directed graph.
  Euler circuit is a path that traverses every edge of a graph, and the path 
  ends on the starting vertex.
 * <pre>
 * The implementation follows from:
 * https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/
 * written by
 * https://in.linkedin.com/in/ashutosh-kumar-9527a7105
 *
   Geeks For Geeks  has a liberal copyright policy at:
   https://www.geeksforgeeks.org/copyright-information/?ref=footer
   an excerpt:
   You are free to:
    Share — copy and redistribute the material in any medium or format
    Adapt — remix, transform, and build upon the material for any purpose.

    Under the following terms:

    Attribution — You must give appropriate credit, provide a link to the 
      license, and indicate if changes were made. You may do so in any r
      easonable manner, but not in any way that suggests the licensor endorses 
      you or your use.

    Hyperlink each article directly back to their user profile page on the source site.
    * 
 * </pre>
 * @author nichole
 */
public class HierholzersEulerCircuit {
    
    /**
     * 
     * @param g adjacency list of directed graph
     * @return 
     */
    public int[] createCircuit(TIntObjectMap<TIntSet> g) {
        if (g.isEmpty()) {
            return new int[0];
        }
        int startNode = 0;
        // OR, in the case that 0 is not connected
        //startNode = g.keySet().iterator().next();
        
        return createCircuit(g, startNode);
    }
    /**
     * 
     * @param g adjacency list of directed graph
     * @param startNode
     * @return 
     */
    public int[] createCircuit(TIntObjectMap<TIntSet> g, int startNode) {
        if (g.isEmpty()) {
            return new int[0];
        }
       
        int i, u;
        
        // make a copy of g to modify
        TIntObjectIterator<TIntSet> iter = g.iterator();
        TIntSet neighbors, neighbors2;
        TIntObjectMap<TIntSet> g2 = new TIntObjectHashMap<TIntSet>();
        for (i = 0; i < g.size(); ++i) {
            iter.advance();
            u = iter.key();
            neighbors = iter.value();
            neighbors2 = new TIntHashSet(neighbors);
            g2.put(u, neighbors2);
        }
                
        Stack<Integer> curPath = new Stack<Integer>();
        TIntList circuit = new TIntArrayList();
        
        // start vertex
        int curV = startNode;
        
        curPath.add(curV);
        
        int nextV;
        
        while (!curPath.isEmpty()) {
            curV = curPath.peek();
            neighbors2 = g2.get(curV);
            if (neighbors2 != null && !neighbors2.isEmpty()) {
                nextV = neighbors2.iterator().next();
                neighbors2.remove(nextV);
                curPath.add(nextV);
            } else {
                circuit.add(curPath.pop());
            }
        }
        
        circuit.reverse();
        
        return circuit.toArray();
    }
}
