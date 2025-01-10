package algorithms.graphs;

import algorithms.matrix.MatrixUtil;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;

import java.util.*;

/**
 Hierholzer’s Algorithm to create a euler circuit from a directed graph.
 Euler circuit is a path that traverses every edge of a graph, and the path
 ends on the starting vertex.  Edges are only included once, but vertexes can be included more than once.

 Note that if the input graph is a tree, the result is the same as
 a tree pre-order traversal when starting from the same source node.

 <pre>
 The implementation follows from:
 https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/
 written by
 https://in.linkedin.com/in/ashutosh-kumar-9527a7105

   Geeks For Geeks  has a liberal copyright policy at:
   https://www.geeksforgeeks.org/copyright-information/?ref=footer
   an excerpt:
   You are free to:
    Share — copy and redistribute the material in any medium or format
    Adapt — remix, transform, and build upon the material for any purpose.

    Under the following terms:

    Attribution — You must give appropriate credit, provide a link to the 
      license, and indicate if changes were made. You may do so in any reasonable manner, 
      but not in any way that suggests the licensor endorses you or your use.

    Hyperlink each article directly back to their user profile page on the source site.
    * 
 * </pre>
 * @author nichole
 */
public class HierholzersEulerCircuit {
    
    /**
     create a euler circuit from a directed graph.
     Euler circuit is a path that traverses every edge of a graph, and the path 
     ends on the starting vertex.
     Note that if given a tree with directed edges that do not return to start,
     the result is a Euler Tour defined as a pre-order traversal of the tree starting
     from startNode=0.
     @param g adjacency list of directed graph
     @return euler circuit
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
     create a euler circuit from a directed graph.
     Euler circuit is a path that traverses every edge of a graph, and the path 
     ends on the starting vertex.
     Note that if given a tree with directed edges that do not return to start,
     the result is a Euler Tour defined as a pre-order traversal of the tree starting
     from startNode.
     @param g adjacency list of directed graph
     @param startNode
     @return euler circuit
     */
    public int[] createCircuit(TIntObjectMap<TIntSet> g, int startNode) {
        if (g.isEmpty()) {
            return new int[0];
        }
               
        // make a copy of g to modify
        TIntObjectMap<TIntSet> g2 = MatrixUtil.copy(g);

        TIntList circuit = new TIntArrayList();

        Stack<Integer> curPath = new Stack<Integer>();
        curPath.add(startNode);

        while (!curPath.isEmpty()) {
            TIntSet neighbors2 = g2.get(curPath.peek());
            if (neighbors2 != null && !neighbors2.isEmpty()) {
                int nextV = neighbors2.iterator().next();
                neighbors2.remove(nextV);
                curPath.add(nextV);
            } else {
                circuit.add(curPath.pop());
            }
        }
        
        circuit.reverse();
        
        return circuit.toArray();
    }

    /**
     create a euler circuit from a directed graph.
     Euler circuit is a path that traverses every edge of a graph, and the path
     ends on the starting vertex.
     Note that if given a tree with directed edges that do not return to start,
     the result is a Euler Tour defined as a pre-order traversal of the tree starting
     from startNode.
     @param g adjacency list of directed graph
     @param startNode
     @return a 2-dimensional array where row 0 is the euler circuit
     and row 1 is the distance of each node from startNode (which has dist=0)
     */
    public int[][] createCircuitAndDepth(Map<Integer, LinkedList<Integer>> g, int startNode) {
        if (g.isEmpty()) {
            return new int[0][0];
        }

        // make a copy of g to modify
        Map<Integer, LinkedList<Integer>> g2 = GraphUtil.copy2(g);

        TIntList circuit = new TIntArrayList();
        Map<Integer, Integer> distMap = new HashMap<>();
        distMap.put(startNode, 0);

        Stack<Integer> curPath = new Stack<Integer>();
        curPath.add(startNode);

        while (!curPath.isEmpty()) {
            LinkedList<Integer> neighbors2 = g2.get(curPath.peek());
            if (neighbors2 != null && !neighbors2.isEmpty()) {
                curPath.add(neighbors2.pollFirst());
            } else {
                int node = curPath.pop();

                // calc dist
                int prevDist = 0;
                if (!circuit.isEmpty()) {
                    prevDist = distMap.get(circuit.get(circuit.size() - 1));
                }
                if (distMap.containsKey(node)) {
                    distMap.put(node, Math.min(distMap.get(node), prevDist + 1));
                } else {
                    distMap.put(node, prevDist + 1);
                }

                circuit.add(node);
            }
        }

        circuit.reverse();
        int[][] out = new int[2][];
        out[0] = circuit.toArray();
        out[1] = new int[out[0].length];
        for (int i = 0; i < out[1].length; ++i) {
            out[1][i] = distMap.get(out[0][i]);
        }
        return out;
    }
}
