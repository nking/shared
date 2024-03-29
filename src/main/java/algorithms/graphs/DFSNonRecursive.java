package algorithms.graphs;

import algorithms.sort.MiscSorter;
import algorithms.util.SimpleLinkedListNode;
import java.util.Arrays;
import java.util.Stack;

/**
   DFS

   searches the full depth of a graph or subgraph when possible first then
      backtracks to the unexplored edges and unexplored nodes repeating until
      all nodes are visited.  unlike BFS, it may contain many predecessor trees, 
      that is a predecessor forest of nodes that are the shortest from the 
      source to each reachable node.  for this reason, DFS searches can need a 
      lot of memory.

   average runtime is approx O(|E|), worst case runtime: O(|V| + |E|)
   worst case space needed: O(|V|)

   implemented using recursion following Cormen, Leiserson, Rivest, and Stein "Introduction To Algorithms"
   then re-factored to non-recursive by using a stack, following advice in 
   https://www.codeproject.com/Articles/418776/How-to-replace-recursive-functions-using-stack-and

   first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright

   A note on recursion from lecture notes of 
   CS473: Fundamental Algorithms (Spring 2011)
   https://courses.engr.illinois.edu/cs473/sp2011/Lectures/09_lec.pdf
   An iterative algorithm B obtained from a recursive algorithm A for a problem Π 
   does the following: for each instance I of Π, it computes a topological sort 
   of G(I) and evaluates sub-problems according to the topological ordering.
   In some cases (not all) the computation of an optimal solution reduces to a 
   shortest/longest path in DAG G(I).
* 
 * @author nichole
 */
public class DFSNonRecursive {
    /**
     * adjacency matrix with connected i to j indicated by the index and each
     *    node in the linked list, respectively.
     * for example, adjacent to node 3 is found via directedEdges[3] as all in the linked list.
     */
    private SimpleLinkedListNode[] g;
    
    /** 
     * holds state for whether a node has been visited.  0 = not visited,
     * 1 = visiting now, 2 = was visited.
    */
    private int[] visited;

    /**
     * time when node is first discovered
     */
    private int[] td;

    /**
     * time when node's adjacency list has all been visited
     */
    private int[] tf;
   
    private int[] predecessor;

    private int time;

    /**
     *
     */
    public DFSNonRecursive() {
        
    }

    /**
     @param directedEdges  adjacency matrix with connected i to j indicated
     * by the index and each node in the linked list, respectively.
     * Note that the key of each node is expected to be the same as it's index
     * in the adjacency matrix.
     * For example, adjacent to node 3 is found via directedEdges[3] as all in 
     * the linked list.
     */
    public void walk(SimpleLinkedListNode[] directedEdges) {
        if (directedEdges == null || directedEdges.length == 0) {
            throw new IllegalArgumentException("directedEdges cannot be null or empty");
        }
        g = directedEdges.clone();
        for (int i = 0; i < g.length; ++i) {
            g[i] = new SimpleLinkedListNode(directedEdges[i]);
        }
        visited = new int[g.length];
        td = new int[g.length];
        tf = new int[g.length];
        predecessor = new int[g.length];
        Arrays.fill(td, -1);
        Arrays.fill(tf, -1);
        Arrays.fill(predecessor, -1);
        time = 0;
        
        for (int u = 0; u < g.length; u++) {
            if (visited[u] == 0) {
                walk(u);
            }
        }
    }
    
    private void walk(int u) {
        
        Stack<Snapshot> stack = new Stack<Snapshot>();
        Snapshot current;
        
        //System.out.println("*load method frame for " + u);
        
        current = new Snapshot(u);
        current.stage = 0;
        stack.push(current);
        
        while(!stack.empty()) {
            
            current = stack.pop();
            
            //System.out.println(current.toString());
            
            switch(current.stage) {
                case 0: { 
                    // before recursion is invoked
                    visited[current.node] = 1;
                    time++;
                    //System.out.println("  0: visiting " + current.node + " to set td=" + time);
                    td[current.node] = time;
                    
                    current.stage = 1;
                    stack.push(current);
                    
                    //System.out.format("  0: push onto stack u=%d\n", current.node);
                            
                    SimpleLinkedListNode next = g[current.node];
                    
                    if (next != null && next.getNumberOfKeys() > 0) {
                        
                        int v = next.getKey();
                        
                        g[current.node].delete(next);
                                                      
                        if (visited[v] == 0) {
                            
                            predecessor[v] = current.node;
                            
                            Snapshot newSnapshot = new Snapshot(v);
                            newSnapshot.stage = 0;
                            stack.push(newSnapshot);

                            //System.out.format("   0: and push onto stack v=%d\n", v);
                            //System.out.println("   0: [v: " + newSnapshot.toString() + "]");
  
                            continue;
                        } else if (predecessor[v] == -1) {
                            // in case the instance graph is not ordered top-down
                            // this is a back edge.
                            predecessor[v] = current.node;
                        }
                    }
                    break;
                }
                case 1: {
                    //System.out.println(" 1: have all child links been visited?  snap="
                    //   + current.toString());
                    
                    SimpleLinkedListNode next = g[current.node];
                    if (next != null && next.getNumberOfKeys() > 0) {
                        
                        int v = next.getKey();
                        
                        //System.out.format(" 1: there is a child link %d\n", v);
                        
                        g[current.node].delete(next);
                        
                        current.stage = 1;
                        stack.push(current);

                        //System.out.format("  0: push onto stack u=%d\n", current.node);
                                                      
                        if (visited[v] == 0) {
                            
                            predecessor[v] = current.node;
                            
                            Snapshot newSnapshot = new Snapshot(v);
                            newSnapshot.stage = 0;
                            stack.push(newSnapshot);

                            //System.out.format("   1: and push onto stack v=%d\n", v);
                            //System.out.println("   1: [v: " + newSnapshot.toString() + "]");
  
                            continue;
                        } else if (predecessor[v] == -1) {
                            //this is a back edge
                            predecessor[v] = current.node;
                        }
                        
                        continue;
                    }
                    
                    visited[current.node] = 2;
                    time++;
                    tf[current.node] = time;
                    //System.out.format(" 1: end visit to %d, set tf=%d\n",
                    //    current.node, time);

                    break;
                }
            }
        }
    }
    
    private class Snapshot {
        
        /**
         * index of current snapshot within DFSIterative instance's arrays.
         */
        protected final int node;
                
        protected int stage = 0;
                        
        public Snapshot(int u) {
            this.node = u;
        }
                
        public Snapshot(Snapshot s) {
            this.stage = s.stage;
            this.node = s.node;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("node=").append(Integer.toString(node))
                .append(", stage=").append(Integer.toString(stage))
                .append(", prev=").append(Integer.toString(predecessor[node]))
                .append(", visited=").append(Integer.toString(visited[node]))
            ;
            return sb.toString();
        }
        
    }
    
    /**
     * get predecessor indexes
     @return get predecessor indexes
     */
    public int[] getPredecessorIndexes() {
        if (predecessor == null) {
            return null;
        }
        return Arrays.copyOf(predecessor, predecessor.length);
    }
    
    /**
     * return the indexes in order of the starts of their traversals
     @return 
     */
    public int[] getOrderedBeginIndexes() {
        if (td == null) {
            return null;
        }
        return sortForIndexes(td);
    }
    
    private int[] sortForIndexes(int[] a) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (g == null) {
            return null;
        }
        assert(a.length == g.length);
        a = Arrays.copyOf(a, a.length);
        int[] idxs = new int[a.length];
        for (int i = 0; i < idxs.length; ++i) {
            idxs[i] = i;
        }
        MiscSorter.sortBy1stArg(a, idxs);
        return idxs;
    }
    
    /**
     * return the indexes in order of the ends of their traversal
     @return 
     */
    public int[] getOrderedEndIndexes() {
        if (tf == null) {
            return null;
        }
        return sortForIndexes(tf);
    }
    
    /**
     *
     @return
     */
    public int[] getTd() {
        if (td == null) {
            return null;
        }
        return td;
    }

    /**
     *
     @return
     */
    public int[] getTf() {
        if (tf == null) {
            return null;
        }
        return tf;
    }
    
}
