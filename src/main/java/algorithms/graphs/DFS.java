package algorithms.graphs;

import algorithms.sort.MiscSorter;
import algorithms.util.SimpleLinkedListNode;
import java.util.Arrays;
import java.util.stream.IntStream;

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

   implemented following Cormen, Leiserson, Rivest, and Stein "Introduction To Algorithms"

   pros of DFS:
      useful in topological sorting
      connected components, espec 2-(edge or vertex) and 3-(edge or vertex).
      finding the bridges of a graph
      Generating words in order to plot the limit set of a group.
      Finding strongly connected components.
      Planarity testing
      Solving puzzles with only one solution, such as mazes. (
      biconnectivity in graphs.
  
   cons of DFS:
      difficult to parallelize.
      requires more memory, in C++ too because it's not using tail recursion
      in its recursion.
  
 first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright

* 
 * @author nichole
 */
public class DFS {

    /**
     * adjacency matrix with connected i to j indicated by the index and each
     *    node in the linked list, respectively.
     * for example, adjacent to node 3 is found via directedEdges[3] as all in the linked list.
     */
    protected final SimpleLinkedListNode[] g;

    /** 
     * holds state for whether a node has been visited.  0 = not visited,
     * 1 = visiting now, 2 = was visited.
    */
    protected final int[] visited;

    /**
     * time when node is first discovered
     */
    protected final int[] td;

    /**
     * time when node's adjacency list has all been visited
     */
    protected final int[] tf;
   
    /**
     *
     */
    protected final int[] predecessor;

    /**
     *
     */
    protected int time = 0;
    
    /**
     @param directedEdges  adjacency matrix with connected i to j indicated 
     * by the index and each node in the linked list, respectively.
     * Note that the key of each node is expected to be the same as it's index
     * in the adjacency matrix.
     * For example, adjacent to node 3 is found via directedEdges[3] as all in 
     * the linked list.
     */
    public DFS(SimpleLinkedListNode[] directedEdges) {
        if (directedEdges == null) {
            throw new IllegalArgumentException("directedEdges cannot be null");
        }
        g = directedEdges.clone();
        for (int i = 0; i < g.length; ++i) {
            g[i] = new SimpleLinkedListNode(directedEdges[i]);
        }
        this.visited = new int[g.length];
        this.td = new int[g.length];
        this.tf = new int[g.length];
        this.predecessor = new int[g.length];
        Arrays.fill(predecessor, -1);
    }

    /**
     *
     */
    public void walk() {

        for (int u = 0; u < g.length; u++) {
            if (visited[u] == 0) {
                visit(u);
            }
        }        
    }

    /**
     * walk only the paths reachable from src
     * @param src
     */
    public void walk(int src) {
        visit(src);
    }
    
    /**
     * alterative pattern for walking code
     @param vertexOrder
     */
    void _walk(int[] vertexOrder) {
        for (int u : vertexOrder) {
            if (visited[u] == 0) {
                visit(u);
            }
        }
    }
    
    private void visit(int u) {
        //System.out.println("load method frame for " + u);
        
        visited[u] = 1;
        time++;
        //System.out.println("  visiting " + u + " to set td=" + time);
        td[u] = time;

        SimpleLinkedListNode next = g[u];
        
        while (next != null && next.getNumberOfKeys() > 0) {
            int v = next.getKey();
            if (visited[v] == 0) {
                predecessor[v] = u;
                visit(v);
            } else if (predecessor[v] == -1) {
                // visited v because of order invoked from walk, but predecessor of v was not visited before v
                predecessor[v] = u;
            }
            next = next.getNext();
        }
        visited[u] = 2;
        time++;
        tf[u] = time;
        //System.out.println("  visited " + u + ") to set tf=" + time);
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
        return sortForIndexes(tf);
    }

    public int[] getEulerTourFromEndTimes() {
        // decending sort
        int[] sIdxs = IntStream.range(0, tf.length).boxed()
                .sorted((a, b) -> {
                    int c = tf[b] - tf[a];
                    if (c != 0) return c;
                    return a - b;
                }).mapToInt(ele->ele).toArray();

        return sIdxs;
    }
    
    /**
     *
     @return
     */
    public int[] getTd() {
        return td;
    }

    /**
     *
     @return
     */
    public int[] getTf() {
        return tf;
    }
    
}
