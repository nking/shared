package algorithms.graphs;

import algorithms.sort.MiscSorter;
import algorithms.util.PairIntArray;
import algorithms.util.SimpleLinkedListNode;
import java.util.Arrays;

/**
 * Find all bridges using DFS.
 * 
 * An edge (u,v) is a bridge if and only if it is a tree edge and (assuming 
 * that u is the parent of v) there is no back edge within v’s subtree that 
 * leads to a vertex whose discovery time is strictly smaller than v’s 
 * discovery time.

   definitions:
   Bridge is any edge whose removal results in a disconnected graph.
   tree edges connect a parent with its child in the DFS tree.
   back edges connect a (non-parent) ancestor with a (non-child) descendant.
   2-Edge Connected: A graph is 2-edge connected if it contains no bridges
    
   average runtime is approx O(|E|), worst case runtime: O(|V| + |E|)
   worst case space needed: O(|V|).

   References:
     * <pre>
     * lecture 4 notes of David Mount for CMSC 451 
     * Design and Analysis of Computer Algorithms (with some corrections for pseudocode indexes).
     * https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect04-edge-connectivity.pdf
     * </pre>
* 
* explanation:
* Suppose that we are currently processing a vertex u in DFSvisit, and we see 
* an edge (u,v) going to a neighbor v of u. 
* If this edge is a back edge (that is, if v is an ancestor of u) 
* then (u,v) cannot be a bridge, because the tree edges between u and v provide 
* a second way to connect these vertices. Therefore, we may limit consideration 
* to when (u,v) is a tree edge, that is, v has not yet been discovered, 
* and so we will invoke DFSvisit(v). While we are doing this, we will keep 
* track of the back edges in the subtree rooted at v.
 *
 * An edge (u,v) is a bridge if and only if it is a tree edge and (assuming 
 * that u is the parent of v) there is no back edge within v’s subtree that 
 * leads to a vertex whose discovery time is strictly smaller than v’s 
 * discovery time.
 */
public class FindAllBridgesDFS {

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
     * Low[u] is the closest to the root that you can get in the tree by 
     * taking any one back edge from either u or any of its descendants.
     */
    protected final int[] tdLow;

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
    public FindAllBridgesDFS(SimpleLinkedListNode[] directedEdges) {
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
        this.tdLow = new int[g.length];
    }

    /**
     * find all bridges using DFS.
     @return the bridges as (u, v) pairs of edges.
     */
    public PairIntArray walk() {
        int u, v;
        for (u = 0; u < g.length; u++) {
            if (visited[u] == 0) {
                visit(u);
            }
        }
        PairIntArray b = new PairIntArray();
        for (v = 0; v < g.length; v++) {
            if ((td[v] == tdLow[v]) && predecessor[v] != -1) {
                // found a bridge
                b.add(predecessor[v], v);
            }
        }
        return b;
    }
    
    /**
     * alternative pattern for walking code
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
        tdLow[u] = td[u];

        SimpleLinkedListNode next = g[u];
        
        while (next != null && next.getNumberOfKeys() > 0) {
            int v = next.getKey();
            if (visited[v] == 0) {
                predecessor[v] = u;
                visit(v);
                // tree edge
                tdLow[u] = Math.min(tdLow[u], tdLow[v]);  // update Low[u]
            } else if (predecessor[u] != v) {
                // back edge
                tdLow[u] = Math.min(tdLow[u], td[v]);       // update Low[u]
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
