package algorithms.graphs;

import algorithms.util.SimpleLinkedListNode;

/**
 * From Cormen et al. "Introduction to Algorithms",
 * Topological sort sorts a DAG (directed acyclic graph) by vertices such that
 * a directed edge uv to vertices u and v results in u before v in a linear
 * ordering.
 *   - call DFS(G) to compute finish times for each vertex v, f[v]
 *   - as each vertex is finished, insert it onto from of a linkedlist
 *   - return linked list of vertices
 * 
 * http://en.wikipedia.org/wiki/Topological_sorting
 * 
 * Good for dependency graphs or scheduling.
 *
 * A topological ordering is possible if and only if the graph has no directed
 * cycles, that is, if it is a directed acyclic graph (DAG).
 *
 * Runtime complexity is <em>O(V + E)</em>.
 * 
 * first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright

 * @author nichole
 */
public class TopologicalSort {

    /**
     * adjacency matrix with connected i->j indicated by the index and each
     *    node in the linked list, respectively.
     * for example, adjacent to node 3 is found via directedEdges[3] as all in the linked list.
     */
    protected final SimpleLinkedListNode[] directedEdges;

    /**
     * 
     * @param dag 
     */
    public TopologicalSort(SimpleLinkedListNode[] dag){
        if (dag == null || dag.length == 0) {
            throw new IllegalArgumentException("dag cannot be null or empty");
        }
        directedEdges = dag.clone();
        for (int i = 0; i < directedEdges.length; ++i) {
            directedEdges[i] = new SimpleLinkedListNode(dag[i]);
        }
    }
    
    /**
     * 
     * implemented following Cormen et al. "Introduction To Algorithms"
     * 
     * @return indexes ordered by finish time of traversal
     */
    public int[] sort() {
         //- call DFS(G) to compute finish times for each vertex v, f[v]
         //- as each vertex is finished, insert it onto front of a linkedlist
         // - return linked list of vertices
         
         //DFS dfs = new DFS(this.directedEdges);
         //dfs.walk();
         DFSIterative dfs = new DFSIterative();
         dfs.walk(directedEdges);
         int[] fIdxs = dfs.getOrderedEndIndexes();
         
         reverse(fIdxs);
        
         /*
         NOTE: some unit tests suggest that some implementations of topological sort
         next use partitioning of connected components and then further 
         sorts the results by the longest subsequences within the results,
         but does not change order for same subsequence length.
         */
         return fIdxs;
    }
    
    private void reverse(int[] a) {
        int idxLo = 0;
        int idxHi = a.length - 1;
        int n = idxHi - idxLo + 1;
        
        int end = idxLo + (n/2);
        
        int count = 0;
        for (int i = idxLo; i < end; i++) {
            int idx2 = idxHi - count;
            int swap = a[i];
            a[i] = a[idx2];
            a[idx2] = swap;
            count++;
        }
    }
    
}
