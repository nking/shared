package algorithms.graphs;

import algorithms.util.SimpleLinkedListNode;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Queue;

/**
 * From Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms",
 * Topological sort sorts a DAG (directed acyclic graph) by vertices such that
 * a directed edge uv to vertices u and v results in u before v in a linear
 * ordering.
 *   - call DFS(G) to compute finish times for each vertex v, f[v]
 *   - as each vertex is finished, insert it into front of a linkedlist
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
     * adjacency matrix with connected i to j indicated by the index and each
     *    node in the linked list, respectively.
     * for example, adjacent to node 3 is found via directedEdges[3] as all in the linked list.
     */
    protected final SimpleLinkedListNode[] directedEdges;

    /**
     * 
     @param dag 
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
     * implemented following Cormen, Leiserson, Rivest, and Stein "Introduction To Algorithms".
     * Runtime complexity is <em>O(V + E)</em>.
     @return indexes ordered by finish time of traversal
     */
    public int[] sort() {
         //- call DFS(G) to compute finish times for each vertex v, f[v]
         //- as each vertex is finished, insert it onto front of a linkedlist
         // - return linked list of vertices
         
         //DFS dfs = new DFS(this.directedEdges);
         //dfs.walk();
         DFSNonRecursive dfs = new DFSNonRecursive();
         dfs.walk(directedEdges);
         int[] fIdxs = dfs.getOrderedEndIndexes();
         
         reverse(fIdxs);
        
         return fIdxs;
    }
    
    private void reverse(int[] a) {
        int idxLo = 0;
        int idxHi = a.length - 1;
        int n = idxHi - idxLo + 1;        
        int swap;
        int idx2 = idxHi;
        for (int i = idxLo; i < (idxLo + (n/2)); i++) {
            swap = a[i];
            a[i] = a[idx2];
            a[idx2] = swap;
            idx2--;
        }
    }

    /**
     * use Kahn's method to solve for the topoligcal sorting.
     <pre>
     https://en.m.wikipedia.org/wiki/Topological_sorting
     https://www.interviewkickstart.com/learn/kahns-algorithm-topological-sorting
     </pre>
     * @return the topologically sorted keys of the DAG, or null if a cycle was detected
     */
    public int[] sortKahn() {

        /*
        (0) initialize an empty array "out" to hold results
        (1) make an array
             of inDegree in O(n) by traversing each vertec
        (2) initialize a queue
             with all vertexes that have inDegree = 0
        (3) while !q.isEmpty()
             index = q.poll()
             write index to out
             for all neighbors of index:
                 reduce the neigbhor inDeg by 1.
                 if their inDeg is now 0, q.offer( neighbor )
        (4) if out is not full, return null,
            else return out
         */

        SimpleLinkedListNode[] dag = directedEdges.clone();
        int n = dag.length;

        int i, v;
        SimpleLinkedListNode next;
        int[] inDeg = new int[n];
        for (i = 0; i < n; ++i) {
            next = dag[i];
            while (next != null && next.getNumberOfKeys() > 0) {
                v = next.getKey();
                inDeg[v]++;
                next = next.getNext();
            }
        }

        // prime the queue with nodes that have no incoming edges
        Queue<Integer> q = new ArrayDeque<>();
        for (i = 0; i < n; ++i) {
            if (inDeg[i] == 0) q.add(i);
        }

        int[] r = new int[n];
        Arrays.fill(r, -1);
        int iR = -1;
        while (!q.isEmpty()) {
            i = q.poll();
            r[++iR] = i;

            next = dag[i];
            while (next != null && next.getNumberOfKeys() > 0) {
                v = next.getKey();
                --inDeg[v];
                if (inDeg[v] == 0) {
                    q.add(v);
                }
                next = next.getNext();
            }
            dag[i] = null;
        }

        if (iR != (n-1)) {
            return null;
        }

        return r;
    }
    
}
