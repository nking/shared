package algorithms.graphs;

import algorithms.bipartite.MinHeapForRT2012;
import algorithms.heapsAndPQs.HeapNode;
import algorithms.misc.MiscMath0;
import algorithms.sort.MiscSorter;
import algorithms.util.SimpleLinkedListNode;
import java.util.Arrays;

/**
 * From pseudocode in "Introduction to Algorithms" by Cormen, Leiserson, Rivest, and Stein
 * (impl is Tarjan's algorithm).
 * 
 * Directed graphs are said to be strongly connected if every vertex is 
 * reachable from every other vertex and there is a path in each 
 * direction between each pair of vertices of the graph.
 * 
 * @author nichole
 */
public class StronglyConnectedComponents2 {
        
    /**
     * Given the adjacency list for a graph whose vertexes are the indexes of 
     * the list, return the strongly connected components (the strongly 
     * connected components are a DAG).
     * 
     * RT complexity for worse case is O(|V| + |E|).
     * 
     @param connected adjacency of connected components
     @return an array of component numbers that correspond to the indexes of 
     * connected.
     * e.g. [0,1,1,0,2] says that vertexes 0 and 3 are strongly connected, 
     * and vertexes 1 and 2 are strongly connected...
     */
    public int[] findStronglyConnectedComponents(SimpleLinkedListNode[] 
        connected) {
        
        /*
        1) call DFS(g) to compute finishing times f[u] for each vertex u
        2) compute g^T (where g^T is g with edges reversed)
        3) call DFS(g^T) but in the main loop of DFS, consider the vertices in
           order of decreasing f[u] (as computed in line (1))
        4) output the vertices of each tree in the depth first forest formed in
           line (3) as a separate strongly connected component.
        */
        DFS dfs = new DFS(connected);
        dfs.walk();
           
        SimpleLinkedListNode[] reversedEdges = reverse(connected);
        
        int[] rTf = Arrays.copyOf(dfs.tf, dfs.tf.length);
        int[] indexes = new int[rTf.length];
        for (int i = 0; i < rTf.length; ++i) {
            indexes[i] = i;
        }
        MiscSorter.sortByDecr(rTf, indexes);
              
        // DFS traversing vertices rTF on graph gT
        DFS dfs2 = new DFS(reversedEdges);
        dfs2._walk(indexes);
        
        /* components can be found by looking at the G' start and end times.
        The components are between them.
        
        For example, see Fig 22.9 from Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms.
        The components are indexes c0={0,1,4}, c1={2,3}, c2={5,6}, c3={7}
        
        td is time when node is first discovered in dfs search
        tf is finish time
        
                indexes= 0   1   2   3   4   5   6   7
                G  td=[  1,  2, 11, 12,  9,  3,  4,  5]
        ===>    G' td=[  1,  3,  7,  8,  2, 11, 12,  15] <===
                G  tf=[ 16, 15, 14, 13, 10,  8,  7,   6]
        ===>    G' tf=[  6,  4, 10,  9,  5, 14,  13, 16] <====

        looking at G' td and tf:
           c0 starts at t=1:6, so that is indexes 0,4,1 
           c1 starts at t=7:10 is indexes 2,3
           c2 starts at t=11:14 is indexes 5,6
           c3 starts at t=15:16 is index 7
        */
        
        //ascending sort by dfs2.td, then select over ascending
        // for a list:  sort by mergesort or quicksort is O(N*log_2(N)),  select is O(1)
        // for a minheap (yfasttrie): sort is N inserts of O(log_2(w)) + O(w-l)
        //                            where w is the number of bits used for a data
        //                            and l is the prefix tree already filled
        //                            leading up to the value node.
        //                           and select is extractMin which is O(log_2(w)) + O(w-l).
        // so if use a data size of 16 bits:
        //     for n<20 use a list, else use a yfastrie
        // if use a data size of 32 bits:
        //     for n<32 use a list, else use a yfasttrie
        
        // using a yfasttrie with key=start time.  
        //   the yfasttrie extractMin() has constant small runtime complexity
        int maxV = MiscMath0.findMax(dfs2.td);
        int approxN = dfs2.td.length;
        int maxNumberOfBits = (int)Math.ceil(Math.log(maxV)/Math.log(2));
        MinHeapForRT2012 heap = new MinHeapForRT2012(maxV, approxN,
            maxNumberOfBits);
        for (int i = 0; i < dfs2.td.length; ++i) {
            // td is time when node is first discovered in dfs search
            HeapNode node = new HeapNode(dfs2.td[i]);
            NodeData d = new NodeData();
            d.idx = i;
            d.tf = dfs2.tf[i];
            node.setData(d);
            heap.insert(node);
        }
   
        int c = -1;
        int endTime = Integer.MIN_VALUE;
        int[] components = new int[dfs2.td.length];
        for (int i = 0; i < dfs2.td.length; ++i) {
            HeapNode node = heap.extractMin();
            NodeData d = (NodeData) node.getData();
            //System.out.printf("ti=%d, tf=%d, idx=%d   endTime=%d c=%d\n", 
            //    node.getKey(), d.tf, d.idx, endTime, c);
            if (node.getKey() <= endTime){
                components[d.idx] = c;
            } else {
                c++;
                endTime = d.tf;
                components[d.idx] = c;
            }
            //System.out.printf("  idx=%d c=%d\n", d.idx, components[d.idx]);
        }
       
        return components;
    }
    
    private static class NodeData {
        int idx;
        int tf;
    }
   
    private SimpleLinkedListNode[] reverse(SimpleLinkedListNode[] connected) {
        
        SimpleLinkedListNode[] reversed = new SimpleLinkedListNode[connected.length];
        for (int i = 0; i < connected.length; ++i) {
            reversed[i] = new SimpleLinkedListNode();
        }
        
        for (int v = 0; v < connected.length; v++) {
            SimpleLinkedListNode uNode = connected[v];
            while (uNode != null && uNode.getNumberOfKeys() > 0) {
                int u = uNode.getKey();
                reversed[u].insert(v);
                uNode = uNode.getNext();
            }
        }
        
        return reversed;
    }
    
}
