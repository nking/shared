package algorithms.msts;

import algorithms.bipartite.MinHeapForRT2012;
import algorithms.heapsAndPQs.HeapNode;
import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.Arrays;

/**
 * minimum spanning tree is the subset of edges in a weighted undirected graph
 * that connect all vertexes for a total minimum cost (sum of edge weights).
 * 
 * Prim's is the same as Dijkstra's excepting 2 things:
    (1) d[v] is the minimum cost of any edge connecting to v
    and (2) the relax function compares the weight of u,v to d[v]
 * 
 * Implemented from pseudo code in Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms and
      from http://en.wikipedia.org/wiki/Prim's_algorithm
     
      Time complexity for different implementations:
     
          Minimum edge weight data structure   Time complexity (total)
          ----------------------------------   -----------------------
          adjacency matrix, searching          O(V^2)
          binary heap and adjacency list       O((V + E) lg2 V) = O(E lg2 V)
          Fibonacci heap and adjacency list    O(E + V lg2 V) 
          YFastTrie and adjacency list         O(E + V)
     
      Prim's algorithm:
     
      Grow a Tree in Running time is O((|N| + |E|)log|N|)
      -- Start by picking any vertex to be the root of the tree.
      -- While the tree does not contain all vertices in the graph find shortest
         edge leaving the tree and add it to the tree.
   
  * a few definitions:
 *     a cut(S, V-S) of an undirected graph G=(V,E) is a partition of V.
 *     
 *     an edge (u,v) in E crosses the cut (S, V-S) if one of its end points is 
 *        in S and the other is in V-S
 *        
 *     a cut respects a set of edges A if no edge in A crosses the cut
 *     
 *     an edge is a light edge crossing a cut if its weight is the minimum of any edge
 *         crossing the cut.   a tie can mean more than one light edge for a cut.
 *         
 *    
 *     the goal is to visit every node in the input graph in a greedy BFS style (that is
 *     find the next best connected edge and continue from there) adding an edge if
 *     the end points are not already in the tree and if the edge is the minimum of the
 *     u neighbors.
 *     
 *     The input graph and edges are determined by cost, so if the edges are undirected,
 *     be sure to create cost[u][v] = value and cost[v][u] = value.

* 
  Following pseudo-code from Introduction to Algorithms,
  by Cormen, Leiserson, Rivest, and Stein

* 
 * this implementation uses a YFastTrie min priority queue and adjacency list.
 * 
* 
 * @author nichole
 */
public class PrimsMinimumSpanningTree<T> {
    
    /**
     * find a minimum spanning tree using Prim's algorithm.
     * @param graph an adjacency list for the graph, where the index of the array
     * is the vertex number and the keys within each vertex's list node are the
     * vertex on the other end of an edge.
     * @param edgeWeights map of edge weights with key being pairs of the
     * vertex numbers.
     * @param r root node of final mst
     * @return a minimum spanning tree of the weighted given graph
     */
    public static TIntObjectMap<SimpleLinkedListNode> mst(
        SimpleLinkedListNode[] graph, TObjectIntMap<PairInt> edgeWeights,
        int r) {
        
        /* MST-Prim(G, w, r):
         * 
         * for each u in V[G]
         *     do d[u] = inf
         *         prev[u] = nil
         * d[r] = 0
         * Q = V[G]
         * while (Q != 0)
         *     do u = extractMin(Q)
         *         for each v in Adj[u]
         *             do if v is in Q and w(u,v) < d[v]
         *                 then prev[v] = u
         *                     d[v] = w(u,v)
         */
        
        //max(weights) is needed to estimate maximum number of bits needed by trie
        int maxValue = findMax(edgeWeights);
        int sentinel = maxValue;
        if (maxValue < Integer.MAX_VALUE) {
            sentinel = maxValue + 1;
        }
        int maxNBits = (int)Math.ceil(Math.log(sentinel/Math.log(2)));
         
        int nE = edgeWeights.size();
        int nV = graph.length;
        
        int[] d = new int[nV];
        Arrays.fill(d, sentinel);
        
        int[] prev = new int[nV];
        Arrays.fill(prev, -1);
        
        HeapNode[] nodes = new HeapNode[nV];
        
        d[r] = 0;
        
        HeapNode node;
        PairInt uv;
        
        //O(|V|)
        //int maxValue, int approxN, int maxNumberOfBits
        MinHeapForRT2012 heap = new MinHeapForRT2012(sentinel, nV, maxNBits);
        for (int v = 0; v < nV; ++v) {
            node = new HeapNode(d[v]);
            node.setData(Integer.valueOf(v));
            heap.insert(node);
            nodes[v] = node;
        }
        
        HeapNode uNode;
        int u;
        SimpleLinkedListNode vNode;
        int v;
        int wUV;
        
        // worse case O((|V| + |E|)*(lg_2 lg_2 (maxNBits)))
        while (heap.getNumberOfNodes() > 0) {
            
            //essentially O(small constant of lg_2 lg_2 (maxNBits))
            uNode = heap.extractMin();
            
            u = (Integer)uNode.getData();
            
            // null the entry in nodes so it isn't used in decrease key
            nodes[u] = null;
            
            vNode = graph[u];
                        
            while (vNode != null && vNode.getNumberOfKeys() > 0) {
            
                v = vNode.getKey();
                
                uv = new PairInt(u, v);
                
                if (!edgeWeights.containsKey(uv)) {
                    throw new IllegalStateException("no weight found for edge " 
                    + u + " to " + v);
                }
                if (nodes[v] == null) {
                    // v is no longer in the heap
                    vNode = vNode.getNext();                   
                    continue;
                }
                
                wUV = edgeWeights.get(uv);
                
                //compare d[v] to only wUV instead of dijkstra's (d[u] + wUV)             
                if (wUV < d[v]) {
                    prev[v] = u;
                    d[v] = wUV;
                    //essentially O(small constant of lg_2 lg_2 (maxNBits))
                    heap.decreaseKey(nodes[v], wUV);
                }
                
                vNode = vNode.getNext();
            }
        }
        
        // read predecessors to add edges to tree 'a'
        
        int nMSTEdges = 0;
        long sum = 0;
        TIntObjectMap<SimpleLinkedListNode> a = new TIntObjectHashMap<SimpleLinkedListNode>();
        for (v = 0; v < nV; ++v) {
            if (prev[v] > -1) {
                if (a.containsKey(prev[v])) {
                    a.get(prev[v]).insert(v);
                } else {
                    a.put(prev[v], new SimpleLinkedListNode(v));
                }
                //System.out.printf("adding %d to %d\n", pi[v], v);
                nMSTEdges++;
                PairInt uv2 = new PairInt(prev[v], v);
                sum += edgeWeights.get(uv2);
            }
        }
        
        System.out.printf("%d edges out of %d in minimumspanning tree.  sum=%d\n",
            nMSTEdges, nE, sum);
        
        return a;
    }

    private static int findMax(TObjectIntMap<PairInt> edgeWeights) {
        
        int max = Integer.MIN_VALUE;
        
        int w;
        TObjectIntIterator<PairInt> iter = edgeWeights.iterator();
        for (int i = 0; i < edgeWeights.size(); ++i) {
            iter.advance();
            w = iter.value();
            if (w > max) {
                max = w;
            }
        }
        
        return max;
    }
}
