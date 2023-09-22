package algorithms.shortestPaths;

import algorithms.misc.MiscMath0;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;
import java.util.ArrayDeque;
import java.util.Queue;

/**
 * given a weighted directed graph with weight function, solves the single
 * source shortest paths.
 * 
 * All edge weights must be non-negative.
 * 
 * The runtime complexity is O(|V| + |E|) 

 * @author nichole
 */
public class BreadthFirstSearch {
    
    /**
     * directed acyclic graph where the main index of the object array
     * is the vertex number, and each vertex has a linked list of outgoing 
     * edges connecting to the vertexes held in each linked list node key.
     */
    protected SimpleLinkedListNode[] g = null;
    
    protected TIntIntMap[] w = null;
    
    /**
     * distances from src to all destination vertexes where index of dist is the vertex index
     */
    protected int[] dist = null;

    /**
     # array in which the indexes are the vertex indexes and the 
       values are the predecessor vertex indexes
     */
    protected int[] predecessor = null;   
    
    /**
     * the src index
     */
    protected final int src;
    
    /**
     *
     @param dAG directed acyclic graph where the main index of the object array
     * is the vertex number, and each vertex has a linked list of outgoing 
     * edges connecting to the vertexes held in each linked list node key.
     * Note that all vertexes, including edge vertexes, must be present as an
     * index of the array dAG, i.e. all vertexes must have numerical value 
     * less than dAG.length.
     @param weights the edge weights in format of outer array index being the
     * first vertex and the map within having a key being the 2nd vertex of the
     * edge with value being the edge weight.
     @param sourceVertex the source vertex index
     */
    public BreadthFirstSearch(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, 
    int sourceVertex) {
        
        if (dAG == null || dAG.length == 0) {
            throw new IllegalArgumentException("dAG cannot be null");
        }
        if (sourceVertex < 0 || sourceVertex >= dAG.length) {
            throw new IllegalArgumentException("sourceIndex cannot be null");
        }

        this.src = sourceVertex;
        
        g = dAG.clone();
        for (int i = 0; i < dAG.length; ++i) {
            g[i] = new SimpleLinkedListNode(dAG[i]);
        }
        w = weights.clone();
        for (int i = 0; i < weights.length; ++i) {
            if (weights[i] != null) {
                w[i] = new TIntIntHashMap(weights[i]);
            }
        }

        dist = new int[g.length];
        predecessor = new int[g.length];
        
        Arrays.fill(dist, Integer.MAX_VALUE);
        Arrays.fill(predecessor, -1);
        
        dist[src] = 0;
        
    }
        
    /**
     *  find the single shortest path in dAG with edge weights w starting from s.
     */
    public void find() {

        int[] visited = new int[g.length];

        Queue<Integer> q = new ArrayDeque<Integer>();

        dist[src] = 0;
        q.add(src);
        visited[src] = 1;
                               
        int u;
        int v;
        SimpleLinkedListNode vNode;
        TIntIntMap uWeights;
        int wUV;
        int dUPlusWUV;

        while (!q.isEmpty()) {

            u = q.poll();
            
            if (visited[u] == 2) {
                continue;
            }

            uWeights = w[u];
            
            vNode = g[u];
            
            while (vNode != null && vNode.getNumberOfKeys() > 0) {
            
                v = vNode.getKey();

                if (!uWeights.containsKey(v)) {
                    throw new IllegalStateException("no weight found for edge " 
                    + u + " to " + v);
                }
                wUV = uWeights.get(v);
            
                dUPlusWUV = wUV;
                if (dist[u] == Integer.MAX_VALUE) {
                    dUPlusWUV = Integer.MAX_VALUE;
                } else {
                    dUPlusWUV += dist[u];
                }

                if (dist[v] > dUPlusWUV) {
                    dist[v] = dUPlusWUV;
                    predecessor[v] = u;
                    if (visited[v] == 0) {
                        visited[v] = 1;
                        q.add(v);
                    }
                }
                
                vNode = vNode.getNext();
            }
            visited[u] = 2;
        }
    }
    /**
     * get shortest path from source to destIndex
     @param destVertex
     @return 
     */
    public int[] getShortestPathToVertex(int destVertex) {
        if (destVertex < 0 || destVertex >= g.length) {
            throw new IllegalArgumentException("destIndex cannot be null");
        }
        if (predecessor == null) {
            throw new IllegalStateException("find must be run before this method can be used");
        }
        
        int[] p = new int[g.length];
        p[p.length - 1] = destVertex;
                
        for (int i = p.length - 2; i > -1; --i) {
            if (destVertex == src) {
                int len = p.length - 1 - i;
                int[] t = new int[len];
                System.arraycopy(p, i + 1, t, 0, len);
                return t;
            } else if (destVertex == -1) {
                throw new IllegalStateException("path did not complete correctly");
            }
            p[i] = predecessor[destVertex];
            destVertex = p[i];
        }
        
        if (p[0] != src) {
            throw new IllegalStateException("path did not complete correctly for destIndex");
        }
        
        return p;
    }
    
    /**
     *
     @param vertexes
     @return
     */
    public int getSumOfPath(int[] vertexes) {
        if (vertexes == null) {
            throw new IllegalArgumentException("vertexes cannot be null");
        }
        if (dist == null) {
            throw new IllegalStateException("find must be run before this method can be used");
        }
        int sum = 0;
        int u, v;
        for (int i = 1; i < vertexes.length; ++i) {
            u = vertexes[i - 1];
            v = vertexes[i];
            
            TIntIntMap uWeights = w[u];
            
            if (!uWeights.containsKey(v)) {
                throw new IllegalStateException("no weight found for edge " 
                    + u + " to " + v);
            }
            
            sum += uWeights.get(v);
        }
        return sum;
    }

}
