package algorithms.shortestPaths;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.Arrays;

/**
 <pre>

 * solves the single source to all destinations shortest paths problem.
 * 
 * It can handle <em>negative edge weights</em>, but if there is a negative cycle, it returns false.
 * It solves the problem using dynamic programming (= dividing the problem into
 * sub-problems and combining them into the solution).
 *
 * It is non-greedy solution that uses the relax function to find the minimum
 * edge for a node.  relax is invoked |V| - 1 times.
 * 
 * Runtime complexity is <em>O(|V||E|)</em>.
 * 
 * implemented following Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms"
 </pre>

 * @author nichole
 */
public class BellmanFord {
    
    private int nV;
    
    /**
     *
     */
    protected int[] dist = null;

    /**
     *
     */
    protected int[] predecessor = null;   
    
    /**
     *
     */
    protected int src = -1;
    
    private static int sentinel = Integer.MAX_VALUE;
    
    private TObjectIntMap<PairInt> edges;
    
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
    public BellmanFord(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, int sourceVertex) {
    
        if (dAG == null || dAG.length == 0) {
            throw new IllegalArgumentException("dAG cannot be null");
        }
        if (sourceVertex < 0 || sourceVertex >= dAG.length) {
            throw new IllegalArgumentException("sourceIndex cannot be null");
        }
        
        init(dAG, weights, sourceVertex);
        
    }
        
    /**
     * find the single shortest paths from start s to all destinations in dAG with edge weights w.
     @return returns false if a negative cycle is present, else returns true 
     * and the results are usable.
     */
    public boolean find() {
        /*
        init
        for i=1 to |V[G]| - 1
           for each edge (u,v) in E[G]
              relax(u,v,w)
        for each edge (u,v) in E[G]
           if d(v) > d(u) + w(u,v)
              return false
        return true
        */
        
        PairInt uv;
        int wUV, u, v;
        TObjectIntIterator<PairInt> iter;
        
        for (int ii = 0; ii < nV; ++ii) {
            iter = edges.iterator();
            for (int i = 0; i < edges.size(); ++i) {
                iter.advance();
                uv = iter.key();
                u = uv.getX();
                v = uv.getY();
                wUV = iter.value();
                relax(u, v, wUV);
            }
            
            /*System.out.printf("    d=[%d, %d, %d, %d, %d]\n   p=[%d, %d, %d, %d, %d]\n",
                    dist[0], dist[1], dist[2], dist[3], dist[4],
                    predecessor[0], predecessor[1], predecessor[2], predecessor[3], predecessor[4]);
            System.out.flush();*/
        }
        
        return checkForNegativeCycle();
    }
    
    private void relax(int u, int v, int wUV) {
    
        /*
        if dv > du + wuv
            dv = du + wuv
            pred[v]=u
        */
        
        int dUPlusWUV = wUV;
        if (dist[u] == sentinel) {
            dUPlusWUV = sentinel;
        } else {
            dUPlusWUV += dist[u];
        }

        if (dist[v] > dUPlusWUV) {
            dist[v] = dUPlusWUV;
            predecessor[v] = u;
        }
        
    }

    private void init(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, int sourceVertex) {
        
        nV = dAG.length;
        
        src = sourceVertex;
            
        dist = new int[nV];
        predecessor = new int[nV];
        
        Arrays.fill(dist, sentinel);
        Arrays.fill(predecessor, -1);
        
        dist[src] = 0;
               
        edges = new TObjectIntHashMap<PairInt>();
        for (int u = 0; u < nV; ++u) {
            SimpleLinkedListNode vNode = dAG[u];            
            TIntIntMap uWgts = weights[u];
            while (vNode != null && vNode.getNumberOfKeys() > 0 && uWgts != null) {
                int v = vNode.getKey();
                if (uWgts.containsKey(u)) {
                    throw new IllegalStateException("no weight found for edge "
                        + u + " to " + v);
                }
                int wUV = uWgts.get(v);
                edges.put(new PairInt(u, v), wUV);
                vNode = vNode.getNext();
            }
        }
    }
    
    /**
     * get shortest path from source to destIndex
     @param destVertex
     @return 
     */
    public int[] getShortestPathToVertex(int destVertex) {
        if (destVertex < 0 || destVertex >= nV) {
            throw new IllegalArgumentException("destIndex cannot be null");
        }
        if (predecessor == null) {
            throw new IllegalStateException("find must be run before this method can be used");
        }
        
        int[] p = new int[nV];
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
        PairInt uv;
        for (int i = 1; i < vertexes.length; ++i) {
            u = vertexes[i - 1];
            v = vertexes[i];
            uv = new PairInt(u, v);
            if (!edges.containsKey(uv)) {
                throw new IllegalStateException("no weight found for edge " 
                    + u + " to " + v);
            }
                        
            sum += edges.get(uv);
        }
        return sum;
    }

    private boolean checkForNegativeCycle() {
        
        TObjectIntIterator<PairInt> iter = edges.iterator();
        PairInt uv;
        int wUV, u, v;
        for (int i = 0; i < edges.size(); ++i) {
            iter.advance();
            uv = iter.key();
            u = uv.getX();
            v = uv.getY();
            wUV = iter.value();

            if (dist[u] == sentinel) {
                throw new IllegalStateException("dist[" + u + "] was not determined");
            }

            int dUPlusWUV = wUV;
            if (dist[u] == sentinel) {
                dUPlusWUV = sentinel;
            } else {
                dUPlusWUV += dist[u];
            }

            if (dist[v] > dUPlusWUV) {
                return false;
            }
        }
        
        return true;
    }
}
