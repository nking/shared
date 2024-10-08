package algorithms.shortestPaths;

import java.util.Arrays;

/**
 * from Cormen, Leiserson, Rivest, and Stein "Intro to Algorithms"
 * 
 * finds shortest paths between all pairs of vertices in a directed graph
 * with dynamic programming.
 * 
 * The runtime complexity is <em>O(|V|^3)</em>.
 *
 * Note: consider using DFS single source-all destination search on each
   vertex, making a worst case runtime complexity of O(|V|*(|V| + |E|)).
 * 
 *  
 * W, D
 * 
 * W is an nxn matrix of w[i][j] where each entry is the edge weight if any between i and j.  
 *     has values 0 if i==j and inf where there is no connection.
 * D is an nxn matrix of d[i][j] where each entry is the weight of the shortest path between i and j.
 *  
 *     d[i][j]_k = w[i][j]                                   if k = 0
 *               = min( d[i][k]_(k-1) + d[k][j]_(k-1) )      if k .geq. 1
 *
 * 
 * @author nichole
 */
public class FloydWarshalAllPairs {
    
    /**
     *
     */
    protected int[][] dist = null;
    
    /**
     *
     */
    protected int[][] prev = null;
    
    /**
     *
     */
    protected boolean debug = false;
    
    /**
     *
     */
    public FloydWarshalAllPairs() {
    }
    
    /**
     *
     @param useDebug
     */
    public void setDebug(boolean useDebug) {
        this.debug = useDebug;
    }
    
    /**
     * find the shortest paths between pairs of vertexes.
     @param w a square matrix of w[i][j] where each entry is the edge 
       weight if any between i and j.  It has values 0 if i==j and inf 
       where there is no connection.
     */
    public void findShortestPaths(int[][] w) {
                
        int n = w.length;
        
        dist = new int[n][];
        prev = new int[n][n];
        for (int i = 0; i < n; i++) {
            dist[i] = Arrays.copyOf(w[i], w[i].length);
        }
        
        /*
         * prev[i][j] =    nil    if (i==j) || (w[i][j]==inf)
         *            =    i      if (i!=j) || (w[i][j] < inf)
         */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ( (i == j) || (w[i][j] == Integer.MAX_VALUE) ) {
                    prev[i][j] = Integer.MIN_VALUE;
                } else {
                    prev[i][j] = i;
                }
            }
        }
        
        for (int j = 0; j < n; j++) {
            
            if (debug) {
                System.out.println("k=" + j);
                for (int i = 0; i < n; i++) {
                    System.out.println("dist i=" + i + " : " + Arrays.toString(dist[i]));
                }
                for (int i = 0; i < n; i++) {
                    System.out.println("prev i=" + i + " : " + Arrays.toString(prev[i]));
                }
            }
            
            for (int i = 0; i < n; i++) {
                for (int k = 0; k < n; k++) {
                    int s0 = dist[i][k];
                    int s1 = dist[i][j] + dist[j][k];
                    if (i == k) {
                        dist[i][k] = 0;
                    } else if ((s0 <= s1) || ((dist[i][j] == Integer.MAX_VALUE)
                        || (dist[j][k] == Integer.MAX_VALUE))) {
                        dist[i][k] = s0;
                    } else {
                        dist[i][k] = s1;
                        prev[i][k] = prev[j][k];
                    } 
                }
            }
        }
        
        if (debug) {
            System.out.println("final:");
            for (int i = 0; i < n; i++) {
                System.out.println("dist i=" + i + " : " + Arrays.toString(dist[i]));
            }
            for (int i = 0; i < n; i++) {
                System.out.println("prev i=" + i + " : " + Arrays.toString(prev[i]));
            }
        }
    }
    
    /**
     *
     @param srcIndex
     @param destIndex
     @return
     */
    public int[] getPath(int srcIndex, int destIndex) {
        if (srcIndex > (prev.length - 1)) {
            throw new IllegalArgumentException("srcIndex is out of bounds of prev array");
        }
        if (destIndex > (prev.length - 1)) {
            throw new IllegalArgumentException("destIndex is out of bounds of prev array");
        }
        int[] predecessors = prev[srcIndex];
        
        // assuming no cycles when creating maximum size of return array:
        int[] destnodes = new int[prev.length];
        int pos = destnodes.length - 1;
        
        int idx = destIndex;
        while (pos > -1) {
            if (idx == Integer.MAX_VALUE) {
                return new int[0];
            }
            destnodes[pos] = idx;
            pos--;
            if (idx == srcIndex) {
                break;
            }
            idx = predecessors[idx];
        }
        // move up indexes if pos > 0
        if (pos > 0) {
            int moveUp = pos + 1;

            for (int i = (pos + 1); i < destnodes.length; i++) {
                destnodes[i - moveUp] = destnodes[i];
            }
            int n = destnodes.length - (pos + 1);
            destnodes = Arrays.copyOf(destnodes, n);
        }
        
        return destnodes;
    }
}
