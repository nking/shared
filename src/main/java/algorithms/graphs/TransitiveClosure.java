package algorithms.graphs;

import java.util.Arrays;

/**
 * from Cormen et al. "Intro to Algorithms"
 * 
 * The transitive closure of a graph is the set of pairs of nodes (u, v) such that
   there is a path from u to v of length one or more. 

 * 
 * The runtime complexity is <em>O(V^3)</em>.
 
 <pre>
   TODO: 
      implement Smart TC algorithm:
          “Transitive Closure and Recursive Datalog Implemented on Clusters” by 
          Foto N. Afrati, and Jeffrey D. Ullman,June 17, 2011
      also, consider TC determined from the condensation graph of the 
          Strongly Connected Components (SCC):
          Nuutila 1995, Purdom 1970 in Boost library:
          https://www.boost.org/doc/libs/1_66_0/libs/graph/doc/transitive_closure.html
</pre>
 * @author nichole
 */
public class TransitiveClosure {
    
    //TODO: can conserve space by using a byte array of n bits
    protected boolean[][] t = null;
        
    protected boolean debug = false;
    
    public TransitiveClosure() {
    }
    
    public void setDebug(boolean useDebug) {
        this.debug = useDebug;
    }
    
    /**
     * 
     * @param w
     */
    public void calc(boolean[][] w) {
                
        int n = w.length;
        
        t = new boolean[n][n];
        for (int i = 0; i < n; i++) {
            t[i] = Arrays.copyOf(w[i], w[i].length);
        }
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j || w[i][j]) {
                    t[i][j] = true;
                }
            }
        }
      
        for (int k = 0; k < n; k++) {
            
            if (debug) {
                System.out.println("k=" + k);
                for (int i = 0; i < n; i++) {
                    System.out.println("t i=" + i + " : " + Arrays.toString(t[i]));
                }
            }
            
            boolean s0, s1;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    s0 = t[i][j];
                    s1 = t[i][k] && t[k][j];
                    t[i][j] = s0 | s1;
                }
            }
        }
        
        if (debug) {
            System.out.println("final:");
            for (int i = 0; i < n; i++) {
                System.out.println("t i=" + i + " : " + Arrays.toString(t[i]));
            }
        }
    }
    
}
