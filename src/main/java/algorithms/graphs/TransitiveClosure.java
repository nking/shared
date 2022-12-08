package algorithms.graphs;

import algorithms.VeryLongBitString;
import algorithms.misc.MiscMath0;
import java.util.Arrays;

/**
 * from Cormen, Leiserson, Rivest, and Stein "Intro to Algorithms"
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
    
    /**
     *
     */
    protected boolean debug = false;
    
    /**
     *
     */
    public TransitiveClosure() {
    }
    
    /**
     *
     @param useDebug
     */
    public void setDebug(boolean useDebug) {
        this.debug = useDebug;
    }
    
    /**
     * runtime complexity is O(n^3) where n is w.length
     @param w a square adjacency matrix for a DAG with |V|=w.length.
     @return 
     */
    public boolean[][] calc(boolean[][] w) {
                
        int n = w.length;
        
        boolean[][] t = new boolean[n][n];
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
            
            if (debug) {
                for (int i = 0; i < n; i++) {
                    System.out.println("    t i=" + i + " : " + Arrays.toString(t[i]));
                }
            }
        }
        
        return t;
    }
    
    /**
     *
     @param w
     @return
     */
    public static VeryLongBitString[] convert(boolean[][] w) {
        int n = w.length;

        int nBits = MiscMath0.numberOfBitsWOB(n);
        VeryLongBitString[] bs = new VeryLongBitString[n];
        for (int i = 0; i < n; ++i) {
            bs[i] = new VeryLongBitString(nBits);
        }
        for (int i = 0; i < w.length; ++i) {
            for (int j = 0; j < w[i].length; ++j) {
                if (w[i][j]) {
                    bs[i].setBit(j);
                }
            }
        }
        return bs;
    }
    
    /**
     * runtime complexity is O(n^3) where n is w.length.
     @param w a square adjacency matrix for a DAG with |V|=w.length.
     @return 
     */
    public VeryLongBitString[] calc(VeryLongBitString[] w) {
                
        int n = w.length;
        
        VeryLongBitString[] t = new VeryLongBitString[n];
        for (int i = 0; i < n; i++) {
            t[i] = w[i].copy();
        }
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j || w[i].isSet(j)) {
                    t[i].setBit(j);
                }
            }
        }
      
        for (int k = 0; k < n; k++) {
            
            if (debug) {
                System.out.println("k=" + k);
                for (int i = 0; i < n; i++) {
                    System.out.println("t i=" + i + " : " + t[i].toString());
                }
            }
            
            boolean s0, s1;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    s0 = t[i].isSet(j);//t[i][j];
                    s1 = t[i].isSet(k) && t[k].isSet(j);//t[i][k] && t[k][j];
                    if (s0 | s1) {
                        t[i].setBit(j);
                    } else {
                        t[i].clearBit(j);
                    }
                }
            }
            
            if (debug) {
                for (int i = 0; i < n; i++) {
                    System.out.println("    t i=" + i + " : " + t[i].toString());
                }
            }
        }
        
        return t;
    }
}
