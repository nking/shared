package algorithms.shortestPaths;

import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

/**
 *
 * @author nichole
 */
public class BellmanFordTest extends TestCase {
    
    public BellmanFordTest() {
    }
    
    public void test0() {
        
        /*
        from Cormen at al. Fig 24.4
        
             edge     weight
             s t    6   
             s y    7
             t x    5
             t y    8
             t z    -4
             x t    -2
             y x    -3    
             y z    9     
             z s    2
             z x    7
        
        */
        final int s = 0; final int t=1; final int x=2; final int y=3; final int z=4;
        //final int s = 4; final int t=0; final int x=1; final int y=2; final int z=3;
        SimpleLinkedListNode[] g = new SimpleLinkedListNode[5];
        TIntIntMap[] w = new TIntIntMap[5];
         for (int i = 0; i < g.length; ++i) {
             g[i] = new SimpleLinkedListNode();
             w[i] = new TIntIntHashMap();
         }
         g[s].insert(y); g[s].insert(t); 
         w[s].put(t,6); w[s].put(y,7);
         
         g[t].insert(z); g[t].insert(y); g[t].insert(x);  
         w[t].put(x,5); w[t].put(y,8); w[t].put(z,-4);
         
         g[x].insert(t);
         w[x].put(t,-2);
         
         g[y].insert(z); g[y].insert(x); 
         w[y].put(x,-3); w[y].put(z,9);
         
         g[z].insert(x); g[z].insert(s); 
         w[z].put(s,2); w[z].put(x,7);
         
         int src = s;
         
         BellmanFord bf = new BellmanFord(g, w, src);
         boolean noNegativeCycles = bf.find();
         assertTrue(noNegativeCycles);
         
         int[] p;
         int dist, dest;
         
         dest = t;
         p = bf.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{s, y, x, t}, p));
         dist = bf.getSumOfPath(p);
         assertEquals(2, dist);
         
         dest = 2;
         p = bf.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{s, y, x}, p));
         dist = bf.getSumOfPath(p);
         assertEquals(4, dist);
         
         dest = 3;
         p = bf.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{s, y}, p));
         dist = bf.getSumOfPath(p);
         assertEquals(7, dist);
         
         dest = 4;
         p = bf.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{s, y, x, t, z}, p));
         dist = bf.getSumOfPath(p);
         assertEquals(-2, dist);
    }
}
