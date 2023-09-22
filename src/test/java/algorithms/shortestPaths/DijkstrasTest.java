package algorithms.shortestPaths;

import algorithms.util.FormatArray;
import algorithms.util.LinkedListCostNode;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DijkstrasTest extends TestCase {
    
    public DijkstrasTest() {
    }
    
    public void test0() {
        /*
        from Cormen, Leiserson, Rivest, and Stein Fig 24.6
           edges   weight
             s t     10
             s y     5
             t x     1
             t y     2
             x z     4
             y t     3
             y x     9
             y z     2
             z s     7
             z x     4
        
           s  t  x  y  z
           0  1  2  3  4
        */
        SimpleLinkedListNode[] g = new SimpleLinkedListNode[5];
         for (int i = 0; i < g.length; ++i) {
             g[i] = new SimpleLinkedListNode();
         }
         g[0].insert(1); g[0].insert(3);
         g[1].insert(2); g[1].insert(3);
         g[2].insert(4);
         g[3].insert(1); g[3].insert(2); g[3].insert(4);
         g[4].insert(0); g[4].insert(2);
         TIntIntMap[] w = new TIntIntMap[5];
         w[0] = new TIntIntHashMap(); w[0].put(1, 10); w[0].put(3, 5);
         w[1] = new TIntIntHashMap(); w[1].put(2, 1); w[1].put(3, 2);
         w[2] = new TIntIntHashMap(); w[2].put(4, 4);
         w[3] = new TIntIntHashMap(); w[3].put(1, 3); w[3].put(2, 9); w[3].put(4, 2);
         w[4] = new TIntIntHashMap(); w[4].put(0, 7); w[4].put(2, 4);
         
         Dijkstras dijkstras = new Dijkstras(g, w, 0);
         dijkstras.find();
         
         int[] p;
         int dist, dest;
         // s  t  x  y  z
         // 0  1  2  3  4

         dest = 3;
         p = dijkstras.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 3}, p));
         dist = dijkstras.getSumOfPath(p);
         assertEquals(5, dist);
         
         dest = 1;
         p = dijkstras.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 3, 1}, p));
         dist = dijkstras.getSumOfPath(p);
         assertEquals(8, dist);
         
         dest = 2;
         p = dijkstras.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 3, 1, 2}, p)); //  s y t x
         dist = dijkstras.getSumOfPath(p);
         assertEquals(9, dist);

         dest = 4;
         p = dijkstras.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 3, 4}, p));
         dist = dijkstras.getSumOfPath(p);
         assertEquals(7, dist);
         
         dest = 3;
         p = dijkstras.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 3}, p));
         dist = dijkstras.getSumOfPath(p);
         assertEquals(5, dist);

         BreadthFirstSearch bfs = new BreadthFirstSearch(g, w, 0);
         bfs.find();
         //System.out.printf("Dijk.dist=%s\n", FormatArray.toString(dijkstras.dist, "%s"));
         //System.out.printf("bfs.dist=%s\n", FormatArray.toString(bfs.dist, "%s"));
         //System.out.printf("Dijk.pred=%s\n", FormatArray.toString(dijkstras.predecessor, "%s"));
         //System.out.printf("bfs.pred=%s\n", FormatArray.toString(bfs.predecessor, "%s"));

         dest = 2;
         p = bfs.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 3, 1, 2}, p)); //  s y t x
         dist = bfs.getSumOfPath(p);
         assertEquals(9, dist);

         // test BeamSearch
         /*
         edges   weight
             s t     10   0:1
             s y     5    0:3
             t x     1    1:2
             t y     2    1:3
             x z     4    2:4
             y t     3    3:1
             y x     9    3:2
             y z     2    3:4
             z s     7    4:0
             z x     4    4:2

           s  t  x  y  z
           0  1  2  3  4
          */
         LinkedListCostNode[] adjList = new LinkedListCostNode[21];
         for (int i = 0; i < 21; i++) {
              adjList[i] = new LinkedListCostNode();
         }
         adjList[0].insert(1, 10);
         adjList[0].insert(3, 5);
         adjList[1].insert(2, 1);
         adjList[1].insert(3, 2);
         adjList[2].insert(4, 4);
         adjList[3].insert(1, 3);
         adjList[3].insert(2, 9);
         adjList[3].insert(4, 2);
         adjList[4].insert(0, 7);
         adjList[4].insert(2, 4);
         BeamSearch bms = new BeamSearch(adjList, 0, 5);
         bms.search();
         p = bms.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 3, 1, 2}, p)); //  s y t x
         dist = bms.getSumOfPath(p);
         assertEquals(9, dist);
    }
    
}
