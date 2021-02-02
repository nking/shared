package algorithms.shortestPaths;

import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class AStarTest extends TestCase {
    
    public AStarTest() {
    }
    
    public void testSearch0() throws Exception {
        /*
        from https://en.wikipedia.org/awiki/A*_search_algorithm
        */
        int src = 0; int a=1; int b=2; int c=3; int d=4; int e=5;  int dest=6;
        
        SimpleLinkedListNode[] dAG = new SimpleLinkedListNode[dest + 1];
        TIntIntMap[] weights = new TIntIntHashMap[dAG.length];
        int[] heuristics = new int[dAG.length];
        
        for (int i = 0; i < dAG.length; ++i) {
            dAG[i] = new SimpleLinkedListNode();
            weights[i] = new TIntIntHashMap();
        }
        dAG[src].insert(a);
        weights[src].put(a, (int)(2*1.5));
        
        dAG[src].insert(d);
        weights[src].put(d, (2*2));
        
        dAG[a].insert(b);
        weights[a].put(b, (2*2));
        
        dAG[b].insert(c);
        weights[b].put(c, (2*3));
        
        dAG[c].insert(dest);
        weights[c].put(dest, (2*4));
        
        dAG[d].insert(e);
        weights[d].put(e, (2*3));
        
        dAG[e].insert(dest);
        weights[e].put(dest, (2*2));
        
        heuristics[a] = (2*4);
        heuristics[b] = (2*2);
        heuristics[c] = (2*4);
        heuristics[d] = (int)(2*4.5);
        heuristics[e] = (2*2);
        
        heuristics[src] = Integer.MAX_VALUE;
        heuristics[dest] = 0;
        
        AStar srch = new AStar(dAG, weights, heuristics, src, dest);
        srch.find();
        int[] path = srch.getShortestPathToVertex(dest);
        int dist = srch.getDistanceFromSrc(dest);
        int pathCost = srch.getSumOfPath(path);
        
        int[] expectedPath = new int[]{src, d, e, dest};
        
        assertEquals(4, path.length);
        assertEquals(7*2, dist);
        assertEquals(dist, pathCost);
        assertTrue(Arrays.equals(expectedPath, path));
        /*
         *     5
         *
         *     4 [8]-----------[6]
         *        |             |
         *     3  |       [3]   |
         *        |       /     |
         *     2 [7]   [2]-[4]--.
         *       /    /
         *     1[5] [1]
         *      | /                     [9]
         *-1  [0]/  1   2   3   4   5
         *
         *    -1
         */

        //SimpleLinkedListNode[] dAG, TIntIntMap[] weights, 
        //int[] heuristics, int sourceVertex, int destVertex
        
        
    }
}
