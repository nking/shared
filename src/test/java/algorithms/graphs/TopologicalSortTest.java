package algorithms.graphs;

import algorithms.util.SimpleLinkedListNode;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TopologicalSortTest extends TestCase {

    public void testSort() {
        
        /*
        from Cormen, Leiserson, Rivest, and Stein Fig 22.7
        
        11/16 undershorts     socks 17/18    
               |         \     |               watch 9/10
               V          V    V
          12/15 pants ----> shoes 13/14
               |
               |     shirt  1/8
               V     /   V
          6/7 belt <    tie  2/5
                  \      |
                   \     V
                    > jacket 3/4
        
        ordered by tf (reversed post-order traversal):
             r        r                       r       r
        // socks, undershorts, pants, shoes, watch, shirt, belt, tie, jacket
        tf  18        16        15      14    10      8      7    5     4
        
        ordered by ti (pre-order traversal):
             1      2     3     6      9        11        12      13    17
          shirt, tie,  jacket, belt, watch, undershorts, pants, shoes, socks
            r                          r        r                        r
        */
        Map<String, Integer> objs = new HashMap<String, Integer>();
        objs.put("shirt", 0);      // belt, tie
        objs.put("belt", 1);       // jacket
        objs.put("tie", 2);        // jacket
        objs.put("jacket", 3);     //
        objs.put("watch", 4);      //
        objs.put("undershorts", 5);//  pants, shoes
        objs.put("pants", 6);      //  shoes
        objs.put("shoes", 7);      //
        objs.put("socks", 8);      //  shoes
        
        SimpleLinkedListNode[] connected = new SimpleLinkedListNode[9];
        for (int i = 0; i < connected.length; ++i) {
            connected[i] = new SimpleLinkedListNode();
        }
        connected[objs.get("shirt")].insert(objs.get("belt")); connected[objs.get("shirt")].insert(objs.get("tie"));
        connected[objs.get("belt")].insert(objs.get("jacket"));
        connected[objs.get("tie")].insert(objs.get("jacket"));
        connected[objs.get("undershorts")].insert(objs.get("pants")); connected[objs.get("undershorts")].insert(objs.get("shoes"));
        connected[objs.get("pants")].insert(objs.get("shoes")); connected[objs.get("pants")].insert(objs.get("shoes"));
        connected[objs.get("socks")].insert(objs.get("shoes"));
        // socks, undershorts, pants, shoes, watch, shirt, belt, tie, jacket
        int[] expResult = new int[]{objs.get("socks"), objs.get("undershorts"),
            objs.get("pants"), objs.get("shoes"), objs.get("watch"),
            objs.get("shirt"), objs.get("belt"), objs.get("tie"),
            objs.get("jacket")
        };

        TopologicalSort ts = new TopologicalSort(connected);
        
        int[] result = ts.sort();
        
        //System.out.println("expected = " + Arrays.toString(expResult));
        //System.out.println("result = " + Arrays.toString(result));
        
        assertTrue(Arrays.equals(expResult, result));

        int[] result2 = ts.sortKahn();
        // another valid ordering
        int[] expResult2 = new int[]{0, 4, 5, 8, 2, 1, 6, 3, 7};
        assertTrue(Arrays.equals(expResult2, result2));
    }
    
    /**
     * Test of sort method, of class TopologicalSort.
     */
    public void testSort_SimpleDAG() {

        System.out.println("testSort_SimpleDAG");
        
        // constructing tests from MIT open courseware
        // network_optimization/MIT15_082JF10_av03.pdf
        
        SimpleLinkedListNode[] connected = new SimpleLinkedListNode[8];

        /*              <5> --------> <0>
         *            >  |       >     |
         *         .     V   .         V
         *      <4> --> <1>           <7> ---> <2>
         *      >  .                 >        .>
         *     .     .            .      .
         *    .         >      .   .
         *  <6> ------> <3> .
         */
        /*
        expected=[6, 4, 3, 5, 1, 0, 7, 2]
        result=  [6, 4, 5, 3, 1, 0, 7, 2]
        6 4 5
        3
        1
        0 7 2
        */
        
        for (int i = 0; i < 8; ++i) {
            connected[i] = new SimpleLinkedListNode();
        }
        
        connected[0].insert(7);

        connected[1].insert(0);

        connected[3].insert(2);
        connected[3].insert(7);

        connected[4].insert(1);
        connected[4].insert(5);
        connected[4].insert(3);

        connected[5].insert(0);
        connected[5].insert(1);

        connected[6].insert(3);
        connected[6].insert(4);

        connected[7].insert(2);

        int[] expResult = new int[]{6, 4, 3, 5, 1, 0, 7, 2};

        TopologicalSort ts = new TopologicalSort(connected);
        
        int[] result = ts.sort();
        int[] result2 = ts.sortKahn();
        
        /*
        expected=[6, 4, 3, 5, 1, 0, 7, 2]
        result=  [6, 4, 5, 3, 1, 0, 7, 2]
        6 4 5
        3
        1
        0 7 2
        */
        //System.out.println("expected=" + Arrays.toString(expResult));
        //System.out.println("result=  " + Arrays.toString(result));

        //assertTrue(Arrays.equals(expResult, result));
        assertTrue(Arrays.equals(expResult, result2));
    }
   
    public void testSort2() {

        System.out.println("testSort2");
        
        // constructing test from Cormen, Leiserson, Rivest, and Stein's "Introduction to Algorithms"
        SimpleLinkedListNode[] connected = new SimpleLinkedListNode[9];

        /*    *0 \         *3
         *    ||    \      ||        *8
         *    \/       \-> \/
         *    *1 --------> *4
         *    ||
         *    ||    /*5
         *    ||   / ||
         *    \/<-/  ||
         *    *2     \/
         *      \    *6
         *       \   ||
         *       \   ||
         *        \->\/
         *           *7
         */
        for (int i = 0; i < connected.length; i++) {
            connected[i] = new SimpleLinkedListNode();
        }
        connected[0].insert(1);
        connected[0].insert(4);
        connected[1].insert(2);
        connected[1].insert(4);
        connected[2].insert(7);
        connected[3].insert(4);
        connected[5].insert(2);
        connected[5].insert(6);
        connected[6].insert(7);
        
        /*    *0 \         *3
         *    ||    \      ||        *8
         *    \/       \-> \/
         *    *1 --------> *4
         *    ||
         *    ||    /*5
         *    ||   / ||
         *    \/<-/  ||
         *    *2     \/
         *      \    *6
         *       \   ||
         *       \   ||
         *        \->\/
         *           *7

        Book solution:
        3   0   1   4    8    5    2    6    7
        ------------>         --------->----->
            --->---->         ---->---------->
            -------->
                8, 5, 6, 3, 0, 1, 2, 7, 4
        result=  [8, 
                  5, 
                  6, 
                  3, 0, 1, 4, 
                  2, 
                  7]
         */

        int[] expResult = new int[]{3, 0, 1, 4, 8, 5, 2, 6, 7};

        TopologicalSort ts = new TopologicalSort(connected);
        
        int[] result = ts.sort();
        
        //System.out.println("expected=" + Arrays.toString(expResult));
        System.out.println("result=  " + Arrays.toString(result));

        int[] result2 = ts.sortKahn();
        assertNotNull(result2);

        //assertTrue(Arrays.equals(expResult, result));
    }
}
