package thirdparty.edu.princeton.cs.algs4;

import algorithms.misc.Misc0;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TLongArrayList;
import java.util.Random;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class RedBlackBSTTest extends TestCase {
    
    public RedBlackBSTTest(String testName) {
        super(testName);
    }
    
    public void testKeyOperations0() throws Exception {
    
        System.out.println("testKeyOperations");
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.currentTimeMillis();
        //seed = 1500070815033L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        RedBlackBST<Long, Integer> bt = new RedBlackBST<Long, Integer>();
        
        int n = 100;
        
        int count = 0;
        
        TLongList nodes = new TLongArrayList(2*n);
        
        for (long i = 0; i < n/2; ++i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i += skip;
                continue;
            }
            nodes.add(i);
            bt.put(i, (int)i);
            assertTrue(bt.contains(i));
            assertEquals(nodes.get(0), bt.min().longValue());
            assertEquals(i, bt.max().longValue());
            assertEquals(nodes.size(), bt.size());
            count++;
        }
        for (long i = (n - 1); i >= (n/2); --i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i -= skip;
                continue;
            }
            nodes.add(i);
            bt.put(i, (int)i);
            assertTrue(bt.contains(i));
            assertEquals(nodes.get(0), bt.min().longValue());
            assertEquals(nodes.size(), bt.size());
            count++;
        }
        nodes.sort();
        int n2 = bt.size();
        assertEquals(count, n2);
        assertEquals(count, nodes.size());
        
        for (int nIter = 0; nIter < 3; ++nIter){
            n2 = bt.size();
            assertEquals(n2, nodes.size());

            for (int i = 0; i < n2 - 1; ++i) {
                long idx = nodes.get(i);
                long foundIndex = bt.contains(idx) ? idx : -1;
                assertTrue(foundIndex > -1);

                assertEquals(nodes.size(), bt.size());
                //System.out.println("\n* " + idx + " expected next=" + expected);

            }
            
            
            if ((nIter & 1) == 1) {
                //randomly remove some nodes
                for (int i = 0; i < n2/4; ++i) {
                    int idx = rand.nextInt(nodes.size());
                    long v = nodes.get(idx);
                    assertTrue(bt.contains(v));
                    
                    //bt.printPreOrderTraversal();
                    //System.out.println("delete " + v + " idx=" + idx);
                    
                    bt.delete(v);
                    assertFalse(bt.contains(v));
                    nodes.removeAt(idx);
                    assertEquals(nodes.size(), bt.size());
                    
                    //bt.printPreOrderTraversal();
                    //System.out.println("AFTER delete " + v + " idx=" + idx);
                }
                
            } else {
                //randomly add some nodes
                for (int i = 0; i < n2/4; ++i) {
                    int idx = nodes.size() + rand.nextInt(2*nodes.size());
                    // the nodes contains is linear search, so could use a temp
                    // set if this test gets large one day
                    if (!nodes.contains(idx)) {
                        bt.put((long)idx, idx);
                        nodes.add(idx);
                        assertTrue(bt.contains((long)idx));
                    }
                }
            }
            
            n2 = nodes.size();
            nodes.sort();
           
            assertEquals(nodes.get(nodes.size()/2), 
                bt.select(nodes.size()/2).longValue());            
                    
            long max = nodes.get(nodes.size() - 1);
            assertTrue(bt.contains(max));
            //System.out.println("will delete max=" + max + " from this tree:");
            //bt.printPreOrderTraversal();
            bt.deleteMax();
            //System.out.println(" after delete:");
            //bt.printPreOrderTraversal();
            assertFalse(bt.contains(max));
            nodes.removeAt(nodes.size() - 1);

            long min = nodes.get(0);
            assertTrue(bt.contains(min));
            //System.out.println("will delete min=" + min);
            bt.deleteMin();
            assertFalse(bt.contains(min));
            nodes.removeAt(0);
        }            
    }
    
}
