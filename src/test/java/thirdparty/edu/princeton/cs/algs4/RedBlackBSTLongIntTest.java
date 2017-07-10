package thirdparty.edu.princeton.cs.algs4;

import algorithms.misc.Misc0;
import gnu.trove.list.TIntList;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.set.hash.TLongHashSet;
import java.security.SecureRandom;
import java.util.Random;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class RedBlackBSTLongIntTest extends TestCase {
    
    public RedBlackBSTLongIntTest(String testName) {
        super(testName);
    }
    
    public void testKeyOperations0() throws Exception {
    
        System.out.println("testKeyOperations");
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.currentTimeMillis();
        //seed = 1499675478087L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        RedBlackBSTLongInt bt = new RedBlackBSTLongInt();
        
        int n = 100;
        
        int count = 0;
        
        TLongList nodes = new TLongArrayList(2*n);
        
        for (int i = 0; i < n/2; ++i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i += skip;
                continue;
            }
            nodes.add(i);
            bt.put(i, i);
            assertTrue(bt.contains(i));
            assertEquals(nodes.get(0), bt.min());
            assertEquals(i, bt.max());
            count++;
        }
        for (int i = (n - 1); i >= (n/2); --i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i -= skip;
                continue;
            }
            nodes.add(i);
            bt.put(i, i);
            assertTrue(bt.contains(i));
            assertEquals(nodes.get(0), bt.min());
            count++;
        }
        nodes.sort();
        int n2 = bt.size();
        assertEquals(count, n2);
        assertEquals(count, nodes.size());
        
        long[] kOutput = new long[2];
        int[] vOutput = new int[2];
        
        for (int nIter = 0; nIter < 3; ++nIter){
            n2 = bt.size();
            assertEquals(n2, nodes.size());

            for (int i = 0; i < n2 - 1; ++i) {
                long idx = nodes.get(i);
                long foundIndex = bt.contains(idx) ? idx : -1;
                assertTrue(foundIndex > -1);

                long expected = nodes.get(i + 1);
                //System.out.println("\n* " + idx + " expected next=" + expected);

                bt.higher(idx, kOutput);
                assertTrue(kOutput[0] != -1);
                long next = kOutput[1];
                //System.out.println(idx + "   next=" + next);
                assertEquals(expected, next);

                if (next > (idx + 1)) {
                    // test ceiling of idx+1
                    bt.ceiling(idx + 1, kOutput);
                    assertTrue(kOutput[0] != -1);
                    long ceil = kOutput[1];
                    assertEquals(expected, ceil);
                }
            }

            for (int i = 1; i < n2; ++i) {
                long idx = nodes.get(i);
                long foundIndex = bt.contains(idx) ? idx : -1;
                assertTrue(foundIndex > -1);

                bt.get(idx, vOutput);
                assertTrue(vOutput[0] != -1);
                assertEquals(idx, vOutput[1]);
                
                long expected = nodes.get(i - 1);
                //System.out.println("\n* " + idx + " expected prev=" + expected);

                bt.lower(idx, kOutput);
                assertTrue(kOutput[0] != -1);
                long prev = kOutput[1];
                //System.out.println(idx + "   prev=" + prev);
                assertEquals(expected, prev);

                if (prev < (idx - 1)) {
                    // test floor of idx-1
                    bt.floor(idx - 1, kOutput);
                    assertTrue(kOutput[0] != -1);
                    long floor = kOutput[1];
                    assertEquals(expected, floor);
                }
            }
        
            if ((nIter & 1) == 1) {
                //randomly remove some nodes
                for (int i = 0; i < n2/4; ++i) {
                    int idx = rand.nextInt(nodes.size());
                    long v = nodes.get(idx);
                    assertTrue(bt.contains(v));
                    bt.delete(v);
                    assertFalse(bt.contains(v));
                    nodes.removeAt(idx);
                    assertEquals(nodes.size(), bt.size());
                }
            } else {
                //randomly add some nodes
                for (int i = 0; i < n2/4; ++i) {
                    int idx = nodes.size() + rand.nextInt(2*nodes.size());
                    // the nodes contains is linear search, so could use a temp
                    // set if this test gets large one day
                    if (!nodes.contains(idx)) {
                        bt.put(idx, idx);
                        nodes.add(idx);
                        assertTrue(bt.contains(idx));
                    }
                }
            }
            
            n2 = nodes.size();
            nodes.sort();
           
            assertEquals(nodes.get(nodes.size()/2), bt.select(nodes.size()/2));            
                    
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
            bt.deleteMin();
            assertFalse(bt.contains(min));
            nodes.removeAt(0);
        }        
    }
    
}
