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
public class RedBlackBSTLongInt2Test extends TestCase {
    
    public RedBlackBSTLongInt2Test(String testName) {
        super(testName);
    }
    
    public void testKeyOperations00() throws Exception {
    
        System.out.println("testKeyOperations");
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.currentTimeMillis();
        seed = 1499931908167L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        RedBlackBSTLongInt2 bt = new RedBlackBSTLongInt2();
        
        int n = 100;
        
        int count = 0;
        
        TLongList nodes = new TLongArrayList(2*n);
        
        long[] kOutput = new long[2];
        
        for (int i = 0; i < n/2; ++i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i += skip;
                continue;
            }
            nodes.add(i);
            bt.put(i, i);
            assertTrue(bt.rootIsSet);
            assertTrue(bt.contains(i));
            bt.min(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(nodes.get(0), kOutput[1]);
            
            bt.max(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(i, kOutput[1]);
            
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
            assertTrue(bt.rootIsSet);
            assertTrue(bt.contains(i));
            bt.min(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(nodes.get(0), kOutput[1]);
            count++;
        }
        nodes.sort();
        int n2 = bt.size();
        assertEquals(count, n2);
        assertEquals(count, nodes.size());
        
        
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

                assertTrue(bt.rootIsSet);
                
 //System.out.println("nIter=" + nIter + " i=" + i + " idx=" + idx);
 
                if (next > (idx + 1)) {
                    // test ceiling of idx+1
                    bt.ceiling(idx + 1, kOutput);
                    assertTrue(kOutput[0] != -1);
                    long ceil = kOutput[1];
                    assertEquals(expected, ceil);
                }
            }
            
            bt.higher(nodes.get(nodes.size() - 1), kOutput);
            assertTrue(kOutput[0] == -1);
            bt.lower(nodes.get(0), kOutput);
            assertTrue(kOutput[0] == -1);
            
            for (int i = 1; i < n2; ++i) {
                long idx = nodes.get(i);
                long foundIndex = bt.contains(idx) ? idx : -1;
                assertTrue(foundIndex > -1);

                bt.get(idx, vOutput);
                assertTrue(vOutput[0] != -1);
                assertEquals(idx, vOutput[1]);
                
                assertTrue(bt.rootIsSet);
                
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
           
            n2 = nodes.size();
            nodes.sort();
           
            assertEquals(nodes.get(nodes.size()/2), bt.select(nodes.size()/2));            
                  
            assertEquals(nodes.size(), bt.size());
            
            System.out.println("before deleteMin bt.size=" + bt.size());
            
            
            long max = nodes.get(nodes.size() - 1);
            assertTrue(bt.contains(max));
            assertTrue(bt.rootIsSet);
            System.out.println("will delete max=" + max + " from this tree:");
            //bt.printPreOrderTraversal();
            bt.deleteMax();
            System.out.println("after deleteMax bt.size=" + bt.size());
            assertTrue(bt.rootIsSet);
            //System.out.println(" after delete:");
            //bt.printPreOrderTraversal();
            assertFalse(bt.contains(max));
            nodes.removeAt(nodes.size() - 1);
            assertEquals(nodes.size(), bt.size());
            
            long min = nodes.get(0);
            System.out.println("deleting min=" + min + 
                "  bt.rootIsSet=" + bt.rootIsSet + " bt.size=" + bt.size());
            assertTrue(bt.rootIsSet);
            boolean contains = bt.contains(min);
            assertTrue(contains);
            bt.deleteMin();
            System.out.println("after deleteMin bt.size=" + bt.size());
            assertFalse(bt.contains(min));
            nodes.removeAt(0);
            assertEquals(nodes.size(), bt.size());
        }            
    }
    
    public void estKeyOperations0() throws Exception {
    
        System.out.println("testKeyOperations");
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.currentTimeMillis();
        seed = 1499931908167L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        RedBlackBSTLongInt2 bt = new RedBlackBSTLongInt2();
        
        int n = 100;
        
        int count = 0;
        
        TLongList nodes = new TLongArrayList(2*n);
        
        long[] kOutput = new long[2];
        
        for (int i = 0; i < n/2; ++i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i += skip;
                continue;
            }
            nodes.add(i);
            bt.put(i, i);
            assertTrue(bt.rootIsSet);
            assertTrue(bt.contains(i));
            bt.min(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(nodes.get(0), kOutput[1]);
            
            bt.max(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(i, kOutput[1]);
            
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
            assertTrue(bt.rootIsSet);
            assertTrue(bt.contains(i));
            bt.min(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(nodes.get(0), kOutput[1]);
            count++;
        }
        nodes.sort();
        int n2 = bt.size();
        assertEquals(count, n2);
        assertEquals(count, nodes.size());
        
        
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

                assertTrue(bt.rootIsSet);
                
 System.out.println("nIter=" + nIter + " i=" + i + " idx=" + idx);
 
                if (next > (idx + 1)) {
                    // test ceiling of idx+1
                    bt.ceiling(idx + 1, kOutput);
                    assertTrue(kOutput[0] != -1);
                    long ceil = kOutput[1];
                    assertEquals(expected, ceil);
                }
            }
            
            bt.higher(nodes.get(nodes.size() - 1), kOutput);
            assertTrue(kOutput[0] == -1);
            bt.lower(nodes.get(0), kOutput);
            assertTrue(kOutput[0] == -1);
            
            for (int i = 1; i < n2; ++i) {
                long idx = nodes.get(i);
                long foundIndex = bt.contains(idx) ? idx : -1;
                assertTrue(foundIndex > -1);

                bt.get(idx, vOutput);
                assertTrue(vOutput[0] != -1);
                assertEquals(idx, vOutput[1]);
                
                assertTrue(bt.rootIsSet);
                
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
                    assertTrue(bt.rootIsSet);
                    int idx = rand.nextInt(nodes.size());
                    long v = nodes.get(idx);
                    assertTrue(bt.contains(v));
                    
                    //bt.printPreOrderTraversal();
                    System.out.println("delete " + v + " idx=" + idx);
                    
                    bt.delete(v);
                    assertTrue(bt.rootIsSet);
                    assertFalse(bt.contains(v));
                    nodes.removeAt(idx);
                    assertEquals(nodes.size(), bt.size());
                    
         
        //DEBUGGING one method at a time
        if (true) {
            return;
        } 
        
                }
        
        //DEBUGGING one method at a time
        if (true) {
            return;
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
                        assertTrue(bt.rootIsSet);
                        assertTrue(bt.contains(idx));
                    }
                }
            }
            
            n2 = nodes.size();
            nodes.sort();
           
            assertEquals(nodes.get(nodes.size()/2), bt.select(nodes.size()/2));            
                  
            assertEquals(nodes.size(), bt.size());
            
            System.out.println("before deleteMin bt.size=" + bt.size());
            
            long max = nodes.get(nodes.size() - 1);
            assertTrue(bt.contains(max));
            assertTrue(bt.rootIsSet);
            System.out.println("will delete max=" + max + " from this tree:");
            //bt.printPreOrderTraversal();
            bt.deleteMax();
            System.out.println("after deleteMax bt.size=" + bt.size());
            assertTrue(bt.rootIsSet);
            //System.out.println(" after delete:");
            //bt.printPreOrderTraversal();
            assertFalse(bt.contains(max));
            nodes.removeAt(nodes.size() - 1);
            assertEquals(nodes.size(), bt.size());

            long min = nodes.get(0);
            System.out.println("deleting min=" + min + 
                "  bt.rootIsSet=" + bt.rootIsSet + " bt.size=" + bt.size());
            assertTrue(bt.rootIsSet);
            boolean contains = bt.contains(min);
            assertTrue(contains);
            bt.deleteMin();
            System.out.println("after deleteMin bt.size=" + bt.size());
            assertFalse(bt.contains(min));
            nodes.removeAt(0);
            assertEquals(nodes.size(), bt.size());
        }            
    }
    
}
