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
        }
        
        for (int i = 1; i < n2; ++i) {
            long idx = nodes.get(i);
            long foundIndex = bt.contains(idx) ? idx : -1;
            assertTrue(foundIndex > -1);
            
            long expected = nodes.get(i - 1);
            //System.out.println("\n* " + idx + " expected prev=" + expected);
            
            bt.lower(idx, kOutput);
            assertTrue(kOutput[0] != -1);
            long prev = kOutput[1];
            //System.out.println(idx + "   prev=" + prev);
            assertEquals(expected, prev);
        }
    }
    
    public void estKeyOperations() throws Exception {
    
        System.out.println("testKeyOperations");
        
        RedBlackBSTLongInt bt = new RedBlackBSTLongInt();
        
        int n = 100;
        
        TIntList nodes = new TIntArrayList(2*n);
        
        for (int i = 0; i < n/2; ++i) {
            nodes.add(i);
            bt.put(i, i);
            assertTrue(bt.contains(i));
            assertEquals(0, bt.min());
            assertEquals(i, bt.max());
        }
        for (int i = (n - 1); i >= (n/2); --i) {
            nodes.add(i);
            bt.put(i, i);
            assertTrue(bt.contains(i));
            assertEquals(0, bt.min());
            assertEquals((n - 1), bt.max());
        }
        
        assertEquals(n, bt.size());
                
        assertEquals(n - 1, bt.max());
        
        assertEquals(0L, bt.min());
        
        int[] vOutput = new int[2];
        long[] kOutput = new long[2];
        
        for (int i = 0; i < n; ++i) {
            long foundIndex = bt.contains(i) ? i : -1;
            assertTrue(foundIndex > -1);
            
            /*
            smallest element in the tree with key greater
            than node.key.            
            */
            
            if (i < (n - 1)) {
                bt.higher(i, kOutput);
                assertTrue(kOutput[0] != -1);
                long next = kOutput[1];
                System.out.println(i + " next=" + next);
                assertEquals((i + 1), next);
            }
            /*
            the largest element in the tree with key smaller 
            than node.key
            */
            /*
            if (i > 0) {
                bt.lower(i, kOutput);
                assertTrue(kOutput[0] != -1);
                long prev = kOutput[1];
                System.out.println(i + " prev=" + prev);
                assertEquals((i - 1), prev);
            }*/
        }
        
        // remove some nodes randomly
        TIntSet rm = new TIntHashSet();
        
        SecureRandom sr = Misc0.getSecureRandom();
        long seed = System.currentTimeMillis();
        //seed = 1465940667831L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        bt.delete(nodes.get(0));
        long nod = bt.contains(nodes.get(0)) ? nodes.get(0) : -1;
        assertEquals(-1, nod);
        rm.add(nodes.get(0));
        
        for (int i = 0; i < n/2; ++i) {
            int idx = sr.nextInt(n);
            Integer r = nodes.get(idx);
            if (!rm.contains(r)) {
                bt.delete(r.longValue());
                rm.add(r);
            }
        }
        
        assertEquals((n - rm.size()), bt.size());
                
        boolean minChecked = false;
        
        for (int i = 0; i < n; ++i) {
            Integer index = nodes.get(i);
            long foundIndex = bt.contains(index.longValue()) ? index.longValue() : -1;
            if (rm.contains(index)) {
                assertEquals(-1, foundIndex);
            } else {
                if (!minChecked) {
                    long min = bt.min();
                    assertEquals(index.intValue(), min);
                    minChecked = true;
                }
                assertEquals(index.intValue(), foundIndex);
                if (index.intValue() < (n - 1)) {
                    bt.higher(index, kOutput);
                    assertTrue(kOutput[0] != -1);
                    long next = kOutput[1];
                    int expected = index.intValue() + 1;
                    while (expected < n) {
                        if (rm.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals(expected, next);
                            break;
                        }
                    }
                }
            }
        }

        // ==== then add n more nodes and repeat assertions
        for (int i = n; i < 2*n; ++i) {
            nodes.add(i);
            bt.put(i, i);
        }

        assertEquals((n - rm.size()) + n, bt.size());

        Integer maxExpected = Integer.valueOf(2*n - 1);
        while (rm.contains(maxExpected)) {
            maxExpected = Integer.valueOf(maxExpected.intValue() - 1);
        }
        
        assertEquals(maxExpected.intValue(), bt.max());
        
        for (int i = n; i < 2*n; ++i) {
            Integer index = nodes.get(i);
            long foundIndex = bt.contains(index) ? index.longValue() : -1;
            assertEquals(index.intValue(), foundIndex);
            if (index.intValue() < (n - 1)) {
                bt.higher(index, kOutput);
                assertTrue(kOutput[0] != -1);
                long next = kOutput[1];
                assertEquals(index.intValue() + 1, next);
            }
        }

        for (int i = 0; i < n/2; ++i) {
            int idx = sr.nextInt(n);
            Integer r = nodes.get(idx);
            if (!rm.contains(r)) {
                rm.add(r);
                bt.delete(r.longValue());
                assertFalse(bt.contains(r.longValue()));
            }
        }
        
        assertEquals((2*n - rm.size()), bt.size());

        minChecked = false;
        
        for (int i = 0; i < 2*n; ++i) {
            Integer index = nodes.get(i);
            long foundIndex  = bt.contains(index.longValue()) 
                ? index.longValue() : -1;
            if (rm.contains(index)) {
                assertEquals(-1, foundIndex);
            } else {  
                if (!minChecked) {
                    long min = bt.min();
                    assertEquals(index.intValue(), min);
                    minChecked = true;
                }
                assertEquals(index.intValue(), foundIndex);

                if (index.intValue() < (n - 1)) {
                    bt.higher(index.longValue(), kOutput);
                    assertTrue(kOutput[0] != -1);
                    long next = kOutput[1];
                    int expected = index.intValue() + 1;
                    while (expected < n) {
                        if (rm.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals(expected, next);
                            break;
                        }
                    }
                }
            }
        }
    }
}
