package thirdparty.ods;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BinaryTrieLongTest extends TestCase {
    
    public BinaryTrieLongTest() {
    }

    public static class BTN2 extends BinaryTrieNode<Integer> {
    };
    
    public void test0() {
        
        System.out.println("test0");
        
		BTN2 node = new BTN2();
		
        Longizer<Integer> it = new Longizer<Integer>() {
            @Override
            public long longValue(Integer x) {
                return x;
            }
        };
		
        BinaryTrieLong<BTN2, Integer> bt
            = new BinaryTrieLong<BTN2, Integer>(node, it, 4);
        
        assertEquals(0, bt.size());
        
        boolean added = bt.add(Integer.valueOf(2));
        assertTrue(added);
        
        assertEquals(1, bt.size());
                
        added = bt.add(Integer.valueOf(3));
        assertTrue(added);
                
        assertEquals(2, bt.size());
        
        added = bt.add(Integer.valueOf(5));
        assertTrue(added);
                
        assertEquals(3, bt.size());
        
        added = bt.add(Integer.valueOf(4));
        assertTrue(added);
       
        assertEquals(4, bt.size());
        
        //bt.debugNodes();
        
        /*
        add "2", 0010
                                               r
                               0(pad)                   null
                       0(pad)       
                           1                  
                          0            
        -----------
        add "2", 0010
        add "3", 0011
                                               r
                               0(pad)                   null
                       0(pad)       
                           1                  
                         0   1                  
        ------------
        add "2", 0010
        add "3", 0011
        add "5", 0101
                                               r
                               0(pad)                   null
                       0(pad)              1     
                           1           0    
                         0   1           1
        ------------
        add "2", 0010
        add "3", 0011
        add "5", 0101
        add "4", 0100
                                               r
                               0(pad)                   null
                       0(pad)              1     
                           1           0    
                         0   1       0   1
        */
        
        assertEquals(3, bt.find(Integer.valueOf(3)).intValue());
        assertEquals(2, bt.find(Integer.valueOf(2)).intValue());
        assertEquals(4, bt.find(Integer.valueOf(4)).intValue());
        assertEquals(5, bt.find(Integer.valueOf(5)).intValue());
        assertTrue(bt.remove(3));
        assertEquals(3, bt.size());
        assertNull(bt.find(Integer.valueOf(3)));
        assertEquals(2, bt.find(Integer.valueOf(2)).intValue());
        assertEquals(4, bt.find(Integer.valueOf(4)).intValue());
        assertEquals(5, bt.find(Integer.valueOf(5)).intValue());

        assertNull(bt.find(Integer.valueOf(0)));
    
        bt.debugNodes();
        System.out.println("bt.n=" + bt.n);
    }
        
    public void test1() throws Exception {
    
        System.out.println("test1");
        
		BinaryTrieNode<Integer> node = new BinaryTrieNode<Integer>();
		
        Longizer<Integer> it = new Longizer<Integer>() {
            @Override
            public long longValue(Integer x) {
                return x;
            }
        };
		
        BinaryTrieLong<BinaryTrieNode<Integer>, Integer> bt
            = new BinaryTrieLong<BinaryTrieNode<Integer>, 
                Integer>(node, it);
        
        int n = 100;
        
        List<Integer> nodes = new ArrayList<Integer>(2*n);
        
        for (int i = 0; i < n/2; ++i) {
            Integer ii = Integer.valueOf(i);
            nodes.add(ii);
            bt.add(ii);
        }
        for (int i = (n - 1); i >= (n/2); --i) {
            Integer ii = Integer.valueOf(i);
            nodes.add(ii);
            bt.add(ii);
        }
        
        bt.debugNodes();
        System.out.println("bt.n=" + bt.n);
        
        assertEquals(n, bt.size());
        
        for (int i = 0; i < n; ++i) {
            Integer index = nodes.get(i);
            Integer foundIndex = bt.find(index);
            assertNotNull(foundIndex);
            assertEquals(index.intValue(), foundIndex.intValue());
            
            /*
            smallest element in the tree with key greater
            than node.key.            
            */
            
            if (index.intValue() < (n - 1)) {
                Integer next = bt.successor(index);
                assertEquals((index.longValue() + 1), 
                    it.longValue(next));
            }
            
            /*
            the largest element in the tree with key smaller 
            than node.key
            */
            if (index.intValue() > 0) {
                Integer prev = bt.predecessor(index);
                assertEquals((index.longValue() - 1), 
                    it.longValue(prev));
            }
        }
        
        // remove some nodes randomly
        Set<Integer> rm = new HashSet<Integer>();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1465940667831L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        bt.remove(nodes.get(0));
        Integer nod = bt.find(nodes.get(0));
        assertNull(nod);
        rm.add(nodes.get(0));
        
        for (int i = 0; i < n/2; ++i) {
            int idx = sr.nextInt(n);
            Integer r = nodes.get(idx);
            if (!rm.contains(r)) {
                //System.out.println("==> removing " + r);                
                assertTrue(bt.remove(r));
                rm.add(r);
            }
        }
        
        assertEquals((n - rm.size()), bt.size());
   
        //System.out.println("--2-- n="+n);
        //bt.debugNodes();
        
        for (int i = 0; i < n; ++i) {
            Integer index = nodes.get(i);
            Integer foundIndex = bt.find(index);
            if (rm.contains(index)) {
                assertNull(foundIndex);
            } else {
                assertEquals(index.intValue(), foundIndex.intValue());
                int expectedNext = index.intValue() + 1;
                while (rm.contains(expectedNext)) {
                    expectedNext++;
                }
                if ((expectedNext < n) && (index.intValue() < (n - 1))) {
                    Integer next = bt.successor(index);
                    //System.out.println("successor of index(" + index + ")=" + next);
                    assertEquals(expectedNext, it.longValue(next));
                }
            }
        }

        // ==== then add n more nodes and repeat assertions
        for (int i = n; i < 2*n; ++i) {
            Integer ii = Integer.valueOf(i);
            nodes.add(ii);
            bt.add(ii);
        }

        assertEquals((n - rm.size()) + n, bt.size());

        Integer maxExpected = Integer.valueOf(2*n - 1);
        while (rm.contains(maxExpected)) {
            maxExpected = Integer.valueOf(maxExpected.intValue() - 1);
        }
        
        for (int i = n; i < 2*n; ++i) {
            Integer index = nodes.get(i);
            Integer foundIndex = bt.find(index);
            assertEquals(index.intValue(), foundIndex.intValue());
            if (index.intValue() < (n - 1)) {
                Integer next = bt.successor(index);
                assertEquals(index.longValue() + 1, 
                    it.longValue(next));
            }
        }

        for (int i = 0; i < n/2; ++i) {
            int idx = sr.nextInt(n);
            Integer r = nodes.get(idx);
            if (!rm.contains(r)) {
                rm.add(r);
                bt.remove(r);
                assertNull(bt.find(r));
            }
        }
        
        assertEquals((2*n - rm.size()), bt.size());

        //System.out.println("--3--");
        //bt.debugNodes();
        
        for (int i = 0; i < 2*n; ++i) {
            Integer index = nodes.get(i);
            Integer foundIndex  = bt.find(index);
            if (rm.contains(index)) {
                assertNull(foundIndex);
            } else {              
                assertEquals(index.intValue(), foundIndex.intValue());
                int expectedNext = index.intValue() + 1;
                while (rm.contains(expectedNext)) {
                    expectedNext++;
                }
                if ((expectedNext < n) && (index.intValue() < (n - 1))) {
                    Integer next = bt.successor(index);
                    //System.out.println("successor of index(" + index + ")=" + next);
                    assertEquals(expectedNext, it.longValue(next));
                }
            }
        }
    }
    
    public void test3() throws Exception {
    
        System.out.println("test3");
        
		BinaryTrieNode<Long> node = new BinaryTrieNode<Long>();
		
        Longizer<Long> it = new Longizer<Long>() {
            @Override
            public long longValue(Long x) {
                return x;
            }
        };
		
        BinaryTrieLong<BinaryTrieNode<Long>, Long> bt
            = new BinaryTrieLong<BinaryTrieNode<Long>, 
                Long>(node, it);
        
        int n = 100;
        
        List<Long> nodes = new ArrayList<Long>(2*n);
        
        long start = (1L<<32);
        
        for (long i = start; i < (start + n/2); ++i) {
            Long ii = Long.valueOf(i);
            nodes.add(ii);
            bt.add(ii);
        }
        for (long i = (start + n - 1); i >= start + (n/2); --i) {
            Long ii = Long.valueOf(i);
            nodes.add(ii);
            bt.add(ii);
        }
       
        bt.debugNodes();
        System.out.println("bt.n=" + bt.n);
        
        for (long i = start; i < (start + n); ++i) {
            Long ii = Long.valueOf(i);
            Long f = bt.find(ii);
            assertEquals(ii, f);
        }
        
        for (long i = start + 1; i < (start + n - 1); ++i) {
            Long ii = Long.valueOf(i);
            
            Long p = bt.predecessor(ii);
            Long s = bt.successor(ii);
            
            assertEquals(i - 1, p.longValue());
            
            assertEquals(i + 1, s.longValue());
        }
        
        assertEquals(start, bt.minimum().longValue());
        
        assertEquals(start + n - 1, bt.maximum().longValue());
    }
}
