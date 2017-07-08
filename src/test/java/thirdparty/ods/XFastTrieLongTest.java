package thirdparty.ods;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;

/**
 *
 * @author nichole
 */
public class XFastTrieLongTest extends TestCase {
    
    public XFastTrieLongTest() {
    }

    public void test0() {
        
        System.out.println("test0");
        
		XFastTrieNodeLong<Integer> node = new XFastTrieNodeLong<Integer>();
		
        Longizer<Integer> it = new Longizer<Integer>() {
            @Override
            public long longValue(Integer x) {
                return x;
            }
        };
        
        int w = 4;
		
        XFastTrieLong<XFastTrieNodeLong<Integer>, Integer> bt
            = new XFastTrieLong<XFastTrieNodeLong<Integer>, Integer>(node, it, w);
        
        assertEquals(0, bt.size());
  
        System.out.println("add 2 0010");
        
        boolean added = bt.add(Integer.valueOf(2));
        assertTrue(added);
        
        assertEquals(1, bt.size());
                
        System.out.println("add 3 0011");
        
        added = bt.add(Integer.valueOf(3));
        assertTrue(added);
                
        assertEquals(2, bt.size());
        
        System.out.println("add 5 0101");
         
        added = bt.add(Integer.valueOf(5));
        assertTrue(added);
                
        assertEquals(3, bt.size());
        
        System.out.println("add 4 0100");
       
        added = bt.add(Integer.valueOf(4));
        assertTrue(added);
       
        assertEquals(4, bt.size());
        
        bt.debugNodes();
        
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
        
        Integer next = bt.successor(4);
        assertEquals(5, it.longValue(next));
        
        next = bt.successor(0);
        assertEquals(2, it.longValue(next));
        
        next = bt.successor(2);
        assertEquals(3, it.longValue(next));
        
        next = bt.successor(3);
        assertEquals(4, it.longValue(next));
        
        next = bt.successor(4);
        assertEquals(5, it.longValue(next));
        
        Integer prev;
       
        prev = bt.predecessor((1<<(w-1)) - 1);
        assertEquals(5, it.longValue(prev));
        
        prev = bt.predecessor(6);
        assertEquals(5, it.longValue(prev));
       
        prev = bt.predecessor(5);
        assertEquals(4, it.longValue(prev));
        
        prev = bt.predecessor(4);
        assertEquals(3, it.longValue(prev));
        
        prev = bt.predecessor(3);
        assertEquals(2, it.longValue(prev));
        
        prev = bt.predecessor(2);
        assertNull(prev);
        
        assertEquals(2, bt.minimum().intValue());
        assertEquals(5, bt.maximum().intValue());
        
        assertTrue(bt.remove(3));
        assertEquals(3, bt.size());
        assertNull(bt.find(Integer.valueOf(3)));
        assertEquals(2, bt.find(Integer.valueOf(2)).intValue());
        assertEquals(4, bt.find(Integer.valueOf(4)).intValue());
        assertEquals(5, bt.find(Integer.valueOf(5)).intValue());

        assertNull(bt.find(Integer.valueOf(0)));
        
        added = bt.add(Integer.valueOf(3));
        assertTrue(added);
        assertEquals(3, bt.find(Integer.valueOf(3)).intValue());
        
        assertTrue(bt.remove(5));
        assertNull(bt.find(Integer.valueOf(5)));
        assertTrue(bt.remove(4));
        assertNull(bt.find(Integer.valueOf(4)));
        
        prev = bt.predecessor((1<<(w-1)) - 1);
        assertEquals(3, it.longValue(prev));
        
        assertEquals(2, bt.minimum().intValue());
        assertEquals(3, bt.maximum().intValue());
    }
        
    public void test1() throws Exception {
    
        System.out.println("test1");
        
		XFastTrieNodeLong<Integer> node = new XFastTrieNodeLong<Integer>();
		
        Longizer<Integer> it = new Longizer<Integer>() {
            @Override
            public long longValue(Integer x) {
                return x;
            }
        };
		
        XFastTrieLong<XFastTrieNodeLong<Integer>, Integer> bt
            = new XFastTrieLong<XFastTrieNodeLong<Integer>, Integer>(node, it);
        
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
        
        assertEquals(n, bt.size());
                
        assertEquals(n - 1, bt.maximum().intValue());
        
        assertEquals(0L, bt.minimum().intValue());
        
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
                assertEquals((index.intValue() + 1), it.longValue(next));
            }
            /*
            the largest element in the tree with key smaller 
            than node.key
            */
            if (index.intValue() > 0) {
                Integer prev = bt.predecessor(index);
                assertEquals((index.intValue() - 1), it.longValue(prev));
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
                assertTrue(bt.remove(r));
                rm.add(r);
            }
        }
        
        assertEquals((n - rm.size()), bt.size());
                
        boolean minChecked = false;
        
        for (int i = 0; i < n; ++i) {
            Integer index = nodes.get(i);
            Integer foundIndex = bt.find(index);
            if (rm.contains(index)) {
                assertNull(foundIndex);
            } else {
                if (!minChecked) {
                    Integer min = bt.minimum();
                    assertEquals(index.intValue(), min.intValue());
                    minChecked = true;
                }
                assertEquals(index.intValue(), foundIndex.intValue());
                if (index.intValue() < (n - 1)) {
                    Integer next = bt.successor(index);
                    int expected = index.intValue() + 1;
                    while (expected < n) {
                        if (rm.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals(expected, it.longValue(next));
                            break;
                        }
                    }
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
        
        assertEquals(maxExpected.intValue(), bt.maximum().intValue());
        
        for (int i = n; i < 2*n; ++i) {
            Integer index = nodes.get(i);
            Integer foundIndex = bt.find(index);
            assertEquals(index.intValue(), foundIndex.intValue());
            if (index.intValue() < (n - 1)) {
                Integer next = bt.successor(index);
                assertEquals(index.intValue() + 1, it.longValue(next));
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

        minChecked = false;
        
        for (int i = 0; i < 2*n; ++i) {
            Integer index = nodes.get(i);
            Integer foundIndex  = bt.find(index);
            if (rm.contains(index)) {
                assertNull(foundIndex);
            } else {  
                if (!minChecked) {
                    Integer min = bt.minimum();
                    assertEquals(index.intValue(), min.intValue());
                    minChecked = true;
                }
                assertEquals(index.intValue(), foundIndex.intValue());

                if (index.intValue() < (n - 1)) {
                    Integer next = bt.successor(index);
                    int expected = index.intValue() + 1;
                    while (expected < n) {
                        if (rm.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals(expected, it.longValue(next));
                            break;
                        }
                    }
                }
            }
        }
    }
   
    public void test3() throws Exception {
    
        System.out.println("test3");
        
        XFastTrieNodeLong<Long> node = new XFastTrieNodeLong<Long>();
		
        Longizer<Long> it = new Longizer<Long>() {
            @Override
            public long longValue(Long x) {
                return x;
            }
        };
        
        int w = 61;
		
        XFastTrieLong<XFastTrieNodeLong<Long>, Long> bt
            = new XFastTrieLong<XFastTrieNodeLong<Long>, Long>(node, it, w);
        
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
