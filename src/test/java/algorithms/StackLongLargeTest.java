package algorithms;

import algorithms.misc.Misc0;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TLongArrayList;
import java.util.Random;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class StackLongLargeTest extends TestCase {
    
    public StackLongLargeTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
        Random rand = Misc0.getSecureRandom();
        
        int n = 500;
        
        StackLongLarge stack = new StackLongLarge(n);
        
        TLongList a = new TLongArrayList();
        
        for (int i = 0; i < n; ++i) {
            long v = rand.nextLong();
            a.add(v);
            stack.push(v);
            assertEquals(v, stack.peek());
            assertEquals(i + 1, stack.size());
            
            assertEquals(i, stack.idxLast);
        }
        
        for (int i = 0; i < n; ++i) {
            int j = n - i - 1;
            
            long v = stack.pop();
            
            //System.out.println("i=" + i + " idxFirst=" 
            //    + stack.idxFirst + " idxLast=" 
            //    + stack.idxLast);
            
            assertEquals(j, stack.size());
            
            assertEquals(a.get(j), v);
            
            if (j > 0) {
                assertEquals(j - 1, stack.idxLast);
            }
        }
        
        assertTrue(stack.isEmpty());
        assertEquals(-1, stack.idxLast);
    
        // add a back in
        // pop 2 items
        // then add a back in
        // then assert contents
        
        for (int i = 0; i < n; ++i) {
            long v = a.get(i);
            stack.push(v);
            assertEquals(v, stack.peek());
            assertEquals(i + 1, stack.size());
        }
        long v = stack.pop();
        assertEquals(a.get(a.size() - 1), v);
        v = stack.pop();
        assertEquals(a.get(a.size() - 2), v);
        TLongList b = new TLongArrayList(a.subList(0, a.size() - 2));
        assertEquals(a.size() - 2, b.size());
        assertEquals(b.size(), stack.size());
        
        for (int i = 0; i < n; ++i) {
            v = a.get(i);
            b.add(v);
            stack.push(v);
            assertEquals(v, stack.peek());
            assertEquals(i + 1 + n - 2, stack.size());
        }
        
        assertEquals(b.size(), stack.size());
        
        // assert contents match b
        for (int i = 0; i < b.size(); ++i) {
            
            int j = b.size() - i - 1;
            
            v = stack.pop();
            
            assertEquals(j, stack.size());
            
            assertEquals(b.get(j), v);
        }
    }
    
}
