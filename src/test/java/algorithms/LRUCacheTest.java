package algorithms;

import java.util.Iterator;
import java.util.Map;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LRUCacheTest extends TestCase {
    
    public LRUCacheTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        int initialCapacity = 3;
        float loadFactor = 0.75f;
        boolean accessOrder = true;
        
        LRUCache<Integer, Integer> lRUCache = new LRUCache<Integer, Integer>(
            initialCapacity, loadFactor);
        
        for (int i = 0; i < initialCapacity; ++i) {
            lRUCache.put(i, i);
        }
        Iterator<Map.Entry<Integer, Integer>> iter = lRUCache.entrySet().iterator();
        Map.Entry<Integer, Integer> entry;
        while (iter.hasNext()) {
            entry = iter.next();
            System.out.printf("lru: key=%d\n", entry.getKey(), entry.getValue());
        }
        for (int i = initialCapacity; i < 2*initialCapacity; ++i) {
            lRUCache.put(i, i);
            if (i == (initialCapacity+1)) {
                // access methods that affect eviction:
                // put, putIfAbsent, get, getOrDefault, compute, computeIfAbsent, computeIfPresent, or merge
                lRUCache.get(2); // wanting to keep key=2 in map
            }
        }
        iter = lRUCache.entrySet().iterator();
        while (iter.hasNext()) {
            entry = iter.next();
            System.out.printf("*lru: key=%d\n", entry.getKey(), entry.getValue());
        }
    }
}
