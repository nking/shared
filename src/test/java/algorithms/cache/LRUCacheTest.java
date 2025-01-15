package algorithms.cache;

import java.util.HashSet;
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
        float loadFactor = (3 + 2.f)/3;
        boolean accessOrder = true;
        
        algorithms.cache.LRUCache<Integer, Integer> lRUCache = new LRUCache<Integer, Integer>(
            initialCapacity, loadFactor);
        
        for (int i = 0; i < initialCapacity; ++i) {
            lRUCache.put(i, i);
        }
        Iterator<Map.Entry<Integer, Integer>> iter = lRUCache.entrySet().iterator();
        Map.Entry<Integer, Integer> entry;
        int max = -1;
        while (iter.hasNext()) {
            entry = iter.next();
            System.out.printf("lru: key=%d\n", entry.getKey(), entry.getValue());
            if (entry.getKey() > max) {
                max = entry.getKey();
            }
        }
        assertEquals(initialCapacity - 1, max);
        //3,4,5, but access 2, so expecting 4,5,2
        for (int i = initialCapacity; i < 2*initialCapacity; ++i) {
            lRUCache.put(i, i);
            if (i == (initialCapacity+1)) {
                // access methods that affect eviction:
                // put, putIfAbsent, get, getOrDefault, compute, computeIfAbsent, computeIfPresent, or merge
                lRUCache.get(2); // wanting to keep key=2 in map
            }
        }
        HashSet<Integer> expected = new HashSet<Integer>();
        expected.add(4);
        expected.add(5);
        expected.add(2);
        iter = lRUCache.entrySet().iterator();
        while (iter.hasNext()) {
            entry = iter.next();
            System.out.printf("*lru: key=%d\n", entry.getKey(), entry.getValue());
            assert(expected.contains(entry.getKey()));
        }
        assertEquals(initialCapacity, lRUCache.size());
    }
}
