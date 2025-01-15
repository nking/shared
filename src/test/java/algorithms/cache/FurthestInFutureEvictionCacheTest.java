package algorithms.cache;

import junit.framework.TestCase;

import java.util.Arrays;

public class FurthestInFutureEvictionCacheTest extends TestCase {

    public void test0() {
        int[] blocks = new int[]{1,5,2,3,7,8,3,1,2,3, 7,4};
        int cacheSize = 4;

        /*
        cache
        1, 5, 2, 3, 7, 8, 3, 1, 2, 3, 7, 4
        1, 5, 2, 3   comp misses = 4
                    * misses = 1, evict 5
        1, 2, 3, 7
                      * misses = 2, evict 7
        1, 2, 3, 8        * hits = 1
                             *  *  * hits=4
                 7                   * misses=3
                 4                     misses = 4
         */
        int[] expected = new int[]{4, 4, 4};

        int[] counts = FurthestInFutureEvictionCache.predict(blocks, cacheSize);
        assertTrue(Arrays.equals(expected, counts));
    }
}
