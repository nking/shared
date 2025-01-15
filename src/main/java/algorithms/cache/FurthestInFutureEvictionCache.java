package algorithms.cache;

import java.util.*;

/**
 * given a list of block ids that will be requested in that order,
 * and given a cache with a size limit, when the cache is full and
 * experiences a miss, this algorithm chooses the cache element to
 * evict as the index in the cache with a next use index furthest in future than the others' in the cache
 */
public class FurthestInFutureEvictionCache {

    /**
     * given the future block ids and the cache size, return the number of hits and misses from following
     * the furthest in future cache evict policy.
     * r.t.c. is O(n * cacheSize) where n is length of block indexes.
     * @param blockIdxs the block numbers that will be requested, in order of time
     * @param cacheSize the number of blocks that can fit in the cache.
     * @return an array of length 3 holding:
     <pre>
     row 0 contains the compulsory cache misses (the number of hits to fill the cache.
     row 1 contains the number of cache misses minus the compulsory.
     row 2 contains the number of cache hits.
    </pre>
     */
    public static int[] predict(int[] blockIdxs, int cacheSize) {

        Map<Integer, LinkedList<Integer>> blockIndexMap = new HashMap<>();
        for (int i = 0; i < blockIdxs.length; ++i) {
            blockIndexMap.putIfAbsent(blockIdxs[i], new LinkedList<>());
            blockIndexMap.get(blockIdxs[i]).add(i);
        }

        // row 0 is compulsory cache misses
        // row 1 is cache misses - compulsory
        // row 2 is cache hits
        int[] counts = new int[3];
        Set<Integer> cache = new HashSet<>(cacheSize);

        int b;
        for (int i = 0; i < blockIdxs.length; ++i) {
            b = blockIdxs[i];
            if (cache.size() < cacheSize) {
                cache.add(b);
                ++counts[0]; // compulsory miss
            } else {
                if (cache.contains(b)) {
                    ++counts[2]; // a hit
                } else {
                    ++counts[1]; // miss

                    // find key to evict
                    int evictIdx = Integer.MIN_VALUE;
                    int evictKey = -1;
                    for (int key : cache) {
                        if (!blockIndexMap.containsKey(key)) {
                            // no future use, so choose it
                            evictKey = key;
                            break;
                        }
                        int idx = blockIndexMap.get(key).peekFirst();
                        if (idx > evictIdx) {
                            evictIdx = idx;
                            evictKey = key;
                        }
                    }
                    assert(evictKey != -1);
                    cache.remove(evictKey);

                    // add b
                    cache.add(b);
                }
            }
            // remove i from blockIndexMap
            assert(blockIndexMap.containsKey(b));
            assert(blockIndexMap.get(b).peekFirst() == i);
            if (blockIndexMap.get(b).size() == 1) {
                blockIndexMap.remove(b);
            } else {
                blockIndexMap.get(b).pollFirst();
            }
        }

        return counts;
    }
}
