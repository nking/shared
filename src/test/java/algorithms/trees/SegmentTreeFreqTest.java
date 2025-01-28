package algorithms.trees;

import junit.framework.TestCase;
import java.util.*;

public class SegmentTreeFreqTest extends TestCase {

    public void test0() {
        int[] a = new int[]{3,1,2,3,1,1,1,2};
        SegmentTreeFreq st = new SegmentTreeFreq(a);

        int[][] queries = new int[][] {
                {0,7}
        };

        List<Map<Integer, Integer>> ans = st.query(queries);
        assertEquals(queries.length, ans.size());
        for (int i = 0; i < queries.length; ++i) {
            Map<Integer, Integer> expected = createCountMap(a, queries[i]);
            Map<Integer, Integer> r = ans.get(i);
            assertEquals(expected.size(), r.size());
            for (Map.Entry<Integer, Integer> entry : r.entrySet()) {
                assertNotNull(expected.remove(entry.getKey()));
            }
            assertTrue(expected.isEmpty());
        }
    }

    protected Map<Integer, Integer> createCountMap(int[] a, int[] q) {
        Map<Integer, Integer> countMap = new HashMap<>();
        for (int i = q[0]; i<= q[1]; ++i) {
            countMap.put(a[i], countMap.getOrDefault(a[i], 0) + 1);
        }
        return countMap;
    }
}
