package algorithms.range;

import algorithms.trees.SegmentTree;
import junit.framework.TestCase;

public class PrefixSumArrayTest extends TestCase {

    public void test0() {

        boolean use0BasedIndexes  = true;
        long[] a = new long[]{5,8,6,3,2,7,2,6};

        PrefixSumArray psa = new PrefixSumArray(a);
        assertEquals(26, psa.sum(2, 7, use0BasedIndexes));

        int[][] updates = new int[][]{
                {2,4,3},
                {1,5,4},
                {5,7,1},
                {3,3,3}
        };
        psa.update(updates, use0BasedIndexes);
        assertEquals(74, psa.sum(0, 7, use0BasedIndexes));
        assertEquals(24, psa.sum(4, 6, use0BasedIndexes));

    }
}
