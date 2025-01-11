package algorithms.range;

import algorithms.trees.SegmentTree;
import junit.framework.TestCase;

public class SegmentTreeTest extends TestCase {

    public void test0() {
        int[] a = new int[]{5,8,6,3,2,7,2,6};
        boolean use0BasedIndexes = false;
        SegmentTree st = new SegmentTree(a, SegmentTree.TYPE.SUM, use0BasedIndexes);

        long sum = st.sum(2, 7);
        assertEquals(26, sum);

        use0BasedIndexes = true;
        st = new SegmentTree(a, SegmentTree.TYPE.SUM, use0BasedIndexes);

        sum = st.sum(2, 7);
        assertEquals(26, sum);

    }
}
