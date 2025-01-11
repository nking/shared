package algorithms.trees;

import junit.framework.TestCase;

public class SegmentTreeTest extends TestCase {

    public void test0() {
        int[] a = new int[]{5,8,6,3,2,7,2,6};
        SegmentTree st = new SegmentTree(a);

        long sum = st.sum(2, 7);
        assertEquals(26, sum);

        sum = st.sum(3, 5);
        assertEquals(12, sum);

        sum = st.sum(0, 2);
        assertEquals(19, sum);

/*
        // change i=0 to 6
        st.updateSet(0, 6);
        sum = st.sum(0, 2);
        assertEquals(20, sum);

        sum = st.sum(2, 7);
        assertEquals(26, sum);

        sum = st.sum(3, 5);
        assertEquals(12, sum);

        sum = st.sum(0, 2);
        assertEquals(19, sum);

        // change i=2 to 15
        //6, 8,15, 3, 2, 7, 2, 6
        //0  1  2  3  4  5  6  7
        st.updateSet(2, 15);
        sum = st.sum(0, 2);
        assertEquals(29, sum);

        sum = st.sum(3, 6);
        assertEquals(14, sum);
*/
    }
}
