package algorithms.trees;

import algorithms.range.PrefixSumArray;
import junit.framework.TestCase;

import java.util.Random;

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

        //=====================
        // change i=0 to 6
        st.updateSet(0, 6);

        sum = st.sum(0, 2);
        assertEquals(20, sum);

        sum = st.sum(2, 7);
        assertEquals(26, sum);

        sum = st.sum(3, 5);
        assertEquals(12, sum);

        // change i=2 to 15
        //6, 8,15, 3, 2, 7, 2, 6
        //0  1  2  3  4  5  6  7
        st.updateSet(2, 15);
        sum = st.sum(0, 2);
        assertEquals(29, sum);

        sum = st.sum(3, 6);
        assertEquals(14, sum);

    }

    public void test1() {
        long seed = System.nanoTime();
        //seed = 262023598709905L;
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        int nTests = 10;
        int nMax = 64;
        int nMin = 3;
        for (int nTest = 0; nTest < nTests;++nTest) {

            int n = nMin + rand.nextInt(nMax - nMin);
            int[] a = new int[n];
            for (int i = 0; i < n; ++i) {
                a[i] = rand.nextInt();
            }

            long[] prefixA = prefix(a);

            SegmentTree st = new SegmentTree(a);
            int i0, i1;
            long sum;
            for (int j = 0; j < nTests; ++j) {
                i0 = rand.nextInt(n); // 9; 10-9=1-1
                i1 = i0 + rand.nextInt(n - i0);

                sum = st.sum(i0, i1);

                assertEquals(prefixA[i1+1] - prefixA[i0], sum);
            }
        }
    }

    protected long[] prefix(int[] a) {
        long[] out = new long[a.length + 1];
        for (int i = 0; i < a.length; ++i) {
            out[i+1] = a[i];
        }
        for (int i = 1; i < out.length; ++i) {
            out[i] += out[i-1];
        }

        return out;
    }
}
