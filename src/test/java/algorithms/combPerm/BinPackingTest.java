package algorithms.combPerm;

import junit.framework.TestCase;

public class BinPackingTest extends TestCase {

    public void test0() {
        int[] a;
        int ans, expAns;
        int cap;

        cap = 10;
        expAns = 2;
        a = new int[]{3, 2, 3, 5, 6};
        ans = BinPacking.minNumberOfBins(a, cap);
        assertEquals(expAns, ans);
    }
}
