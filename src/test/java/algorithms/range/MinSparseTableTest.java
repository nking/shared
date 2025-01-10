package algorithms.range;

import junit.framework.TestCase;

import java.util.Arrays;

public class MinSparseTableTest extends TestCase {

    public void test0() {
        int[] a = new int[]{1,3,4,8,6,1,4,2};
        MinSparseTable min = new MinSparseTable(a);

        int[][] queries = new int[][]{
                {3,3}, {1,2}, {3,4},{1,4}
        };
        int[] ans = min.min(queries);
        int[] expAns = new int[]{8,3,6,3};
        assertTrue(Arrays.equals(expAns, ans));
    }
}
