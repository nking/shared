package algorithms.range;

import junit.framework.TestCase;

import java.util.Arrays;

public class MiscTest extends TestCase {

    public void test0() {

        int[] a = new int[]{0, 1, 2, 3, 4};
        int[][] updates = new int[][] {
                {1,3, -1},
                {3,5, +3}
        };

        int[] expA = new int[] {0, 0, 1, 5, 7};

        Misc.updateAddUsingDifferenceArray(a, updates);

        assertTrue(Arrays.equals(expA, a));
    }
}
