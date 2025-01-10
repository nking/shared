package algorithms.trees;

import junit.framework.TestCase;

import java.util.Arrays;

public class LeastCommonAncestorTest extends TestCase {

    public void test0() {
        int[] a = new int[]{6,9,2,4,7,8,5,8,3,7};

        LeastCommonAncestor lca = new LeastCommonAncestor(a);

        int[] expTree = new int[]{
                2, 0, -1, 8, 6, 4, 3, 6, 2, 8
        };

        assertTrue(Arrays.equals(expTree, lca.tree));

        int left = 6;
        int right = 9;
        int ans, expAns;
        expAns = 8;
        ans = lca.findWithLogN(left, right);
        assertEquals(expAns, ans);
        ans = lca.find(left, right);
        assertEquals(expAns, ans);

        left = 1;
        right = 4;
        expAns = 2;
        ans = lca.findWithLogN(left, right);
        assertEquals(expAns, ans);
        ans = lca.find(left, right);
        assertEquals(expAns, ans);
    }
}
