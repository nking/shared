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

    public void test1() {
        NAryTreeNode root = new NAryTreeNode(2);
        NAryTreeNode ch1 = new NAryTreeNode(3);
        NAryTreeNode ch2 = new NAryTreeNode(4);
        NAryTreeNode ch2_1 = new NAryTreeNode(6);
        ch2.addChild(ch2_1);

        root.addChild(ch1);
        root.addChild(ch2);

        NAryTreeNode node1 = new NAryTreeNode(9);
        NAryTreeNode node2 = new NAryTreeNode(5);

        ch1.addChild(node1);
        ch2_1.addChild(node2);

        /*
                root
            ch1      ch2
          node1          ch2_1
                             node2
         */

        long dist = LeastCommonAncestor.distBetweenNodes(root, node1, node2);
        assertEquals(5L, dist);
        dist = LeastCommonAncestor.distBetweenNodes(root, node1.getData(), node2.getData());
        assertEquals(5L, dist);

        /*
                root
            ch1        ch2
                   node1 ch2_1
                           node2*/
        ch1.getChildren().clear();
        ch2.addChild(node1);
        dist = LeastCommonAncestor.distBetweenNodes(root, node1, node2);
        assertEquals(3L, dist);
        dist = LeastCommonAncestor.distBetweenNodes(root, node1.getData(), node2.getData());
        assertEquals(3L, dist);

        // ====
        dist = LeastCommonAncestor.distBetweenNodes(root, root, root);
        assertEquals(0L, dist);
        dist = LeastCommonAncestor.distBetweenNodes(root, root.getData(), root.getData());
        assertEquals(0L, dist);

        // ====
        dist = LeastCommonAncestor.distBetweenNodes(root, ch2_1, ch2_1);
        assertEquals(0L, dist);
        dist = LeastCommonAncestor.distBetweenNodes(root, ch2_1.getData(), ch2_1.getData());
        assertEquals(0L, dist);

    }
}
