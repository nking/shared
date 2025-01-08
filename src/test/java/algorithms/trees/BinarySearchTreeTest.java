package algorithms.trees;

import junit.framework.TestCase;

import java.util.*;

public class BinarySearchTreeTest extends TestCase {
    /*
                    10
              3
          2     5
              4   6

             10
           4
         2   5
              6
     */
    public void testDel() {
        BinarySearchTree<Integer, BinaryTreeNode<Integer>> bst = new BinarySearchTree<>();
        bst.insert(10);
        BinaryTreeNode<Integer> node3 = bst.insert(3);
        bst.insert(2);
        bst.insert(5);
        bst.insert(4);
        bst.insert(6);

        assertEquals(10, bst.root.data.intValue());
        assertEquals(3, bst.root.left.data.intValue());
        assertEquals(2, bst.root.left.left.data.intValue());
        assertEquals(5, bst.root.left.right.data.intValue());
        assertEquals(4, bst.root.left.right.left.data.intValue());
        assertEquals(6, bst.root.left.right.right.data.intValue());

        assert(bst.contains(10));
        assert(bst.contains(3));
        assert(bst.contains(2));
        assert(bst.contains(5));
        assert(bst.contains(4));
        assert(bst.contains(6));

        assertEquals(6L, bst.root.n);
        assertEquals(5L, bst.root.left.n);
        assertEquals(1L, bst.root.left.left.n);
        assertEquals(3L, bst.root.left.right.n);
        assertEquals(1L, bst.root.left.right.left.n);
        assertEquals(1L, bst.root.left.right.right.n);

        List<Integer> before = bst.inOrderTraversal();
        assertEquals(6, before.size());

        boolean s = bst.delete(node3);
        assertTrue(s);

        List<Integer> after = bst.inOrderTraversal();

        assertTrue(before.size() == after.size()+1);

        int i0 = 0;
        int i1 = 0;
        for (int i = 0; i < before.size(); ++i) {
            if (before.get(i0) == node3.data) {
                ++i0;
                continue;
            }
            assertEquals(before.get(i0), after.get(i1));
            ++i0;
            ++i1;
        }

        assertEquals(5L, bst.root.n);
        assertEquals(4L, bst.root.left.n);
        assertEquals(1L, bst.root.left.left.n);
        assertEquals(2L, bst.root.left.right.n);
        assertEquals(1L, bst.root.left.right.right.n);

        /*
                    10
              3
          2     5
              4   6

             10
           4
         2   5
              6
     */
    }

    public void testRandom() {
        long seed = System.nanoTime();
        Random rand = new Random(seed);
        int nPoints = 1000;

        for (int nTest = 0; nTest < 10; ++nTest) {
            List<BinaryTreeNode<Integer>> nodes = new ArrayList<>();
            BinarySearchTree<Integer, BinaryTreeNode<Integer>> bst = new BinarySearchTree<>();
            while (nodes.size() < nPoints) {
                nodes.add(bst.insert(rand.nextInt()));
            }

            BinaryTreeNode<Integer> tmp;
            // shuffle the nodes, use Fisher-Yates
            for (int i = 0; i < nodes.size(); ++i) {
                int j = rand.nextInt(i + 1);
                if (i != j) {
                    tmp = nodes.get(i);
                    nodes.set(i, nodes.get(j));
                    nodes.set(j, tmp);
                }
            }

            // delete each node
            List<Integer> before, after;
            int count = nodes.size();
            for (BinaryTreeNode<Integer> node : nodes) {

                assert(bst.contains(node.data));

                before = bst.inOrderTraversal();
                assertEquals(count, before.size());
                assertTrue(bst.delete(node));
                after = bst.inOrderTraversal();

                // compare before and after
                assertTrue(before.size() == after.size()+1);
                int i0 = 0;
                int i1 = 0;
                for (int i = 0; i < before.size(); ++i) {
                    if (before.get(i0) == node.data) {
                        ++i0;
                        continue;
                    }
                    assertEquals(before.get(i0), after.get(i1));
                    ++i0;
                    ++i1;
                }

                --count;
            }

            assertNull(bst.root);
        }
    }

}
