package algorithms.trees;

import junit.framework.TestCase;

import java.util.*;

public class AVLTreeTest extends TestCase {
    /*
               10
             3
           2   5
              4 6

             10
           4
         2   5
              6
     */
    public void testDelFromBST() {
        BinarySearchRecursiveTree<Integer, BinaryTreeNode<Integer>> bst = new BinarySearchRecursiveTree<>();
        bst.insert(10);
        bst.insert(3);
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

        validateOrder(bst);

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

        assertNull(bst.root.parent);
        assertEquals(10, bst.root.left.parent.data.intValue());
        assertEquals(3, bst.root.left.left.parent.data.intValue());
        assertEquals(3, bst.root.left.right.parent.data.intValue());
        assertEquals(5, bst.root.left.right.left.parent.data.intValue());
        assertEquals(5, bst.root.left.right.right.parent.data.intValue());

        List<Integer> before = bst.inOrderTraversal();
        assertEquals(6, before.size());

        int delDatq = 3;
        boolean s = bst.delete(delDatq);
        assertTrue(s);

        List<Integer> after = bst.inOrderTraversal();

        assertTrue(before.size() == after.size()+1);

        int i0 = 0;
        int i1 = 0;
        for (int i = 0; i < before.size(); ++i) {
            if (before.get(i0) == delDatq) {
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

        validateOrder(bst);

    }

    public void testRandomBST() {
        long seed = System.nanoTime();
        //seed = 134911845104821L;
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);
        int nPoints = 1000;

        for (int nTest = 0; nTest < 5; ++nTest) {

            List<Integer> nodes = new ArrayList<>();
            BinarySearchRecursiveTree<Integer, BinaryTreeNode<Integer>> bst = new BinarySearchRecursiveTree<>();
            while (nodes.size() < nPoints) {
                int node = rand.nextInt();
                bst.insert(node);
                nodes.add(node);
            }

            int tmp;
            // shuffle the nodes, use Fisher-Yates
            for (int i = 0; i < nodes.size(); ++i) {
                int j = rand.nextInt(i + 1);
                if (i != j) {
                    tmp = nodes.get(i);
                    nodes.set(i, nodes.get(j));
                    nodes.set(j, tmp);
                }
            }

            //BinaryTreeNode<Integer> dbg = bst.search(934275390);
            // n=2, left=930949343

            // delete each node
            List<Integer> before, after;
            int count = nodes.size();
            for (int j = 0; j < nodes.size(); ++j) {
                int node = nodes.get(j);

                assert(bst.contains(node));

                before = bst.inOrderTraversal();
                assertEquals(count, before.size());

                assertEquals(count, bst.root.n);

               // System.out.printf("nTest=%d, j=%d\n", nTest, j); System.out.flush();

                assertTrue(bst.delete(node));
                after = bst.inOrderTraversal();

                assertTrue(bst.root == null || bst.root.parent == null);

                // compare before and after
                assertTrue(before.size() == after.size()+1);
                int i0 = 0;
                int i1 = 0;
                for (int i = 0; i < before.size(); ++i) {
                    if (before.get(i0) == node) {
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

    /*
               10
             3
           2

           3
         2   10

      ins(4) uses rotR(10)
            3           3
         2   10  -->  2   5
            5            4 10
            4

      ins(6) uses rotL(3)
        5
      3  10
     2 4
     */
    public void testAVL() {

        AVLTree<Integer, AVLTreeNode<Integer>> avl = new AVLTree<>();

        avl.insert(10);
        avl.insert(3);
        assertEquals(10, avl.root.data.intValue());
        assertEquals(3, avl.root.left.data.intValue());

        avl.insert(2);
        // rotate 10 to right
        assertEquals(3, avl.root.data.intValue());
        assertEquals(2, avl.root.left.data.intValue());
        assertEquals(10, avl.root.right.data.intValue());

        avl.insert(5);
        avl.insert(4);
        avl.insert(6);

        validateOrder(avl);

    }

    public static  <T extends Comparable<T>, S extends BinaryTreeNode<T>> void
    validateOrder(BinarySearchRecursiveTree<T, S> tree) {

        BinaryTreeNode<T> node = tree.root;

        if (node == null) {
            return;
        }

        // level order

        // The BSTs use convention for interpreting S.compareTo
        // as left is < 0 else right

        // add assert levels and height

        Queue<BinaryTreeNode<T>> q = new ArrayDeque<>();
        q.add(node);
        while (!q.isEmpty()) {
            node = q.poll();
            long n = 1;
            if (node.left != null) {
                q.add(node.left);
                int comp = node.left.data.compareTo(node.data);
                assertTrue(comp < 0);
                n += node.left.n;
            }
            if (node.right != null) {
                q.add(node.right);
                int comp = node.right.data.compareTo(node.data);
                assertTrue(comp >= 0);
                n += node.right.n;
            }
            assertEquals(node.n, n);
        }
    }

    public void testOrderStatistics() {
        /*
                           7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12
        */
        AVLTree<Integer, AVLTreeNode<Integer>> avl = new AVLTree<>();

        // level order inserts
        avl.insert(7);

        avl.insert(3);
        avl.insert(11);

        avl.insert(1);
        avl.insert(5);
        avl.insert(9);
        avl.insert(13);

        avl.insert(0);
        avl.insert(2);
        avl.insert(4);
        avl.insert(6);
        avl.insert(8);
        avl.insert(10);
        avl.insert(12);

        /*
        long rank(T data);
        long rank(BinaryTreeNode<T> t, T data);
         */
        long expR, ansR;
        int data;
        AVLTreeNode<Integer> node;

        data = 0;
        ansR = avl.rank(data);
        expR = 1;
        assertEquals(expR, ansR);

        data = 13;
        ansR = avl.rank(data);
        expR = 14;
        assertEquals(expR, ansR);

        data = 6;
        ansR = avl.rank(data);
        expR = 7;
        assertEquals(expR, ansR);

        data = 0;
        ansR = avl.rank(avl.root, data);
        expR = 1;
        assertEquals(expR, ansR);

        data = 4;
        ansR = avl.rank(avl.search(5), data);
        expR = 1;
        assertEquals(expR, ansR);

        data = 9;
        ansR = avl.rank(avl.search(11), data);
        expR = 2;
        assertEquals(expR, ansR);

        /*
        S select(long rank);
        S select(BinaryTreeNode<T> t, long rank);
         */
        long rank;
        AVLTreeNode<Integer> ansNode;

        rank = 1;
        ansNode = avl.select(rank);
        assertTrue(ansNode.data == 0);

        rank = 14;
        ansNode = avl.select(rank);
        assertTrue(ansNode.data == 13);

        rank = 1;
        ansNode = avl.select(avl.search(5), rank);
        assertTrue(ansNode.data == 4);

        rank = 3;
        ansNode = avl.select(avl.search(5), rank);
        assertTrue(ansNode.data == 6);

        rank = 2;
        ansNode = avl.select(avl.search(11), rank);
        assertTrue(ansNode.data == 9);

        rank = 12;
        ansNode = avl.select(avl.search(11), rank);
        assertNull(ansNode);
    }
}
