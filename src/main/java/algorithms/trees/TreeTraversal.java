package algorithms.trees;

import algorithms.DoublyLinkedList;
import algorithms.heapsAndPQs.HeapNode;

import java.util.*;

/**
 *
 * @author nichole
 */
public class TreeTraversal {

    /**
     root, left subtree, right subtree.
     Note that a pre-order traversal gives the Euler tour for this tree when
     Euler start node is the root.
     <pre>
       e.g.
                            7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12

     visits: 7, 3, 1, 0, 2,5,4,6,11,9,8,10,13,12
     </pre>
     @param root root of tree
     */
    public void preorderRecursive(BinaryTreeNode<Integer> root) {
        if (root != null) {
            System.out.println(root.getData());
            preorderRecursive(root.getLeft());
            preorderRecursive(root.getRight());
        }
    }

    /**
     * left subtree, root, right subtree
     <pre>
                           7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12

     visits: 0,1,2,3,4,5,6,7,8,9,10,11,12,13
     </pre>
     @param root tree root
     */
    public void inorderRecursive(BinaryTreeNode<Integer> root) {
        if (root != null) {
            inorderRecursive(root.getLeft());
            System.out.printf("%d, ", root.getData());
            inorderRecursive(root.getRight());
        }
    }

    /**
     * left subtree, right subtree, root
     <pre>
       e.g.
                           7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12

     visits: 0, 2, 1, 4, 6, 5,3, 8, 10, 9, 12, 13, 11, 7

     </pre>
     @param root tree root
     */
    public void postorderRecursive(BinaryTreeNode<Integer> root) {
        if (root != null) {
            postorderRecursive(root.getLeft());
            postorderRecursive(root.getRight());
            System.out.println(root.getData());
        }
    }

    /**
     * root, left subtree, right subtree.

     Note that a pre-order traversal gives the Euler tour for this tree when
     Euler start node is the root.
     <pre>
                           7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12

     visits: 7, 3, 1, 0, 2,5,4,6,11,9,8,10,13,12

     </pre>
     @param node
     */
    public void preorderIterative(BinaryTreeNode<Integer> node) {
        Stack<BinaryTreeNode<Integer>> stack = new Stack<>();
        while (!stack.isEmpty() || node != null) {
            if (node != null) {
                System.out.printf("%d, ", node.getData());
                stack.push(node);
                node = node.getLeft();
            } else {
                node = stack.pop();
                node = node.getRight();
            }
        }
        System.out.println();
    }

    public int[] getPreorder(BinaryTreeNode<Integer> node) {
        List<Integer> vals = new ArrayList<>();
        Stack<BinaryTreeNode<Integer>> s = new Stack<>();
        while (!s.isEmpty() || node != null) {
            if (node != null) {
                vals.add(node.getData());
                s.push(node);
                node = node.getLeft();
            } else {
                node = s.pop();
                node = node.getRight();
            }
        }
        int[] out = new int[vals.size()];
        for (int i = 0; i < vals.size(); ++i) {
            out[i] = vals.get(i);
        }
        return out;
    }

    /**
     * left subtree, root, right subtree
     <pre>
     e.g.
                           7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12

     visits: 0,1,2,3,4,5,6,7,8,9,10,11,12,13
     </pre>
     @param node
     */
    public void inorderIterative(BinaryTreeNode<Integer> node) {
        Stack<BinaryTreeNode<Integer>> stack = new Stack<>();
        int c = 0;
        while (!stack.isEmpty() || node != null) {
            c++;
            if (node != null) {
                stack.push(node);
                node = node.getLeft();
            } else {
                node = stack.pop();
                System.out.printf("%d, ", node.getData());
                node = node.getRight();
            }
        }
        System.out.printf("  nCalls=%d", c);
        System.out.println();
    }

    static class Snapshot {
        final BinaryTreeNode<Integer> node;
        int stage;
        public Snapshot(BinaryTreeNode<Integer> node, int stage) {
            this.node = node;
            this.stage = stage;
        }
    }
    /* left subtree, root, right subtree
                           7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12

    visits: 0,1,2,3,4,5,6,7,8,9,10,11,12,13
         */

    /**
     *
     @param node
     */

    public void inorderIterative2(BinaryTreeNode<Integer> node) {
        // recursion to iteration:
        //    https://www.codeproject.com/Articles/418776/How-to-replace-recursive-functions-using-stack-and
        Stack<Snapshot> s = new Stack<Snapshot>();
        s.push(new Snapshot(node, 0));
        int c = 0;
        Snapshot currentSnapshot;
        while (!s.isEmpty()) {
            c++;
            currentSnapshot = s.pop();
            // 2 recursive function calls, so 3 stages
            if (currentSnapshot.node != null) {
                switch(currentSnapshot.stage) {
                    case 0: {
                        currentSnapshot.stage++;
                        s.push(currentSnapshot);
                        s.push(new Snapshot(currentSnapshot.node.getLeft(), 0));
                        break;
                    }
                    case 1: {
                       System.out.printf("%d, ", currentSnapshot.node.getData());
                        currentSnapshot.stage++;
                        s.push(currentSnapshot);
                        break;
                    }
                    case 2:{
                        s.push(new Snapshot(currentSnapshot.node.getRight(), 0));
                        break;
                    }
                }
            }
        }
        System.out.printf("  nCalls=%d", c);
        System.out.println();
    }

    /**
     * left subtree, right subtree, root
     <pre>
       e.g.
                           7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12

     visits: 0, 2, 1, 4, 6, 5,3, 8, 10, 9, 12, 13, 11, 7

     </pre>
     this is actually bottom-up post-order iterative while recursive
     is top-down.
     TODO:  read this and convert here:
         https://www.cs.odu.edu/~zeil/cs361/latest/Public/recursionConversion/index.html
     @param node
     */
    public void postorderIterative(BinaryTreeNode<Integer> node) {
        Stack<BinaryTreeNode<Integer>> stack = new Stack<>();
        Stack<BinaryTreeNode<Integer>> stack2 = new Stack<>();
        stack.push(node);
        while (!stack.isEmpty()) {
            node = stack.pop();
            stack2.push(node);
            if (node.getLeft() != null) {
                stack.push(node.getLeft());
            }
            if (node.getRight() != null) {
                stack.push(node.getRight());
            }
        }
        while (!stack2.isEmpty()) {
            node = stack2.pop();
            System.out.printf("%d, ", node.getData());
        }
        System.out.println();
    }

    /**
     * a.k.a. breadth first traversal
     <pre>
       e.g.
                           7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12

     visits: 7, 3, 11, 1, 5, 9, 13, 0, 2, 4, 6, 8, 10, 12

     </pre>
     @param node tree root
     */
    public void levelOrderIterative(BinaryTreeNode<Integer> node) {
        Queue<BinaryTreeNode<Integer>> queue = new ArrayDeque<>();
        while (node != null) {
            System.out.printf("%d, ", node.getData());
            if (node.getLeft() != null) {
                queue.add(node.getLeft());
            }
            if (node.getRight() != null) {
                queue.add(node.getRight());
            }
            node = queue.poll(); // returns null if empty
        }
        System.out.println();
    }

    /**
     *
     @param node
     */
    public static void printLevelOrder(HeapNode node) {

        Queue<HeapNode> queue = new ArrayDeque<>();
        Queue<Integer> nodeLevel = new ArrayDeque<>();
        Queue<Long> nodeParent = new ArrayDeque<>();
        Queue<Character> leftOrRight = new ArrayDeque<>();

        int level = 0;
        Character lOrR = '-';
        Long parent = -1L;

        while (node != null) {

            System.out.printf("L=%d, [%d, %s], %s of key %d\n",
                level, node.getKey(), node.getData() != null ?
                ((Integer)node.getData()).toString() : "-",
                lOrR, parent);

            ++level;

            if (node.getLeft() != null) {
                queue.add(node.getLeft());
                nodeLevel.add(level);
                nodeParent.add(node.getKey());
                leftOrRight.add('L');
            }
            if (node.getRight() != null) {
                queue.add(node.getRight());
                nodeLevel.add(level);
                nodeParent.add(node.getKey());
                leftOrRight.add('R');
            }
            //TODO: rearrange to check empty just once instead of at while loop too
            if (queue.isEmpty() || nodeLevel.isEmpty()) {
                break;
            }
            node = queue.poll(); // returns null if empty
            level = nodeLevel.poll();
            lOrR = leftOrRight.poll();
            parent = nodeParent.poll();
        }
    }


    /**
     * get the reverse level-order traversal of tree node.
     * implemented as post-order traversal but using a queue for the first
     * stack:
     * adapted from https://www.geeksforgeeks.org/inorder-tree-traversal-without-recursion/?ref=gcse
     <pre>
       e.g.
                            7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12

     visits: 12, 10, 8, 6, 4, 2, 0, 13, 9, 5, 1, 11, 3, 7
     </pre>
     @param node
     */
    public void reverseLevelOrderIterative(BinaryTreeNode<Integer> node) {
        Queue<BinaryTreeNode<Integer>> queue = new ArrayDeque<>();//FIFO
        Stack<BinaryTreeNode<Integer>> stack2 = new Stack<>();
        queue.add(node);
        while (!queue.isEmpty()) {
            node = queue.remove();//retrieves first in list
            stack2.push(node);
            if (node.getLeft() != null) {
                queue.add(node.getLeft());
            }
            if (node.getRight() != null) {
                queue.add(node.getRight());
            }
        }
        while (!stack2.isEmpty()) {
            node = stack2.pop();
            System.out.printf("%d, ", node.getData());
        }
        System.out.println();
    }

    /**
     * given a tree represented by node, return a doubly-linked list of nodes
     * visited in a reverse level-order traversal.
     *
     @param node the tree to be traversed from bottom up to this node using
     * reverse level-order traversal.
     * NOTE that any next and prev links in the tree are overwritten by the
     * DoublyLinkedList, so if those need to be preserved, give this method
     * node.copyTree() instead.
     @return a double-linked list of nodes in reverse level order traversal.
     */
    public DoublyLinkedList<NAryTreeNode> getReverseLevelOrderIterative2(NAryTreeNode node) {

        DoublyLinkedList<NAryTreeNode> out = new DoublyLinkedList<NAryTreeNode>();

        Queue<NAryTreeNode> queue = new ArrayDeque<>();//FIFO
        Stack<NAryTreeNode> stack2 = new Stack<>();
        queue.add(node);

        Set<NAryTreeNode> children;

        while (!queue.isEmpty()) {
            node = queue.remove();//retrieves first in list
            stack2.push(node);
            children = node.getChildren();
            if (children == null) {
                continue;
            }
            queue.addAll(children);
        }
        while (!stack2.isEmpty()) {
            node = stack2.pop();
            out.add(node);
        }
        return out;
    }

    /**
     * a.k.a. breadth first traversal
     <pre>
       e.g.
                            7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12

     visits: 7, 3, 11, 1, 5, 9, 13, 0, 2, 4, 6, 8, 10, 12
     </pre>
     @param node n-ary tree root
     @return
     */
    public DoublyLinkedList<NAryTreeNode> getLevelOrderIterative(NAryTreeNode node) {
        Queue<NAryTreeNode> queue = new ArrayDeque<>();
        DoublyLinkedList<NAryTreeNode> out = new DoublyLinkedList<NAryTreeNode>();
        Set<NAryTreeNode> children;
        while (node != null) {
            out.add(node);
            children = node.getChildren();
            if (children == null) {
                continue;
            }
            queue.addAll(children);

            node = queue.poll(); // returns null if empty
        }
        return out;
    }

    /**
     *
     @param node tree root
     @return a double-linked list of nodes in reverse level order traversal.
     */
    public List<BinaryTreeNode<Integer>> getReverseLevelOrderIterative(BinaryTreeNode<Integer> node) {
        List<BinaryTreeNode<Integer>> out = new ArrayList<>();
        Queue<BinaryTreeNode<Integer>> queue = new ArrayDeque<>();//FIFO
        Stack<BinaryTreeNode<Integer>> stack2 = new Stack<>();
        queue.add(node);
        while (!queue.isEmpty()) {
            node = queue.remove();//retrieves first in list
            stack2.push(node);
            if (node.getLeft() != null) {
                queue.add(node.getLeft());
            }
            if (node.getRight() != null) {
                queue.add(node.getRight());
            }
        }
        while (!stack2.isEmpty()) {
            node = stack2.pop();
            out.add(node);
        }
        return out;
    }

    /**
     * find the successor node to node in a binary search tree.
     <pre>
     following CLRS chap 12
     (Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms")
     </pre>
     * @param node
     * @return
     */
    public BinaryTreeNode<Integer> successor(BinaryTreeNode<Integer> node) {
        if (node.getRight() != null) {
            // left most node in right subtree
            return minimum(node.getRight());
        }

        // find lowest ancestor of node whose left child is an ancestor of node.
        BinaryTreeNode<Integer> parent = node.getParent();
        BinaryTreeNode<Integer> child = node;
        while (parent != null && child == parent.getRight()) {
            child = parent;
            parent = parent.getParent();
        }
        return parent;
    }

    /**
     find the predeccesor node to node in a binary search tree.
     <pre>
     following CLRS chap 12
     ( Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms")

     e.g.
     the predecessor of 5 is 4.
     if 4 were null, the predecessor of 5 is 3.
                            7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12
     </pre>
     * @param node
     * @return
     */
    public BinaryTreeNode<Integer> predeccesor(BinaryTreeNode<Integer> node) {
        if (node.getLeft() != null) {
            // max in left subtree
            return maximum(node.getLeft());
        }
        // find max in parent left subtree
        BinaryTreeNode<Integer> parent = node.getParent();
        BinaryTreeNode<Integer> child = node;
        while (parent != null && child == parent.getLeft()) {
            child = parent;
            parent = parent.getParent();
        }
        return parent;
    }

    /**
      find the minimum node under node in a binary search tree.
     <pre>
     following CLRS chap 12
     (Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms")

     e.g.
     the minimum node under 11 is 8

                            7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12
     </pre>
     * @param node
     * @return
     */
    public BinaryTreeNode<Integer> minimum(BinaryTreeNode<Integer> node) {
        while (node.getLeft() != null) {
            node = node.getLeft();
        }
        return node;
    }

    /**
      find the maximum node under node in a binary search tree.
     <pre>
     following CLRS chap 12
     (Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms")

     e.g.
     the maximum node under 3 is 6

                            7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12
     </pre>
     * @param node
     * @return
     */
    public BinaryTreeNode<Integer> maximum(BinaryTreeNode<Integer> node) {
        while (node.getRight() != null) {
            node = node.getRight();
        }
        return node;
    }
}
