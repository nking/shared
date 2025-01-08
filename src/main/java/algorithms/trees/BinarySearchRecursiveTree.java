package algorithms.trees;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

/**
 * a binary tree with additional operations for search, successor, predeccesor, minimum, and maximum.
 *
 * The recursive binary search tree is easier to memorize, and is used as a basis for a balanced binary search
 * tree AVL implementation.
 * @param <T> data type within the BinaryTreeNode
 * @param <S> the BinaryTreeNode or extension
 */
public class BinarySearchRecursiveTree<T extends Comparable<T>, S extends BinaryTreeNode<T>> {

    protected S root = null;

    public BinarySearchRecursiveTree(){}

    /**
     * insert data and return the node.  Note that every delete operation modifies the tree
     * in non-recursive version, so do not rely on the state of any node retained after
     * a bst delete.
     * @param data
     * @return
     */
    public S insert(T data) {
        S node = newNode(data);
        root = _insert(root, node);
        return node;
    }

    protected void updateN(BinaryTreeNode<T> node) {
        if (node == null) return;
        node.n = 1;
        if (node.left != null) {
            node.n += node.left.n;
        }
        if (node.right != null) {
            node.n += node.right.n;
        }
    }

    protected void updateNAndAncestors(BinaryTreeNode<T> node) {
        if (node == null) return;
        updateN(node);
        BinaryTreeNode<T> parent = node.parent;
        while (parent != null) {
            updateN(parent);
            parent = parent.parent;
        }
    }

    /**
     * method to construct a BinaryTree node or extension.  overrride as needed.
     * default is to construct S
     * @param data
     * @return
     */
    protected S newNode(T data) {
        return (S) new BinaryTreeNode<T>(data);
    }

    /**
     * insert node into subtree top.  if result is that node has no parent, then top was null and null is returned,
     * else the parent node of the inserted node is returned.
     * Note that every delete operation modifies the tree
     * in non-recursive version, so do not rely on the state of any node retained after
     * @param top the subtree root node
     * @param node the node to insert into the subtree
     * @return parent of node inserted into tree
     * @param <T> node data type that extends Comparable
     */
    protected S _insert(S top, S node) {
        if (top == null) {
            return node;
        }
        if (node == null) {
            return null;
        }
        int comp = node.data.compareTo(top.data);
        if (comp < 0) {
            top.left = _insert((top.left == null) ? null : (S)top.left, node);
            top.left.parent = top;
            updateN((top.left == null) ? null : (S)top.left);
        } else {
            top.right = _insert((top.right == null) ? null : (S)top.right, node);
            top.right.parent = top;
            updateN((top.right == null) ? null : (S)top.right);
        }

        updateN(top);

        return top;
    }

    public boolean contains(T data) {
        return (search(data) != null);
    }

    protected S search(T data) {
        return _search(root, data);
    }

    protected S _search(S top, T data) {
        if (top == null) {
            return null;
        }
        int comp = data.compareTo(top.data);
        if (comp == 0) {
            return top;
        } else if (comp < 0) {
            return _search((top.left == null) ? null : (S)top.left, data);
        } else {
            return _search((top.right == null) ? null : (S)top.right, data);
        }
    }

    /**
     * delete tree a single node having value data.
     * @param node
     * @return
     */
    public boolean delete(T data) {
        S node = search(data);
        return delete(node);
    }

    /**
     * delete node from this tree.
     * note that the in-order traversal print remains the same, but is missing the deleted node.
     * Also note that delete modifies internal state of 2 nodes and so any references
     * obtained during insert may have mutated.
     * @param node
     * @return
     */
    public boolean delete(BinaryTreeNode<T> node) {
        // delete when node has no children
        if (node.left == null && node.right == null) {
            if (node.parent == null) {
                root = null;
            } else if (node.parent.left != null && node.parent.left.equals(node)) {
                node.parent.left = null;
                updateNAndAncestors(node.parent);
            } else {
                assert(node.parent.right.equals(node));
                node.parent.right = null;
                updateNAndAncestors(node.parent);
            }
        } else if (node.left == null || node.right == null) {
            S subTree = (node.left == null) ? (S)node.right : (S)node.left;
            // delete when node has 1 child
            if (node.parent == null) {
                root = subTree;
                if (subTree != null) {
                    subTree.parent = null;
                }
            } else if (node.parent.left != null && node.parent.left.equals(node)) {
                node.parent.left = subTree;
                subTree.parent = node.parent;
                updateNAndAncestors(node.parent.left);
            } else {
                assert(node.parent.right.equals(node));
                node.parent.right = subTree;
                subTree.parent = node.parent;
                updateNAndAncestors(node.parent.right);
            }
        } else {
            // delete when node has 2 children, choose the min of node
            // to replace it
            BinaryTreeNode<T> min = minimum(node.right); // node.right is not null

            //  min is to the right of node and min.left = null

            // keep node in current position, but swap all of its information
            // with min.
            // then delete min.

            T tmp = node.data;
            node.data = min.data;
            min.data = tmp;

            BinaryTreeNode<T> tmp1 = min.parent;

            boolean deleted = delete(min);

            //updateNAndAncestors(tmp1);
        }

        return true;
    }

    /**
     * find minimum node for subtree top.
     * @param top
     * @return minimum node
     * @param <T> comparable data carried by BinaryTreeNode
     */
    public S minimum(BinaryTreeNode<T> top) {
        if (top == null) return null;
        // is left node
        BinaryTreeNode<T> node = top;
        while (node.left != null) {
            node = (node.left == null) ? null : node.left;
        }
        return (node == null) ? null : (S)node;
    }

    /**
     * find maximum node for subtree top.
     * @param top
     * @return minimum node
     * @param <T> comparable data carried by BinaryTreeNode
     */
    public S maximum(BinaryTreeNode<T> top) {
        if (top == null) return null;
        // is rightmost node
        BinaryTreeNode<T> node = top;
        while (node.right != null) {
            node = node.right;
        }
        return (node == null) ? null : (S)node;
    }

    /**
     * find successor node for node.
     * @param top
     * @return successor node
     * @param <T> comparable data carried by BinaryTreeNode
     */
    public S successor(BinaryTreeNode<T> node) {
        if (node == null) return null;
        /*
                0
           10        20
                  11     30
                        21
         */
        if (node.right != null) {
            return minimum(node.right);
        }
        // find lowest ancestor of node whose left child is an ancestor of node
        BinaryTreeNode<T> p = node.parent;
        while (p != null && p.right != null && p.right.equals(node)) {
            p = p.parent;
        }
        return (p == null) ? null : (S)p;
    }

    /**
     * find predecessor node for node.
     * @param top
     * @return predecessor node
     * @param <T> comparable data carried by BinaryTreeNode
     */
    public S predecessor(BinaryTreeNode<T> node) {
        if (node == null) return null;
        /*
                0
           10        20
                  11     30
                        21
         */
        if (node.left != null) {
            return maximum(node.right);
        }
        // find lowest ancestor of node whose right child is an ancestor of node
        BinaryTreeNode<T> p = node.parent;
        while (p != null && p.left != null && p.left.equals(node)) {
            p = p.parent;
        }
        return (p == null) ? null : (S)p;
    }

    /**
     * rotate node to the right
     <pre>
     e.g.
     given tree or subtree A:
              A
           B      C
        D    E
     rightRotate(A) creates:
              B
           D      A
                 E  C
     </pre>
     * @param node
     * @return
     */
    protected BinaryTreeNode<T> rotateRight(BinaryTreeNode<T> node) {
        BinaryTreeNode<T> tmp = node.left; //*
        BinaryTreeNode<T> p = node.parent;
        node.left = tmp.right; //*

        tmp.right = node;  //*
        node.parent = tmp;
        tmp.parent = p;

        if (node.left != null) {
            node.left.parent = node;
        }

        if (p != null) {
            if (p.left.equals(node)) {
                p.left = tmp;
            } else {
                p.right = tmp;
            }
        }
        return tmp;
    }

    /**
     rotate node to the left
     <pre>
     e.g.
     given tree or subtree B:
              B
           D      A
                 E  C
     leftRotate(B) creates:
              A
           B      C
        D    E
     </pre>
     @param node
     @return
     */
    protected BinaryTreeNode<T> rotateLeft(BinaryTreeNode<T> node) {
        BinaryTreeNode<T> node2 = node.right; //*
        BinaryTreeNode<T> p = node.parent;
        node.right = node2.left; //*

        node2.left = node;  //*
        node.parent = node2;
        node2.parent = p;

        if (node.right != null) {
            node.right.parent = node;
        }

        if (p != null) {
            if (p.left.equals(node)) {
                p.left = node2;
            } else {
                p.right = node2;
            }
        }
        return node2;
    }

    /**
     * return data from in-order traversal of tree.
     * @return
     */
    public List<T> inOrderTraversal() {
        List<T> out = new ArrayList<>();
        BinaryTreeNode<T> node = root;
        if (node == null) {
            return out;
        }
        Stack<BinaryTreeNode<T>> s = new Stack<>();
        while (!s.isEmpty() || node != null) {
            if (node != null) {
                s.push(node);
                node = node.left;
            } else {
                node = s.pop();
                out.add(node.data);
                node = node.right;
            }
        }
        return out;
    }
}
