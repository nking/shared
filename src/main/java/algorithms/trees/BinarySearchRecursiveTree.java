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
public class BinarySearchRecursiveTree<T extends Comparable<T>, S extends BinaryTreeNode<T>>
        extends AbstractBinarySearchTree<T, S> {
        //<T extends Comparable<T>, S extends BinaryTreeNode<T>>
//extends AbstractBinarySearchTree<T, S> {

    // NOTE" internally, uses convention for interpreting S.compareTo
    // as left is < 0 else right

    public BinarySearchRecursiveTree(){}

    /**
     * method to construct a BinaryTree node or extension.  overrride as needed.
     * default is to construct S
     * @param data
     * @return
     */
    public S newNode(T data) {
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
     * insert data and return the node.  Note that every delete operation modifies the tree
     * in non-recursive version, so do not rely on the state of any node retained after
     * a bst delete.
     * @param data
     * @return
     */
    public S insert(T data) {
        // method implemented for sake of method documentation which is a warning.
        return super.insert(data);
    }


}
