package algorithms.trees;

import java.util.*;

/**
 * a binary tree with additional operations for search, successor, predeccesor, minimum, and maximum.
 * The non-recursive binary search tree operations here follow CLRS, chap 17.
 * @param <T> data type within the BinaryTreeNode
 * @param <S> the BinaryTreeNode or extension
 */
public class BinarySearchTree<T extends Comparable<T>, S extends BinaryTreeNode<T>>
        extends AbstractBinarySearchTree<T, S> {
    // NOTE" internally, uses convention for interpreting S.compareTo
    // as left is < 0 else right

    public BinarySearchTree(){}

    @Override
    public S insert(T data) {
        S node = newNode(data);
        S parent = _insert(root, node);
        if (parent == null) {
            root = node;
        }

        // update node n's
        updateN(node);

        // update ancestor N's
        while (parent != null) {
            updateN(parent);
            parent = (parent.parent == null) ? null : (S)parent.parent;
        }
        return node;
    }

    /**
     * method to construct a BinaryTree node or extension.  overrride as needed.
     * default is to construct S
     * @param data data for new node to hold
     * @return the new node
     */
    public S newNode(T data) {
        return (S) new BinaryTreeNode<T>(data);
    }

    /**
     * insert node into subtree top.  if result is that node has no parent, then top was null and null is returned,
     * else the parent node of the inserted node is returned.
     * @param top the subtree root node
     * @param node the node to insert into the subtree
     * @return parent of node inserted into tree
     */
    protected S _insert(S top, S node) {
        S x = top;
        S p = null;
        // descend to leaf node
        while (x != null) {
            p = x;
            if (node.data.compareTo(x.data) < 0) {
                x = (x.left == null) ? null : (S)x.left;
            } else {
                x = (x.right == null) ? null : (S)x.right;
            }
        }
        node.parent = p;
        if (p == null) {
            // parent is null so node becomes top for invoker
            return null;
        }
        // set node to be left or right for parent
        if (node.data.compareTo(p.data) < 0) {
            p.left = node;
        } else {
            p.right = node;
        }
        return p;
    }

    @Override
    protected S _search(S top, T data) {
        if (top == null) {
            return null;
        }
        S node = top;
        int comp;
        while (node != null) {
            comp = data.compareTo(node.data);
            if (comp == 0) {
                return node;
            } else if (comp < 0) {
                node = (node.left == null) ? null : (S)node.left;
            } else {
                node = (node.right == null) ? null : (S)node.right;
            }
        }
        return node;
    }

    protected void transplant(BinaryTreeNode<T> oldNode, BinaryTreeNode<T> replNode) {
        if (oldNode == null) return;
        BinaryTreeNode<T> p = oldNode.parent;
        if (p == null) {
            root = (replNode == null) ? null : (S)replNode;
        } else if (p.left != null && p.left.equals(oldNode)) {
            p.left = replNode;
        } else if (p.right != null && p.right.equals(oldNode)) {
            p.right = replNode;
        }
        if (replNode != null) {
            replNode.parent = p;
        }
    }

    /**
     * delete node from this tree.
     * note that the in-order traversal print remains the same, but is missing the deleted node.
     * @param node node to delete
     * @return true always
     */
    public boolean delete(BinaryTreeNode<T> node) {

        // delete when node has no children
        if (node.left == null && node.right == null) {
            transplant(node, node.right);
            updateN(node.parent);
        } else if (node.left == null) {
            // delete when node has 1 child
            transplant(node, node.right);
            updateN(node.parent);
        } else if (node.right == null) {
            // delete when node has 1 child
            transplant(node, node.left);
            updateN(node.parent);
        } else {
            // delete when node has 2 children, choose the min of node
            // to replace it

            // y as a successor of node.   min doesn't have a left child
            BinaryTreeNode<T> y = minimum((node.right == null) ? null : (S) node.right);
            if (!y.equals(node.right)) {// node.right is not null, so y is not null
                transplant(y, y.right);
                y.right = node.right;
                y.right.parent = y;
            }
            transplant(node, y);
            y.left = node.left;
            y.left.parent = y;

            // update the n fields:
            if (y != null) {
                if (y.left != null) {
                    updateN(y.left);
                }
                if (y.right != null) {
                    updateN(y.right);
                }
                updateN(y);
                BinaryTreeNode<T> parent = y.parent;
                while (parent != null) {
                    updateN(parent);
                    parent = parent.parent;
                }
            }
        }

        return true;
    }
}
