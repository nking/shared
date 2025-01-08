package algorithms.trees;

import java.util.*;

/**
 * a binary tree with additional operations for search, successor, predeccesor, minimum, and maximum.
 * The non-recursive binary search tree operations here follow CLRS, chap 17.
 * @param <T> data type within the BinaryTreeNode
 * @param <S> the BinaryTreeNode or extension
 */
public class BinarySearchTree<T extends Comparable<T>, S extends BinaryTreeNode<T>> {

    protected S root = null;

    public BinarySearchTree(){}

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
     * @param top the subtree root node
     * @param node the node to insert into the subtree
     * @return parent of node inserted into tree
     * @param <T> node data type that extends Comparable
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

    public boolean contains(T data) {
        return (search(data) != null);
    }

    protected S search(T data) {
        if (root == null) {
            return null;
        }
        S node = root;
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
     * @param node
     * @return
     */
    public boolean delete(S node) {

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
            node = node.left;
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
        return (node == null) ? null : (S) node;
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
            return minimum((node.right == null) ? null : (S)node.right);
        }
        // find lowest ancestor of node whose left child is an ancestor of node
        BinaryTreeNode<T> p = (node.parent == null) ? null : (S)node.parent;
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
            return maximum((node.right == null) ? null : (S)node.right);
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
    public S rotateRight(BinaryTreeNode<T> node) {
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
        return (tmp == null) ? null : (S) tmp;
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
    public S rotateLeft(BinaryTreeNode<T> node) {
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
        return (node2 == null) ? null : (S)node2;
    }

    /**
     * return data from in-order traversal of tree.
     * @return
     */
    public List<T> inOrderTraversal() {
        List<T> out = new ArrayList<>();
        S node = root;
        if (node == null) {
            return out;
        }
        Stack<S> s = new Stack<>();
        while (!s.isEmpty() || node != null) {
            if (node != null) {
                s.push(node);
                node = (node.left == null) ? null : (S)node.left;
            } else {
                node = s.pop();
                out.add(node.data);
                node = (node.right == null) ? null : (S)node.right;
            }
        }
        return out;
    }
}
