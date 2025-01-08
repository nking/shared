package algorithms.trees;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

public abstract class AbstractBinarySearchTree
        <T extends Comparable<T>, S extends BinaryTreeNode<T>> implements IBinarySearchTree<T,S>{

    protected S root = null;

    protected abstract S _search(S top, T data);
    public abstract S newNode(T data);
    protected abstract S _insert(S top, S node);
    public abstract boolean delete(BinaryTreeNode<T> node);

    public S insert(T data) {
        S node = newNode(data);
        root = _insert(root, node);
        return node;
    }

    public S search(T data) {
        return _search(root, data);
    }

    public boolean contains(T data) {
        return (_search(root, data) != null);
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
        BinaryTreeNode<T> node2 = node.left; //*  y=node2
        BinaryTreeNode<T> p = node.parent;
        node.left = node2.right; //*
        node2.right = node;  //*

        if (node.left != null) {
            node.left.parent = node;
        }

        node.parent = node2;
        node2.parent = p;
        if (p == null) {
            root = (S) node2;
        } else {
            if (p.left.equals(node)) {
                p.left = node2;
            } else {
                p.right = node2;
            }
        }

        return node2;
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
     D   E
     </pre>
     @param node
     @return
     */
    protected BinaryTreeNode<T> rotateLeft(BinaryTreeNode<T> node) {
        BinaryTreeNode<T> node2 = node.right; //*
        BinaryTreeNode<T> p = node.parent;
        node.right = node2.left; //*
        node2.left = node;  //*
        if (node.right != null) {
            node.right.parent = node;
        }

        node.parent = node2;
        node2.parent = p;

        if (p == null) {
            root = (S) node2;
        } else {
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

}
