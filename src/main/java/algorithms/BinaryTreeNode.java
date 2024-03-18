package algorithms;

import algorithms.DoublyLinkedList.DoublyLinkedNode;

/**
 *
 * @author nichole
 */
public class BinaryTreeNode extends DoublyLinkedNode {

    private int data;
    private BinaryTreeNode left;
    private BinaryTreeNode right;
    private BinaryTreeNode parent;

    /**
     @param data
     */
    public BinaryTreeNode(int data) {
        this.data = data;
    }

    public String toString() {
        return Integer.toString(getData());
    }

    /**
     @return the data
     */
    public int getData() {
        return data;
    }

    /**
     @return the left
     */
    public BinaryTreeNode getLeft() {
        return left;
    }

    /**
     @return the right
     */
    public BinaryTreeNode getRight() {
        return right;
    }

    /**
     @return the parent
     */
    public BinaryTreeNode getParent() {
        return parent;
    }

    /**
     @param data the data to set
     */
    public void setData(int data) {
        this.data = data;
    }

    /**
     @param left the left to set
     */
    public void setLeft(BinaryTreeNode left) {
        this.left = left;
    }

    /**
     @param right the right to set
     */
    public void setRight(BinaryTreeNode right) {
        this.right = right;
    }

    /**
     @param parent the parent to set
     */
    public void setParent(BinaryTreeNode parent) {
        this.parent = parent;
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
    public static BinaryTreeNode rotateRight(BinaryTreeNode node) {
        BinaryTreeNode node2 = node.left; //*
        BinaryTreeNode p = node.parent;
        node.left = node2.right; //*

        node2.right = node;  //*
        node.parent = node2;
        node2.parent = p;

        if (node.left != null) {
            node.left.parent = node;
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
    public static BinaryTreeNode rotateLeft(BinaryTreeNode node) {
        BinaryTreeNode node2 = node.right; //*
        BinaryTreeNode p = node.parent;
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
}
