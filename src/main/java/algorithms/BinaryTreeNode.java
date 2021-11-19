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

    public BinaryTreeNode(int data) {
        this.data = data;
    }

    public String toString() {
        return Integer.toString(getData());
    }

    /**
     * @return the data
     */
    public int getData() {
        return data;
    }

    /**
     * @return the left
     */
    public BinaryTreeNode getLeft() {
        return left;
    }

    /**
     * @return the right
     */
    public BinaryTreeNode getRight() {
        return right;
    }

    /**
     * @return the parent
     */
    public BinaryTreeNode getParent() {
        return parent;
    }

    /**
     * @param data the data to set
     */
    public void setData(int data) {
        this.data = data;
    }

    /**
     * @param left the left to set
     */
    public void setLeft(BinaryTreeNode left) {
        this.left = left;
    }

    /**
     * @param right the right to set
     */
    public void setRight(BinaryTreeNode right) {
        this.right = right;
    }

    /**
     * @param parent the parent to set
     */
    public void setParent(BinaryTreeNode parent) {
        this.parent = parent;
    }
}
