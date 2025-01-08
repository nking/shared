package algorithms.trees;

/**
 *
 * @author nichole
 */
public class BinaryTreeNode<T extends Comparable<T>> {

    protected T data;
    protected BinaryTreeNode<T> left;
    protected BinaryTreeNode<T> right;
    protected BinaryTreeNode<T> parent;

    /**
     * number of nodes in this subtree, that is, self + descendants.
     */
    protected long n = 0;
    /**
     @param data
     */
    public BinaryTreeNode(T data) {
        this.data = data;
        n = 1;
    }

    /**
     @return the data
     */
    public T getData() {
        return data;
    }

    /**
     @return the left
     */
    public BinaryTreeNode<T> getLeft() {
        return left;
    }

    /**
     @return the right
     */
    public BinaryTreeNode<T> getRight() {
        return right;
    }

    /**
     @return the parent
     */
    public BinaryTreeNode<T> getParent() {
        return parent;
    }

    /**
     @param data the data to set
     */
    public void setData(T data) {
        this.data = data;
    }

    /**
     @param left the left to set
     */
    public void setLeft(BinaryTreeNode<T> left) {
        this.left = left;
    }

    /**
     @param right the right to set
     */
    public void setRight(BinaryTreeNode<T> right) {
        this.right = right;
    }

    /**
     @param parent the parent to set
     */
    public void setParent(BinaryTreeNode<T> parent) {
        this.parent = parent;
    }

}
