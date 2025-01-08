package algorithms.trees;

public class AVLTreeNode<T extends Comparable<T>> extends BinaryTreeNode<T> {

    // height from max leaf node
    protected int height = -1;

    /**
     * @param data
     */
    public AVLTreeNode(T data) {
        super(data);
    }

}
