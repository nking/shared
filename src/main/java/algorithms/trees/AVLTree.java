package algorithms.trees;

/**
 * a self balancing binary search tree which maintains a height difference no
 * greater than 1 between the 2 subtrees of a node.
 * The height is the largest distance from a node to a descendant leaf.
 * The implementation follows MIT opencourseware class 6.006 Unit 2.
 * @param <T> comparable parameter type for node data
 * @param <S> base node class that must be used or extended
 */
public class AVLTree<T extends Comparable<T>, S extends AVLTreeNode<T>>
        extends BinarySearchRecursiveTree<T, S>{

    // NOTE" internally, uses convention for interpreting S.compareTo
    // as left is < 0 else right

    public AVLTree() {
        super();
    }

    @Override
    public S newNode(T data) {
        return (S) new AVLTreeNode<T>(data);
    }

    protected void updateHeight(AVLTreeNode<T> node) {
        node.height = Math.max(getHeight(node.left), getHeight(node.right)) + 1;
    }

    protected int getHeight(BinaryTreeNode<T> node) {
        if (node == null) {
            return -1;
        }
        return ((AVLTreeNode<T>)node).height;
    }

    protected void rebalance(AVLTreeNode<T> node) {
        //after each operation, updateN for node.parent and the node
        while (node != null) {
            updateHeight(node);
            if (getHeight(node.left) >= 2 + getHeight(node.right)) {
                // node.left is not null
                if (getHeight(node.left.left) >= getHeight(node.left.right)) {
                    //rotation handles root assignment too when node2.parent=null
                    BinaryTreeNode<T> node2 = rotateRight(node);
                } else {
                    BinaryTreeNode<T> node2 = rotateLeft((AVLTreeNode<T>)node.left);
                    updateNAndAncestors(node);
                    updateHeight(node);
                    BinaryTreeNode<T> node3 = rotateRight(node);
                }
            } else if (getHeight(node.right) >= 2 + getHeight(node.left)) {
                // node.right is not null
                if (getHeight(node.right.right) >= getHeight(node.right.left)) {
                    BinaryTreeNode<T> node2 = rotateLeft(node);
                } else {
                    BinaryTreeNode<T> node2 = rotateRight((AVLTreeNode<T>) node.right);
                    updateNAndAncestors(node);
                    BinaryTreeNode<T> node3 = rotateRight(node);
                }
            }
            updateNAndAncestors(node);
            node = (node.parent == null) ? null : (AVLTreeNode<T>) node.parent;
        }
    }

    /**
     * insert data and return the node.  Note that every delete operation modifies the tree
     * in non-recursive version, so do not rely on the state of any node retained after
     * a bst delete.
     * @param data data for new node to hold
     * @return the inserted node
     */
    @Override
    public S insert(T data) {
        AVLTreeNode<T> node = super.insert(data);
        rebalance(node);
        return (S) node;
    }

    /**
     * delete node from this tree.
     * Also note that delete modifies internal state of 2 nodes and so any references
     * obtained during insert may have mutated.
     * @param node node to delete
     * @return true always
     */
    @Override
    public boolean delete(BinaryTreeNode<T> node) {
        AVLTreeNode<T> p = (node.parent == null) ? null : (AVLTreeNode<T>) node.parent;
        boolean s = super.delete(node);
        rebalance(p);
        return true;
    }

    /**
     * delete tree a single node having value data.
     * @param data value to search for in node to delete
     * @return true always
     */
    @Override
    public boolean delete(T data) {
        AVLTreeNode<T> node = search(data);
        return delete(node);
    }

    protected BinaryTreeNode<T> rotateRight(BinaryTreeNode<T> node) {
        BinaryTreeNode<T> node2 = super.rotateRight(node);
        updateHeight((AVLTreeNode<T>) node);
        updateHeight((AVLTreeNode<T>) node2);
        return node2;
    }

    protected BinaryTreeNode<T> rotateLeft(BinaryTreeNode<T> node) {
        BinaryTreeNode<T> node2 = super.rotateLeft(node);
        updateHeight((AVLTreeNode<T>) node);
        updateHeight((AVLTreeNode<T>) node2);
        return node2;
    }
}
