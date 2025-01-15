package algorithms.trees;

import java.util.List;

/**
 * interface for a binary search tree and 2 order statistic methods
 * @param <T> the data type in a node
 * @param <S> the node type in the tree
 */
public interface IBinarySearchTree
        <T extends Comparable<T>, S extends BinaryTreeNode<T>> {

    S newNode(T data);

    S insert(T data);

    boolean contains(T data);

    S search(T data);

    boolean delete(T data);
    boolean delete(BinaryTreeNode<T> node);

    S minimum(BinaryTreeNode<T> top);

    S maximum(BinaryTreeNode<T> top);

    S successor(BinaryTreeNode<T> node);

    S predecessor(BinaryTreeNode<T> node);

    long rank(T data);

    long rank(BinaryTreeNode<T> t, T data);

    S select(long rank);

    S select(BinaryTreeNode<T> t, long rank);

    List<T> inOrderTraversal();
}
