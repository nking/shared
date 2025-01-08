package algorithms.trees;

import java.util.List;

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

    List<T> inOrderTraversal();
}
