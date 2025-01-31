package algorithms.trees;

import java.util.List;

/**
 * interface for a binary search tree and 2 order statistic methods
 * @param <T> the data type in a node
 * @param <S> the node type in the tree
 */
public interface IBinarySearchTree
        <T extends Comparable<T>, S extends BinaryTreeNode<T>> {

    /**
     * create and return a new node
     * @param data value to place in node
     * @return new node
     */
    S newNode(T data);

    /**
     * insert a new node with data into the tree and return the new node.
     * @param data value to be held by new node
     * @return new node
     */
    S insert(T data);

    /**
     * returns true if tree has a node with value data
     * @param data value to search for in nodes
     * @return true if tree has a node with data
     */
    boolean contains(T data);

    /**
     * search for a node having data
     * @param data value to search for
     * @return a node having data, elese null if no nodes in tree have data.
     */
    S search(T data);

    /**
     * delete node with value data
     * @param data value to search for in tree nodes
     * @return true always,
     * TODO: change to return false if not found
     */
    boolean delete(T data);

    /**
     * delete this node from the tree
     * @param node node to delete
     * @return true always
     * TODO: change to return false if not actively connected to tree
     */
    boolean delete(BinaryTreeNode<T> node);

    /**
     * return minimum valued node of subtree top
     * @param top top of subtree
     * @return minimum of top
     */
    S minimum(BinaryTreeNode<T> top);

    /**
     * return maximum valued node of subtree top
     * @param top top of subtree
     * @return maximum of subtree top
     */
    S maximum(BinaryTreeNode<T> top);

    /**
     * return successor of node
     * @param node tree node
     * @return successor of node else null if none
     */
    S successor(BinaryTreeNode<T> node);

    /**
     * return predecessor of node
     * @param node tree node
     * @return predecessor of node else null if none
     */
    S predecessor(BinaryTreeNode<T> node);

    /**
     * return the rank of a node with this data value
     * @param data value to search for
     * @return the rank of a node that has value data
     */
    long rank(T data);

    /**
     * return the rank of a node with this data value
     * @param t subtree to search for node with data
     * @param data value to search for
     * @return the rank of a node that has value data in subtree t
     */
    long rank(BinaryTreeNode<T> t, T data);

    /**
     * select the node with rank
     * @param rank 1-based rank for node to return from this ordered tree
     * @return node having rank, else null if none
     */
    S select(long rank);

    /**
     * select the node with rank from subtree t
     * @param t subtree root node
     * @param rank 1-based rank for node to return from ordered subtree t
     * @return node having rank in subtree t, else null if none
     */
    S select(BinaryTreeNode<T> t, long rank);

    /**
     * an in-order traversal of tree
     * @return nodes from in-order traversal
     */
    List<T> inOrderTraversal();
}
