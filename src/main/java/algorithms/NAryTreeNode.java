package algorithms;

import algorithms.DoublyLinkedList.DoublyLinkedNode;
import java.util.HashSet;
import java.util.Set;

/**
 * a node which has a number of children which is not restricted to 2.
 * @author nichole
 */
public class NAryTreeNode extends DoublyLinkedNode {

    private int data;
    private Set<NAryTreeNode> children = new HashSet<NAryTreeNode>();
    private NAryTreeNode parent;

    public NAryTreeNode(int data) {
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
    public Set<NAryTreeNode> getChildren() {
        return children;
    }

    /**
     * @return the parent
     */
    public NAryTreeNode getParent() {
        return parent;
    }

    /**
     * @param data the data to set
     */
    public void setData(int data) {
        this.data = data;
    }

    /**
     * add a child node to this node.  a side-effect that should be idempotent is
     * that the method also sets the child's parent to this instance.
     * @param child a child node to set
     */
    public void addChild(NAryTreeNode child) {
        children.add(child);
        child.setParent(this);
    }

    /**
     * @param parent the parent to set
     */
    public void setParent(NAryTreeNode parent) {
        this.parent = parent;
    }

}
