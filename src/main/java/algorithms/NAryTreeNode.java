package algorithms;

import algorithms.DoublyLinkedList.DoublyLinkedNode;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * a node which has a number of children which is not restricted to 2.
 * @author nichole
 */
public class NAryTreeNode extends DoublyLinkedNode {

    private int data;
    private final Set<NAryTreeNode> children = new HashSet<NAryTreeNode>();
    private NAryTreeNode parent;

    /**
     *
     @param data
     */
    public NAryTreeNode(int data) {
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
    public Set<NAryTreeNode> getChildren() {
        return children;
    }

    /**
     @return the parent
     */
    public NAryTreeNode getParent() {
        return parent;
    }

    /**
     @param data the data to set
     */
    public void setData(int data) {
        this.data = data;
    }

    /**
     * add a child node to this node.  a side-effect that should be idempotent is
     * that the method also sets the child's parent to this instance.
     @param child a child node to set
     */
    public void addChild(NAryTreeNode child) {
        children.add(child);
        child.setParent(this);
    }

    /**
     @param parent the parent to set
     */
    public void setParent(NAryTreeNode parent) {
        this.parent = parent;
    }

    /**
     * given a root node, copy the tree and return it.  The attributes of 
     * NAryTreeNode copied are data, parent and children.  The next and
     * prev attributes are not copied currently, but could be added upon need.
     * The runtime complexity is O(|V|).
     @param root copy the tree from root
     @return the copied tree root
     */
    public static NAryTreeNode copyTree(NAryTreeNode root) {
        
        // start the copy from bottom up using a reverse level order traversal.
        /*
        example reverse level order traversal is 15, 14, 13, 12, 11, 10, 9, 8,...
         
                                *0
                        *1                2
                  3                              *4
            *5        *6                        7  8  *9
         10 11 12    13 14                             15
        */
        /* copying these properties of NAryTreeNode:
              data, parent, children
           not copying these properties:
              prev, next
        */        
        TIntObjectMap<NAryTreeNode> copiedNodesMap = new TIntObjectHashMap<NAryTreeNode>();
        
        TreeTraversal tt = new TreeTraversal();
        DoublyLinkedList<NAryTreeNode> revLevOrder = tt.getReverseLevelOrderIterative2(root);
       
        NAryTreeNode node = revLevOrder.peekFirst();
        NAryTreeNode c, cp;
        while (node != null) {
            c = copiedNodesMap.get(node.getData());
            if (c == null){
                c = new NAryTreeNode(node.getData());
                copiedNodesMap.put(c.getData(), c);
            }
            if (node.getParent() != null) {
                cp = copiedNodesMap.get(node.getParent().getData());
                if (cp == null) {
                    cp = new NAryTreeNode(node.getParent().getData());
                    copiedNodesMap.put(cp.getData(), cp);
                }
                cp.getChildren().add(c);
                c.setParent(cp);
            }
            
            node = (NAryTreeNode) node.next;
        }
        return copiedNodesMap.get(root.getData());
    }
}
