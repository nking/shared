package algorithms;

import algorithms.DoublyLinkedList.DoublyLinkedNode;

/**
 * a node for DoublyLinkedList holding data for the vertex index.
 * @author nichole
 */
public class VertexNode extends DoublyLinkedNode {
    
    /**
     *
     */
    public final int vertex;
    
    /**
     *
     @param v
     */
    public VertexNode(int v) {
        vertex = v;
    }
}
