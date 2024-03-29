package algorithms;

import algorithms.heapsAndPQs.HeapNode;

/**
 * a doubly linked list with a sentinel in between the last and first item.
 * created specifically to hold Fibonacci Heap nodes.
 *
 * <pre>
 * Runtime complexity:
 * insert is O(1)
 * search is O(N)
 * delete if given node is O(1), else O(N)
 *</pre>
 * @author nichole
 */
public class DoubleLinkedCircularList {

    //appears between head and tail
    private final HeapNode sentinel;
    
    /**
     */
    public final static long sentinelKey = Long.MIN_VALUE;
    
    /**
     */
    public final static long noValue = Long.MIN_VALUE + 1;
    
    /**
     */
    public final static long minValue = Long.MIN_VALUE + 2;

    /**
     *
     */
    protected long number = 0;
    
    /**
     */
    public DoubleLinkedCircularList() {
        sentinel = new HeapNode(sentinelKey);
        resetSentinel();
    }

    /**
     @return 
     */
    public HeapNode getSentinel() {
        return sentinel;
    }
    
    /**
     */
    public final void resetSentinel() {
        this.sentinel.setLeft(sentinel);
        this.sentinel.setRight(sentinel);
    }
    
    /**
     */
    public void resetNumber() {
        this.number = 0;
    }

    /**
    * insert new key into circular doubly linked list,
    * runtime is O(1),
    * 
    * Example:
    *
        sentinel to 2nd inserted to 1st inserted to [ sentinel.right ]
        sentinel from 2nd inserted from 1st inserted from [ sentinel.right ]
    *
    * <pre>
    * subsequent traversal by FIFO should use :
    *    sentinel and proceed left n=number of items
    * 
    * subsequent traversal by LIFO should use :
    *    sentinel and proceed right n=number of times
    * 
    * runtime complexity is O(1).
    * </pre>
    @param node node to insert
    @return inserted child node instance
    */
    public HeapNode insert(HeapNode node) {
        if (node == null) {
            throw new IllegalArgumentException("node cannot be null");
        }
        if (node.getKey() == noValue) {
            throw new IllegalArgumentException("node must have key set");
        }
        
        // nodes are inserted to the right of the sentinel.
                
        HeapNode rightOfSentinel = sentinel.getRight();
        
        node.setRight(rightOfSentinel);
        rightOfSentinel.setLeft(node);
        sentinel.setRight(node);
        node.setLeft(sentinel);
        number++;
        
        return node;
    }

    /**
     * remove first found node with key.
     * runtime complexity is O(n).
     @param key key of node to remove
     @return true if key was found and removed
     */
    public boolean remove(long key) {
        HeapNode cn = search(key);
        remove(cn);
        return (cn != null);
    }
    
    /**
     * remove node by connecting it's relationships to one another and removing
     * self.
     * runtime complexity is O(1).
     @param node  node to remove from heap
     */
    public void remove(HeapNode node) {
    	if (node == null) {
    		return;
    	}
        HeapNode right = node.getRight();
        HeapNode left = node.getLeft();
        right.setLeft(left);
        left.setRight(right);
        
        // reset node's ends to a sentinel.  the user's of the class use that logic.
        node.setRight(new HeapNode(sentinelKey));
        node.getRight().setRight(node);
        node.getRight().setLeft(node);
        node.setLeft(node.getRight());
                
        number--;
    }
    
    /**
     * insert insertNode to a place after existingNode.  The method does not
     * preserve left right relationships of insertNode, but preserves those of 
     * existingNode.  It expects that existingNode is part of this instance's
     * members and updates the number of items, for later use in traversals.
     * <pre>
     * Internally the insertNode is to the left of existingNode using the 
     * convention of this class.
     * 
     * subsequent traversal by FIFO should use :
    *    sentinel and proceed left n=number of items
    * 
    * subsequent traversal by LIFO should use :
    *    sentinel and proceed right n=number of times
     * 
     * runtime complexity is O(1).
     * </pre>
     @param existingNode node already present in heap
     @param insertNode the new node to be inserted into the heap
     */
    public void insertAfter(HeapNode existingNode, HeapNode insertNode) {
        
        if (insertNode.getLeft() != null || insertNode.getRight() != null) {
            throw new IllegalArgumentException("insertNode's existing left or "
                + " right are written over, so remove those before using this "
                + "method for clearer correct use");
        }
                
        HeapNode left = existingNode.getLeft();
        
        left.setRight(insertNode);
        insertNode.setLeft(left);
        existingNode.setLeft(insertNode);
        insertNode.setRight(existingNode);
        
        number++;
    }
    
    /**
     @return
     */
    public long getNumberOfNodes() {
        return number;
    }
   
    /**
     * runtime complexity is up to O(n), so if this method is often used,
     * should choose another data structure for the logic.
     * 
     @param key key of node to search for in heap
     @return the node having key 
     */
    public HeapNode search(long key) {
        
        HeapNode cn = sentinel.getRight();
        
        while ((cn.getKey() != sentinel.getKey()) && (cn.getKey() != key)) {
            cn = cn.getRight();
        }
        
        return ((cn.getKey() != noValue) &&  (cn.getKey() != sentinel.getKey()))
            ? cn : null;
    }
}
