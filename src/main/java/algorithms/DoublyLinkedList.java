package algorithms;

import algorithms.DoublyLinkedList.DoublyLinkedNode;

/**
 * implementation of a doubly linked list.  note that it is not a circularly
 * linked list.  specifically, list.head.prev == null and list.tail.next == null.
 * @author nichole
 * @param <T> node type that is or extends DoublyLinkedNode
 */
public class DoublyLinkedList<T extends DoublyLinkedNode> {

    public static class DoublyLinkedNode {
        public DoublyLinkedNode next = null;
        public DoublyLinkedNode prev = null;
    }
        
    protected T head = null;
    protected T tail = null;
    protected int n = 0;
    
    /**
     * add to end of doubly-linked list.
     * runtime complexity is O(1).
     * @param node node to add
     */
    public void add(T node) {
        if (head == null) {
            assert(tail == null);
            head = node;
            head.prev = null;
            tail = node;
            tail.next = null;
        } else {
            assert(tail != null);
            T t = tail;
            assert(tail.next == null);
            t.next = node;
            node.prev = t;
            node.next = null;
            tail = node;
        }
        n++;
    }
    
    /**
     * insert at end of doubly-linked list.  this method simply invokes add(node).
     * runtime complexity is O(1).
     * @param node node to add to end of list
     */
    public void addLast(T node) {
        add(node);
    }
    
    /**
     * insert at beginning of doubly-linked list.
     * runtime complexity is O(1).
     * @param node node to add to beginning of list
     */
    public void addFirst(T node) {
        if (head == null) {
            assert(tail == null);
            head = node;
            node.prev = null;
            node.next = null;
            tail = node;
        } else {
            assert(tail != null);
            T h = head;
            assert(h.prev == null);
            node.prev = null;
            node.next = h;
            h.prev = node;
            head = node;
        }
        n++;
    }
    
    /**
     * unlink node from the doubly-linked list.
     * note, currently it is the user's responsibility to make sure the node is
     * in the list before invoking this as the method does not perform a 
     * contains check to keep the runtime complexity small.
     * runtime complexity is O(1).
     * @param node the instance of T already in this instance.
     */
    public void unlink(T node) {
        if (node == null) {
            return;
        }
        if (node.prev == null) {
            if (!node.equals(head)) {
                throw new IllegalStateException("node.prev == null, so it "
                + "should be equal to this.head but is not.   The state of node "
                + "and this instance are not consistent.");
            }
            head = (T) head.next;
            if (head == null) {
                tail = null;
            } else {
                head.prev = null;
            }
            // node is now unlinked from this instance
            node.next = null;
            node.prev = null;
        } else if (node.next == null) {
            if (!node.equals(tail)) {
                throw new IllegalStateException("node.next == null, so it "
                + "should be equal to this.tail but is not.   The state of node "
                + "and this instance are not consistent.");
            }
            tail = (T) tail.prev;
            // tail cannot be null because node.prev was not == null
            tail.next = null;
            // node is now unlinked from this instance
            node.next = null;
            node.prev = null;
        } else {
            T prv = (T) node.prev;
            T nxt = (T) node.next;
            prv.next = nxt;
            // nxt cannot be null because node.next was not == null
            nxt.prev = prv;
            // node is now unlinked from this instance
            node.next = null;
            node.prev = null;
        }
        n--;
    }
    
    /**
     * search for an identical node in this doubly-linked list instance.  Note that the runtime
     * complexity worse case is O(n) where n is the number of items in this
     * doubly-linked list instance.
     * @param node an instance of type T which extends DoublyLinkedNode.
     * @return node the node in this doubly-linked list in which .equals(node) returns
     * true, else null for no node found.
     */
    public T search(T node) {
        T current = head;
        while (current != null) {
            if (current.equals(node)) {
                return current;
            }
            current = (T) current.next;
        }
        return null;
    }
    
    /**
     * search for an identical node in this doubly-linked list instance.  Note that the runtime
     * complexity worse case is O(n) where n is the number of items in this
     * doubly-linked list instance.
     * @param node the node in this doubly-linked list in which .equals(node) returns
     * true, else null for no node found.
     * @return node to search for in list
     */
    public boolean contains(T node) {
        T found = search(node);
        return (found != null);
    }
    
    /**
     * unlink the last element from the list and return it.  
     * for use when want to treat this doubly linked list as LIFO, assuming
     * items were inserted using add().
     * @return last node if any
     */
    public T removeLast() {
        T node = tail;
        unlink(node);
        return node;
    }
    
    /**
     * unlink the first element from the list and return it.  
     * for use when want to treat this doubly linked list as FIFO, assuming
     * items were inserted using add().
     * @return first node that was in the list
     */
    public T removeFirst() {
        T node = head;
        unlink(node);
        return node;
    }
    
    /**
     * return the first element in the list without removing it from the list.
     * @return returns first item in the list (note: it is not a copy of the item, intentionally)
     */
    public T peekFirst() {
        return head;
    }
    
    /**
     * return the last element in the list without removing it from the list.
     * @return the last item in the list. (not: it is not a copy of the item, intentionally).
     */
    public T peekLast() {
        return tail;
    }
    
    
    public void clear() {
        head = null;
        tail = null;
        n = 0;
    }
    
    public int size() {
        return n;
    }
    
    public boolean isEmpty() {
        return (n == 0);
    }
}
