package algorithms;

import algorithms.DoublyLinkedList.DoublyLinkedNode;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DoublyLinkedListTest extends TestCase {
    
    public DoublyLinkedListTest(String testName) {
        super(testName);
    }
    
    public static class TNode extends DoublyLinkedNode {
        protected int data;
        public TNode(int data) {
            this.data = data;
        }
    }

    /**
     * Test of add method, of class DoublyLinkedList.
     */
    public void testAdd() {
        int n = 10;
        int i;
        TNode[] nodes = new TNode[n];
        for (i = 0; i < n; ++i) {
            nodes[i] = new TNode(i);
        }
        DoublyLinkedList<TNode> linkedList = new DoublyLinkedList<TNode>();
        for (i = 0; i < n; ++i) {
            if ((i&1)==1) {
                linkedList.add(nodes[i]);
            } else {
                linkedList.addLast(nodes[i]);
            }
        }
        for (i = 0; i < n; ++i) {
            assertTrue(linkedList.contains(nodes[i]));
        }
        for (i = 0; i < n; ++i) {
            TNode t = linkedList.removeFirst();
            assertTrue(t.equals(nodes[i]));
            assertEquals(n-i-1, linkedList.size());
        }
        assertTrue(linkedList.isEmpty());
        
        // ===================
        //9,8,7,6,5,4,3,2,1,0
        for (i = 0; i < n; ++i) {
            linkedList.addFirst(nodes[i]);
        }
        assertTrue(linkedList.peekFirst().equals(nodes[9]));
        assertTrue(linkedList.peekLast().equals(nodes[0]));
        for (i = 0; i < n; ++i) {
            TNode t = linkedList.removeFirst();
            assertTrue(t.equals(nodes[n-i-1]));
            assertEquals(n-i-1, linkedList.size());
        }
        assertTrue(linkedList.isEmpty());
        
        linkedList.add(nodes[2]);
        assertTrue(linkedList.contains(nodes[2]));
        assertFalse(linkedList.contains(nodes[0]));
        linkedList.unlink(nodes[2]);
        assertFalse(linkedList.contains(nodes[2]));
    }

}
