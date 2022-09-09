package algorithms;

import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class NAryTreeNodeTest extends TestCase {
    
    public NAryTreeNodeTest(String testName) {
        super(testName);
    }

    public void testModifiers() {
        /*
        tree example is from lecture slides of Principal lecturer: Dr Thomas Sauerwald
        Advanced Algorithms, University of Cambridge.
        VII. Approximation Algorithms: Covering Problems
           https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/vertexcover.pdf
           who reference Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms"
        
                                *0
                        *1                2
                  3                              *4
            *5        *6                        7  8  *9
         10 11 12    13 14                             15
        */
        int n = 16;
        NAryTreeNode[] nodes = new NAryTreeNode[n];
        int i;
        for (i = 0; i < n; ++i) {
            nodes[i] = new NAryTreeNode(i);
        }
        nodes[0].addChild(nodes[1]);
        nodes[1].addChild(nodes[3]);
        nodes[3].addChild(nodes[5]);
        nodes[3].addChild(nodes[6]);
        nodes[5].addChild(nodes[10]);
        nodes[5].addChild(nodes[11]);
        nodes[5].addChild(nodes[12]);
        nodes[6].addChild(nodes[13]);
        nodes[6].addChild(nodes[14]);
        
        nodes[0].addChild(nodes[2]);
        nodes[2].addChild(nodes[4]);
        nodes[4].addChild(nodes[7]);
        nodes[4].addChild(nodes[8]);
        nodes[4].addChild(nodes[9]);
        nodes[9].addChild(nodes[15]);
        
        for (i = 0; i < n; ++i) {
            assertEquals(i, nodes[i].getData());
        }
        assertTrue(nodes[9].getChildren().contains(nodes[15]));
        assertTrue(nodes[15].getParent().equals(nodes[9]));
        assertNull(nodes[0].getParent());
        assertTrue(nodes[14].getChildren().isEmpty());
        assertFalse(nodes[0].getChildren().isEmpty());
        
        NAryTreeNode copiedNode = NAryTreeNode.copyTree(nodes[0]);
        TreeTraversal tt = new TreeTraversal();
        DoublyLinkedList<NAryTreeNode> copiedLevOrdTraversal = 
            tt.getLevelOrderIterative(copiedNode);
        assertEquals(n, copiedLevOrdTraversal.size());
        NAryTreeNode current = copiedLevOrdTraversal.peekFirst();
        Set<NAryTreeNode> children;
        TIntSet copiedChildrenData;
        while (current != null) {
            if (current.getParent() != null) {
                assertEquals(nodes[current.getData()].getParent().getData(),
                    current.getParent().getData());
            }
            copiedChildrenData = getData(current.getChildren());
            
            children = nodes[current.getData()].getChildren();
            
            assertEquals(children.size(), copiedChildrenData.size());

            for (NAryTreeNode child : children) {
                assertTrue(copiedChildrenData.contains(child.getData()));
            }
            
            current = (NAryTreeNode) current.next;
        }
    }    

    private TIntSet getData(Set<NAryTreeNode> children) {
        TIntSet set = new TIntHashSet();
        for (NAryTreeNode child : children) {
            set.add(child.getData());
        }
        return set;
    }

}
