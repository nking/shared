package algorithms;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TreeTraversalTest extends TestCase {

    
    public void test2() {

        /*
                0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
         */
        BinaryTreeNode[] nodes = new BinaryTreeNode[14];
        int i;
        for (i = 0; i < nodes.length; ++i) {
            nodes[i] = new BinaryTreeNode(i);
        }
        nodes[0].setLeft(nodes[1]);
        nodes[0].setRight(nodes[2]);
        nodes[1].setLeft(nodes[3]);
        nodes[1].setRight(nodes[4]);
        nodes[2].setLeft(nodes[5]);
        nodes[2].setRight(nodes[6]);
        nodes[3].setLeft(nodes[7]);
        nodes[4].setRight(nodes[10]);
        nodes[5].setRight(nodes[8]);
        nodes[6].setRight(nodes[11]);
        nodes[8].setRight(nodes[9]);
        nodes[11].setLeft(nodes[12]);
        nodes[11].setRight(nodes[13]);
        
        nodes[1].setParent(nodes[0]);
        nodes[2].setParent(nodes[0]);
        nodes[3].setParent(nodes[1]);
        nodes[4].setParent(nodes[1]);
        nodes[7].setParent(nodes[3]);
        nodes[5].setParent(nodes[2]);
        nodes[6].setParent(nodes[2]);
        nodes[8].setParent(nodes[5]);
        nodes[9].setParent(nodes[8]);
        nodes[10].setParent(nodes[4]);
        nodes[11].setParent(nodes[6]);
        nodes[12].setParent(nodes[11]);
        nodes[13].setParent(nodes[11]);

        TreeTraversal btt = new TreeTraversal();
        
        System.out.println("pre-order:");
        btt.preorderRecursive(nodes[0]);
        btt.preorderIterative(nodes[0]);
        
        System.out.println("in-order:");
        btt.inorderRecursive(nodes[0]);
        btt.inorderIterative(nodes[0]);
        btt.inorderIterative2(nodes[0]);

        System.out.println("post-order:");
        btt.postorderRecursive(nodes[0]);
        btt.postorderIterative(nodes[0]);
        
        System.out.println("level-order:");
        btt.levelOrderIterative(nodes[0]);

        System.out.println("reverse-level-order:");
        btt.reverseLevelOrderIterative(nodes[0]);

    }
    
    public void testNAryReversedLevelOrder() {
        /*
        tree example is from lecture slides of Principal lecturer: Dr Thomas Sauerwald
        Advanced Algorithms, University of Cambridge.
        VII. Approximation Algorithms: Covering PRoblmes
           https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/vertexcover.pdf
           who reference Cormen et al. "Introduction to Algorithms"
        
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
        
        System.out.println("NAry-tree reversed level order traversal");
        TreeTraversal tt = new TreeTraversal();
        DoublyLinkedList<NAryTreeNode> rlo = tt.getReverseLevelOrderIterative2(nodes[0]);
        assertEquals(n, rlo.size());
        NAryTreeNode current = rlo.peekFirst();
        for (i = 0; i < rlo.size(); ++i) {
            System.out.printf("%d, ", current.getData());
            assertEquals(n-i-1, current.getData());
            current = (NAryTreeNode) current.next;
        }
        System.out.println();
    }

}
