package algorithms.trees;

import algorithms.DoublyLinkedList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TreeTraversalTest extends TestCase {

    @SuppressWarnings({"unchecked", "rawtypes"})
    public void test2() {

        /*
                0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
         */
        algorithms.trees.BinaryTreeNode<Integer>[] nodes = new algorithms.trees.BinaryTreeNode[14];
        int i;
        for (i = 0; i < nodes.length; ++i) {
            nodes[i] = new algorithms.trees.BinaryTreeNode<Integer>(i);
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

        algorithms.trees.TreeTraversal btt = new algorithms.trees.TreeTraversal();
        
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
        VII. Approximation Algorithms: Covering Problems
           https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/vertexcover.pdf
           who reference Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms"
        
                                *0
                        *1                2
                  3                              *4
            *5        *6                        7  8  *9
         10 11 12    13 14                             15
        visited: 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
        */
        int n = 16;
        algorithms.trees.NAryTreeNode[] nodes = new algorithms.trees.NAryTreeNode[n];
        int i;
        for (i = 0; i < n; ++i) {
            nodes[i] = new algorithms.trees.NAryTreeNode(i);
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
        
        //System.out.println("NAry-tree reversed level order traversal");
        algorithms.trees.TreeTraversal tt = new algorithms.trees.TreeTraversal();
        DoublyLinkedList<NAryTreeNode> rlo = tt.getReverseLevelOrderIterative2(nodes[0]);
        assertEquals(n, rlo.size());
        algorithms.trees.NAryTreeNode current = rlo.peekFirst();
        for (i = 0; i < rlo.size(); ++i) {
            //System.out.printf("%d, ", current.getData());
    //        assertEquals(n-i-1, current.getData());
            current = (algorithms.trees.NAryTreeNode) current.next;
        }
        //System.out.println();
    }


    public void testNAryLevelOrder() {
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
        
        visited: 0, 2, 1, 4, 3, 7, 8, 9, 5, 6, 15, 12, 11, 10, 14, 13,
        */
        TIntSet[] expectedLevels = new TIntSet[5];
        expectedLevels[0] = new TIntHashSet(new int[]{0});
        expectedLevels[1] = new TIntHashSet(new int[]{1, 2});
        expectedLevels[2] = new TIntHashSet(new int[]{3, 4});
        expectedLevels[3] = new TIntHashSet(new int[]{5,6, 7,8,9});
        expectedLevels[4] = new TIntHashSet(new int[]{10, 11, 12, 13, 14, 15});
        
        int n = 16;
        algorithms.trees.NAryTreeNode[] nodes = new algorithms.trees.NAryTreeNode[n];
        int i;
        for (i = 0; i < n; ++i) {
            nodes[i] = new algorithms.trees.NAryTreeNode(i);
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
        
        //System.out.println("NAry-tree level order traversal");
        algorithms.trees.TreeTraversal tt = new algorithms.trees.TreeTraversal();
        DoublyLinkedList<algorithms.trees.NAryTreeNode> rlo = tt.getLevelOrderIterative(nodes[0]);
        assertEquals(n, rlo.size());
        algorithms.trees.NAryTreeNode current = rlo.peekFirst();
        TIntSet expectedLevel = null;
        int j;
        for (i = 0; i < rlo.size(); ++i) {
            System.out.printf("%d, ", current.getData());
            for (j = 0; j < expectedLevels.length; ++j) {
                if (expectedLevels[j].isEmpty()) {
                    continue;
                }
                expectedLevel = expectedLevels[j];
                break;
            }
            assertNotNull(expectedLevel);
            assertTrue(expectedLevel.remove(current.getData()));
            current = (NAryTreeNode) current.next;
        }
        //System.out.println();
        for (j = 0; j < expectedLevels.length; ++j) {
            assertTrue(expectedLevels[j].isEmpty());
        }
    }

    /**
     test data for tree:
                            7
                3                        11
           1         5              9         13
        0   2     4   6               10    12
     * @return
     */

    public void testPredecessor() {
        algorithms.trees.TreeTraversal traversal = new algorithms.trees.TreeTraversal();
        algorithms.trees.BinaryTreeNode<Integer> root = getTree0();
        algorithms.trees.BinaryTreeNode<Integer> r;

        // get pred of node 9
        r = traversal.predeccesor(root.getRight().getLeft());
        assertEquals(7, r.data.intValue());
        // get pred of node 5
        r = traversal.predeccesor(root.getLeft().getRight());
        assertEquals(4, r.getData().intValue());

        // get pred of root
        r = traversal.predeccesor(root);
        assertEquals(6, r.getData().intValue());

        // pred for 0
        r = traversal.predeccesor(root.getLeft().getLeft().getLeft());
        assertNull(r);
    }

    /**
     test data for tree:
                            7
                3                        11
           1         5              9         13
        0   2     4   6               10    12
     * @return
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    public void testSuccessor() {
        algorithms.trees.TreeTraversal traversal = new TreeTraversal();
        algorithms.trees.BinaryTreeNode<Integer> root = getTree0();
        algorithms.trees.BinaryTreeNode<Integer> r;

        // get successor of node 9
        r = traversal.successor(root.getRight().getLeft());
        assertEquals(10, r.getData().intValue());

        // get successor of node 10
        r = traversal.successor(root.getRight().getLeft().getRight());
        assertEquals(11, r.getData().intValue());

        // remove node 6
        root.getLeft().getRight().setRight(null);
        // get successor of node 5
        r = traversal.successor(root.getLeft().getRight());
        assertEquals(7, r.getData().intValue());

        // get successor of root
        r = traversal.successor(root);
        assertEquals(9, r.getData().intValue());

        // get successor of 13
        r = traversal.successor(root.getRight().getRight());
        assertNull(r);
    }

    /**
                            7
                3                        11
           1         5              9         13
        0   2     4   6               10    12
     * @return
     */
    public algorithms.trees.BinaryTreeNode<Integer> getTree0() {
        algorithms.trees.BinaryTreeNode<Integer> root = new algorithms.trees.BinaryTreeNode<Integer>(7);
        addLeft(root, 3);
        addLeft(root.getLeft(),1);
        addRight(root.getLeft(),5);
        addLeft(root.getLeft().getLeft(),0);
        addRight(root.getLeft().getLeft(),2);
        addLeft(root.getLeft().getRight(),4);
        addRight(root.getLeft().getRight(),6);

        addRight(root, 11);
        addLeft(root.getRight(),9);
        addRight(root.getRight(),13);
        //addLeft(root.getRight().getLeft(), 8); removing x for unit test branch
        addRight(root.getRight().getLeft(), 10);
        addLeft(root.getRight().getRight(),12);
        return root;
    }
    protected algorithms.trees.BinaryTreeNode<Integer> addLeft(algorithms.trees.BinaryTreeNode<Integer> node, int val) {
        algorithms.trees.BinaryTreeNode<Integer> node2 = new algorithms.trees.BinaryTreeNode<Integer>(val);
        node.setLeft(node2);
        node2.setParent(node);
        return node;
    }
    protected algorithms.trees.BinaryTreeNode<Integer> addRight(algorithms.trees.BinaryTreeNode<Integer> node, int val) {
        algorithms.trees.BinaryTreeNode<Integer> node2 = new BinaryTreeNode<Integer>(val);
        node.setRight(node2);
        node2.setParent(node);
        return node;
    }
}
