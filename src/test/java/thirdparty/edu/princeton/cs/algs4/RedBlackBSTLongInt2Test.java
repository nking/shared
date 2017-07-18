package thirdparty.edu.princeton.cs.algs4;

import algorithms.misc.Misc0;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TLongArrayList;
import java.util.Random;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class RedBlackBSTLongInt2Test extends TestCase {
    
    public RedBlackBSTLongInt2Test(String testName) {
        super(testName);
    }
    
    public void testPutAndRotate() {
        
        int n = 5;
        RedBlackBSTLongInt2 bt = new RedBlackBSTLongInt2();
                
        TLongList nodes = new TLongArrayList(2*n);
        
        long[] kOutput = new long[2];
        
        for (int i = 0; i < n; ++i) {
            nodes.add(i);
            bt.put(i, i);
            assertTrue(bt.rootIsSet);
            assertTrue(bt.contains(i));
        }
        bt.printSmallTree(bt.root, n-1);
        
        /*
         *                      3
         *           root.left       root.right
         *               1                4
         *     left.left   left.right
         *         0           2
         * 
                  node=key=3 val=3 color=0 size=5 p= l=1 r=4
        [junit]   node=key=1 val=1 color=1 size=3 p=3 l=0 r=2
        [junit]   node=key=0 val=0 color=0 size=1 p=1 l= r=
        [junit]   node=key=2 val=2 color=0 size=1 p=1 l= r=
        [junit]   node=key=4 val=4 color=0 size=1 p=3 l= r=
        
        For a right rotate, left of it has to be red
        
                     RIGHT ROTATE(3)
        
                               1
                         0           3
                                   2   4
        node=key=1 val=1 color=0 size=5 p= l=0 r=3
        node=key=0 val=0 color=0 size=1 p=1 l= r=
        node=key=3 val=3 color=1 size=3 p=1 l=2 r=4
        node=key=2 val=2 color=0 size=1 p=3 l= r=
        node=key=4 val=4 color=0 size=1 p=3 l= r=
        
        For a left rotate, right of it has to be red
        
                     LEFT ROTATE(1)
         
                               3
                         1          4
                      0    2 
        
         node=key=3 val=3 color=0 size=5 p= l=1 r=4
         node=key=1 val=1 color=1 size=3 p=3 l=0 r=2
         node=key=0 val=0 color=0 size=1 p=1 l= r=
         node=key=2 val=2 color=0 size=1 p=1 l= r=
         node=key=4 val=4 color=0 size=1 p=3 l= r=
        */
                
        //bt.printPreOrderTraversal();
        
        //System.out.println("ROTATE-RIGHT(3)");
        
        long key = bt.rotateRight(3);
        
        bt.root = key;
        
        System.out.println("ROTATE-RIGHT");
        bt.printSmallTree(bt.root, n-1);
        
        //System.out.println("return key=" + key);
        
        //bt.printPreOrderTraversal();
        
        long[] treeNodes = bt.getPreOrderTraversalIterative(bt.root, 0);
        
        for (long node : treeNodes) {
            int c = (int)node;
            assertEquals(c, bt.nodeMap.getNodeValue(node));
            switch(c) {
                case 0:
                    assertEquals(0, bt.nodeMap.getNodeColor(node));
                    assertEquals(1, bt.nodeMap.getParent(node));
                    assertEquals(1, bt.nodeMap.getNodeSize(node));
                    assertFalse(bt.nodeMap.leftIsSet(node));
                    assertFalse(bt.nodeMap.rightIsSet(node));
                    break;
                case 1:
                    assertEquals(0, bt.nodeMap.getNodeColor(node));
                    assertFalse(bt.nodeMap.parentIsSet(node));
                    assertEquals(5, bt.nodeMap.getNodeSize(node));
                    assertEquals(0, bt.nodeMap.getLeft(node));
                    assertEquals(3, bt.nodeMap.getRight(node));
                    break;
                case 2:
                    assertEquals(0, bt.nodeMap.getNodeColor(node));
                    assertEquals(3, bt.nodeMap.getParent(node));
                    assertEquals(1, bt.nodeMap.getNodeSize(node));
                    assertFalse(bt.nodeMap.leftIsSet(node));
                    assertFalse(bt.nodeMap.rightIsSet(node));
                    break;
                case 3:
                    assertEquals(1, bt.nodeMap.getNodeColor(node));
                    assertEquals(1, bt.nodeMap.getParent(node));
                    assertEquals(3, bt.nodeMap.getNodeSize(node));
                    assertEquals(2, bt.nodeMap.getLeft(node));
                    assertEquals(4, bt.nodeMap.getRight(node));
                    break;
                case 4:
                    assertEquals(0, bt.nodeMap.getNodeColor(node));
                    assertEquals(3, bt.nodeMap.getParent(node));
                    assertEquals(1, bt.nodeMap.getNodeSize(node));
                    assertFalse(bt.nodeMap.leftIsSet(node));
                    assertFalse(bt.nodeMap.rightIsSet(node));
                    break;
                default:
                    break;
            }
        }
        
        key = bt.rotateLeft(1);
        
        bt.root = key;
        
        System.out.println("ROTATE-LEFT");
        bt.printSmallTree(bt.root, n-1);
        
        //System.out.println("return key=" + key);
        
        //bt.printPreOrderTraversal();
        
        treeNodes = bt.getPreOrderTraversalIterative(bt.root, 0);
        
        for (long node : treeNodes) {
            int c = (int)node;
            assertEquals(c, bt.nodeMap.getNodeValue(node));
            switch(c) {
                case 0:
                    assertEquals(0, bt.nodeMap.getNodeColor(node));
                    assertEquals(1, bt.nodeMap.getParent(node));
                    assertEquals(1, bt.nodeMap.getNodeSize(node));
                    assertFalse(bt.nodeMap.leftIsSet(node));
                    assertFalse(bt.nodeMap.rightIsSet(node));
                    break;
                case 1:
                    assertEquals(1, bt.nodeMap.getNodeColor(node));
                    assertEquals(3, bt.nodeMap.getParent(node));
                    assertEquals(3, bt.nodeMap.getNodeSize(node));
                    assertEquals(0, bt.nodeMap.getLeft(node));
                    assertEquals(2, bt.nodeMap.getRight(node));
                    break;
                case 2:
                    assertEquals(0, bt.nodeMap.getNodeColor(node));
                    assertEquals(1, bt.nodeMap.getParent(node));
                    assertEquals(1, bt.nodeMap.getNodeSize(node));
                    assertFalse(bt.nodeMap.leftIsSet(node));
                    assertFalse(bt.nodeMap.rightIsSet(node));
                    break;
                case 3:
                    assertEquals(0, bt.nodeMap.getNodeColor(node));
                    assertFalse(bt.nodeMap.parentIsSet(node));
                    assertEquals(5, bt.nodeMap.getNodeSize(node));
                    assertEquals(1, bt.nodeMap.getLeft(node));
                    assertEquals(4, bt.nodeMap.getRight(node));
                    break;
                case 4:
                    assertEquals(0, bt.nodeMap.getNodeColor(node));
                    assertEquals(3, bt.nodeMap.getParent(node));
                    assertEquals(1, bt.nodeMap.getNodeSize(node));
                    assertFalse(bt.nodeMap.leftIsSet(node));
                    assertFalse(bt.nodeMap.rightIsSet(node));
                    break;
                default:
                    break;
            }
        }
        
        // test flip nodes
        bt.flipColors(1);
        assertEquals(0, bt.nodeMap.getNodeColor(1));
        assertEquals(1, bt.nodeMap.getNodeColor(0));
        assertEquals(1, bt.nodeMap.getNodeColor(2));
        
        bt.flipColors(1);
        assertEquals(1, bt.nodeMap.getNodeColor(1));
        assertEquals(0, bt.nodeMap.getNodeColor(0));
        assertEquals(0, bt.nodeMap.getNodeColor(2));
    
        
    }
    
    public void testKeyOperations00() throws Exception {
    
        System.out.println("testKeyOperations");
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.currentTimeMillis();
        seed = 1500232861625L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        RedBlackBSTLongInt2 bt = new RedBlackBSTLongInt2();
        
        int n = 100;
        
        int count = 0;
        
        TLongList nodes = new TLongArrayList(2*n);
        
        long[] kOutput = new long[2];
        
        for (int i = 0; i < n/2; ++i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i += skip;
                continue;
            }
            nodes.add(i);
            bt.put(i, i);
            assertTrue(bt.rootIsSet);
            assertTrue(bt.contains(i));
            bt.min(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(nodes.get(0), kOutput[1]);
            
            bt.max(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(i, kOutput[1]);
            
            count++;
        }
        
        for (int i = (n - 1); i >= (n/2); --i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i -= skip;
                continue;
            }
            nodes.add(i);
            bt.put(i, i);
            assertTrue(bt.rootIsSet);
            assertTrue(bt.contains(i));
            bt.min(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(nodes.get(0), kOutput[1]);
            count++;
        }
        nodes.sort();
        int n2 = bt.size();
        assertEquals(count, n2);
        assertEquals(count, nodes.size());
        
        int[] vOutput = new int[2];
        
        for (int nIter = 0; nIter < 3; ++nIter){
            n2 = bt.size();
            assertEquals(n2, nodes.size());

            bt.printPreOrderTraversal();
            
            for (int i = 0; i < n2 - 1; ++i) {
                long idx = nodes.get(i);
                long foundIndex = bt.contains(idx) ? idx : -1;
                assertTrue(foundIndex > -1);

                long expected = nodes.get(i + 1);
                System.out.println("\n* " + idx + " expected next=" + expected);

                bt.higher(idx, kOutput);
                assertTrue(kOutput[0] != -1);
                long next = kOutput[1];
                System.out.println(idx + "   next=" + next);
                assertEquals(expected, next);

                assertTrue(bt.rootIsSet);
                
 //System.out.println("nIter=" + nIter + " i=" + i + " idx=" + idx);
 
                if (next > (idx + 1)) {
                    // test ceiling of idx+1
                    bt.ceiling(idx + 1, kOutput);
                    assertTrue(kOutput[0] != -1);
                    long ceil = kOutput[1];
                    assertEquals(expected, ceil);
                }
            }
          
     
            bt.higher(nodes.get(nodes.size() - 1), kOutput);
            assertTrue(kOutput[0] == -1);
            bt.printPreOrderTraversal();
            System.out.println("lower(" + nodes.get(0) + ")");
            bt.lower(nodes.get(0), kOutput);
            assertTrue(kOutput[0] == -1);
            
            for (int i = 1; i < n2; ++i) {
                long idx = nodes.get(i);
                long foundIndex = bt.contains(idx) ? idx : -1;
                assertTrue(foundIndex > -1);

                bt.get(idx, vOutput);
                assertTrue(vOutput[0] != -1);
                assertEquals(idx, vOutput[1]);
                
                assertTrue(bt.rootIsSet);
                
                long expected = nodes.get(i - 1);
                System.out.println("\n* " + idx + " expected prev=" + expected);

                bt.lower(idx, kOutput);
                assertTrue(kOutput[0] != -1);
                long prev = kOutput[1];
                //System.out.println(idx + "   prev=" + prev);
                assertEquals(expected, prev);

                if (prev < (idx - 1)) {
                    // test floor of idx-1
                    bt.floor(idx - 1, kOutput);
                    assertTrue(kOutput[0] != -1);
                    long floor = kOutput[1];
                    assertEquals(expected, floor);
                }
            }
        
            n2 = nodes.size();
            nodes.sort();
           
            assertEquals(nodes.get(nodes.size()/2), bt.select(nodes.size()/2));            
                  
            assertEquals(nodes.size(), bt.size());
            
            //System.out.println("before deleteMax bt.size=" + bt.size());
            
            long max = nodes.get(nodes.size() - 1);
            assertTrue(bt.contains(max));
            assertTrue(bt.rootIsSet);
            //System.out.println("will delete max=" + max + " from this tree:");
            //bt.printPreOrderTraversal();
            bt.deleteMax();
            //System.out.println("after deleteMax bt.size=" + bt.size());
            assertTrue(bt.rootIsSet);
            //System.out.println(" after delete:");
            //bt.printPreOrderTraversal();
            assertFalse(bt.contains(max));
            nodes.removeAt(nodes.size() - 1);
            assertEquals(nodes.size(), bt.size());
            
            long min = nodes.get(0);
            //System.out.println("deleting min=" + min + 
            //    "  bt.rootIsSet=" + bt.rootIsSet + " bt.size=" + bt.size());
            assertTrue(bt.rootIsSet);
            boolean contains = bt.contains(min);
            assertTrue(contains);
            bt.deleteMin();
            //System.out.println("after deleteMin bt.size=" + bt.size());
            assertFalse(bt.contains(min));
            nodes.removeAt(0);
            assertEquals(nodes.size(), bt.size());
            
        }            
    }
    
    public void testKeyOperations0() throws Exception {
    
        System.out.println("testKeyOperations");
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.currentTimeMillis();
        //seed = 1500269991378L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        RedBlackBSTLongInt2 bt = new RedBlackBSTLongInt2();
        
        int n = 100;
        
        int count = 0;
        
        TLongList nodes = new TLongArrayList(2*n);
        
        long[] kOutput = new long[2];
        
        for (int i = 0; i < n/2; ++i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i += skip;
                continue;
            }
            nodes.add(i);
            bt.put(i, i);
            assertTrue(bt.rootIsSet);
            assertTrue(bt.contains(i));
            System.out.println("assert min");
            bt.min(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(nodes.get(0), kOutput[1]);
            
            System.out.println("assert max");
            bt.max(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(i, kOutput[1]);
            
            count++;
        }
        
        for (int i = (n - 1); i >= (n/2); --i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i -= skip;
                continue;
            }
            nodes.add(i);
            bt.put(i, i);
            assertTrue(bt.rootIsSet);
            assertTrue(bt.contains(i));
            System.out.println("assert min");
            bt.min(kOutput);
            assertTrue(kOutput[0] != -1);
            assertEquals(nodes.get(0), kOutput[1]);
            count++;
        }
        nodes.sort();
        int n2 = bt.size();
        assertEquals(count, n2);
        assertEquals(count, nodes.size());
        
        
        int[] vOutput = new int[2];
        
        for (int nIter = 0; nIter < 3; ++nIter){
            n2 = bt.size();
            assertEquals(n2, nodes.size());

            for (int i = 0; i < n2 - 1; ++i) {
                long idx = nodes.get(i);
                long foundIndex = bt.contains(idx) ? idx : -1;
                assertTrue(foundIndex > -1);

                long expected = nodes.get(i + 1);
                System.out.println("\n* " + idx + " expected next=" + expected);

                bt.higher(idx, kOutput);
                assertTrue(kOutput[0] != -1);
                long next = kOutput[1];
                System.out.println(idx + "   next=" + next);
                assertEquals(expected, next);

                assertTrue(bt.rootIsSet);
                
//System.out.println("nIter=" + nIter + " i=" + i + " idx=" + idx);
 
                if (next > (idx + 1)) {
                    // test ceiling of idx+1
                    System.out.println("assert ceiling");
                    bt.ceiling(idx + 1, kOutput);
                    assertTrue(kOutput[0] != -1);
                    long ceil = kOutput[1];
                    assertEquals(expected, ceil);
                }
            }
            
            System.out.println("assert higher");
            bt.higher(nodes.get(nodes.size() - 1), kOutput);
            assertTrue(kOutput[0] == -1);
            System.out.println("assert lower");
            bt.lower(nodes.get(0), kOutput);
            assertTrue(kOutput[0] == -1);
            
            for (int i = 1; i < n2; ++i) {
                long idx = nodes.get(i);
                long foundIndex = bt.contains(idx) ? idx : -1;
                assertTrue(foundIndex > -1);

                System.out.println("assert get");
                bt.get(idx, vOutput);
                assertTrue(vOutput[0] != -1);
                assertEquals(idx, vOutput[1]);
                
                assertTrue(bt.rootIsSet);
                
                long expected = nodes.get(i - 1);
                System.out.println("\n* " + idx + " expected prev=" + expected);

                bt.lower(idx, kOutput);
                assertTrue(kOutput[0] != -1);
                long prev = kOutput[1];
                System.out.println(idx + "   prev=" + prev);
                assertEquals(expected, prev);

                if (prev < (idx - 1)) {
                    // test floor of idx-1
                    System.out.println("assert floor");
                    bt.floor(idx - 1, kOutput);
                    assertTrue(kOutput[0] != -1);
                    long floor = kOutput[1];
                    assertEquals(expected, floor);
                }
            }
           
            if ((nIter & 1) == 1) {
                //randomly remove some nodes
                for (int i = 0; i < n2/4; ++i) {
                    assertTrue(bt.rootIsSet);
                    int idx = rand.nextInt(nodes.size());
                    long v = nodes.get(idx);
                    assertTrue(bt.contains(v));
                    
                    //bt.printPreOrderTraversal();
                    System.out.println("delete " + v + " idx=" + idx);
                    
                    bt.delete(v);
                    assertTrue(bt.rootIsSet);
                    assertFalse(bt.contains(v));
                    nodes.removeAt(idx);
                    assertEquals(nodes.size(), bt.size());
                }
           
            } else {
                //randomly add some nodes
                for (int i = 0; i < n2/4; ++i) {
                    int idx = nodes.size() + rand.nextInt(2*nodes.size());
                    // the nodes contains is linear search, so could use a temp
                    // set if this test gets large one day
                    if (!nodes.contains(idx)) {
                        bt.put(idx, idx);
                        nodes.add(idx);
                        assertTrue(bt.rootIsSet);
                        assertTrue(bt.contains(idx));
                    }
                }
            }
            
            n2 = nodes.size();
            nodes.sort();
           
            assertEquals(nodes.get(nodes.size()/2), bt.select(nodes.size()/2));            
                  
            assertEquals(nodes.size(), bt.size());
            
            System.out.println("before deleteMax bt.size=" + bt.size());
            
            long max = nodes.get(nodes.size() - 1);
            assertTrue(bt.contains(max));
            assertTrue(bt.rootIsSet);
            //System.out.println("will delete max=" + max + " from this tree:");
            //bt.printPreOrderTraversal();
            bt.deleteMax();
            //System.out.println("after deleteMax bt.size=" + bt.size());
            assertTrue(bt.rootIsSet);
            //System.out.println(" after delete:");
            //bt.printPreOrderTraversal();
            assertFalse(bt.contains(max));
            nodes.removeAt(nodes.size() - 1);
            assertEquals(nodes.size(), bt.size());

            long min = nodes.get(0);
            System.out.println("deleting min=" + min + 
                "  bt.rootIsSet=" + bt.rootIsSet + " bt.size=" + bt.size());
            assertTrue(bt.rootIsSet);
            boolean contains = bt.contains(min);
            assertTrue(contains);
            bt.deleteMin();
            //System.out.println("after deleteMin bt.size=" + bt.size());
            assertFalse(bt.contains(min));
            nodes.removeAt(0);
            assertEquals(nodes.size(), bt.size());
        }            
    }
    
}
