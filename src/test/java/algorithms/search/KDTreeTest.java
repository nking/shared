package algorithms.search;

import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.HashSet;
import java.util.Set;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

public class KDTreeTest extends TestCase {

    public void testMemoryEstimates() {
        long est;
        
        long MB = 1024 * 1024;
        
        System.out.println("mem is in MB:");
        
        int[] x = new int[]{
              5000, 1024, 512, 256, 128, 64
        };
        int[] y = new int[]{
              7000, 1024, 512, 256, 128, 64
        };
        
        double rt;
        for (int lenIdx = 0; lenIdx < x.length; ++lenIdx) {
            
            float f = 1.f;
            
            for (int i = 0; i < 2; ++i) {
            
                if (i == 1) {
                    f = 0.1f;
                }
                
                int w = (int)Math.round(Math.log(x[lenIdx] * y[lenIdx])/Math.log(2));
                
                int n = Math.round(f * x[lenIdx] * y[lenIdx]);
                
                est = KDTree.estimateSizeOnHeap(n);
                
                rt = Math.log(n)/Math.log(2);
                
                System.out.format(
                    "width=%4d, height=%4d n=%10d mem=%10d  w=%2d  rt=%4d\n",
                    x[lenIdx], y[lenIdx], n, est/MB, w, (int)Math.round(rt));
                
            }
        }

    }

    public void test() {
	int[] x = new int[] {5, 5, 9, 3, 4, 1, 7, 2};
    	int[] y = new int[] {5, 5, 6, 6, 9, 1, 9, 9};
    	
    	int lastUsableIndex = KDTree.reduceToUniqueWithMoveUp(x, y);
    	
	assertTrue(lastUsableIndex == 6);
		
	int[] expectedx = new int[] {5, 9, 3, 4, 1, 7, 2};
    	int[] expectedy = new int[] {5, 6, 6, 9, 1, 9, 9};
    	
    	for (int i = 0; i < expectedx.length; i++) {
            assertTrue(expectedx[i] == x[i]);
            assertTrue(expectedy[i] == y[i]);
    	}
    }

    public void testTree1() {        
    	
    	int[] x = new int[] {6, 5, 9, 3, 4, 1, 7, 2};
    	int[] y = new int[] {1, 5, 6, 6, 9, 1, 9, 9};
    	
    
    	KDTree kdtree = new KDTree(x, y);
    	kdtree.printTree();
    	/*
    	 * buildTree(0, 0, 7)[left 1]: buildTree(1, 0, 3)[left 2]: buildTree(2, 0, 1)[left 3]: buildTree(3, 0, 0)  1 1

[right 3]: buildTree(3, 1, 1)  3 6


[right 2]: buildTree(2, 2, 3)[left 3]: buildTree(3, 2, 2)  2 9

[right 3]: buildTree(3, 3, 3)  4 9



[right 1]: buildTree(1, 4, 7)[left 2]: buildTree(2, 4, 5)[left 3]: buildTree(3, 4, 4)  5 5

[right 3]: buildTree(3, 5, 5)  6 1


[right 2]: buildTree(2, 6, 7)[left 3]: buildTree(3, 6, 6)  7 9

[right 3]: buildTree(3, 7, 7)  9 6
    	 */
    	/*
    	 * {6, 5, 9, 3, 4, 1, 7, 2};
    	 * {1, 5, 6, 6, 9, 1, 9, 9};
    	 * 
    	 * depth = 0 so sort by x
    	 * {1, 2, 3, 4, 5, 6, 7, 9};
    	 * {1, 9, 6, 9, 5, 1, 9, 6}
    	 * 
    	 * left depth=1 sort by y
    	 *    {1, 2, 3, 4}
    	 *    {1, 9, 6, 9}
    	 *    
    	 *    {1, 3, 2, 4}
    	 *    {1, 6, 9, 9}
    	 *    
    	 *    left, depth=2, sort by x
    	 *    {1, 3}
    	 *    {1, 6}
    	 *    return left.left .left and .right from the immed lines above <===
    	 *    
    	 *    right, depth=2, sort by x
    	 *    {2, 4}
    	 *    {9, 9}
    	 *    return left.right .left and .right from the immed lines above <===
    	 *    
    	 * right depth=1 sort by y
    	 *    {5, 6, 7, 9}
    	 *    {5, 1, 9, 6}
    	 * 
    	 *    {6, 5, 9, 7}
    	 *    {1, 5, 6, 9}
    	 *    
    	 *    left, depth = 2 sort by x
    	 *    {5, 6}
    	 *    {5, 1}
    	 *    return right.left .left and .right from immed lines above <===
    	 *    
    	 *    right, depth = 2 sort by x 
    	 *    {7, 9}
    	 *    {9, 6}
    	 *    return right.right .left and .right from immed lines above <===
    	 * 
    	 :LEFT 4:LEFT 6:LEFT 10(1,1)
		 :LEFT 4:LEFT 6:RIGHT 10(3,6)
		 :LEFT 4:RIGHT 6:LEFT 20(2,9)
		 :LEFT 4:RIGHT 6:RIGHT 20(4,9)
		 :RIGHT 4:LEFT 5:LEFT 50(5,5)
		 :RIGHT 4:LEFT 5:RIGHT 50(6,1)
		 :RIGHT 4:RIGHT 5:LEFT 70(7,9)
		 :RIGHT 4:RIGHT 5:RIGHT 70(9,6)
    	 */
    	assertTrue(kdtree.root.left.left.left.x == 1);
    	assertTrue(kdtree.root.left.left.left.y == 1);
    	assertTrue(kdtree.root.left.left.right.x == 3);
    	assertTrue(kdtree.root.left.left.right.y == 6);
    	
    	assertTrue(kdtree.root.left.right.left.x == 2);
    	assertTrue(kdtree.root.left.right.left.y == 9);
    	assertTrue(kdtree.root.left.right.right.x == 4);
    	assertTrue(kdtree.root.left.right.right.y == 9);
    	
    	assertTrue(kdtree.root.right.left.left.x == 5);
    	assertTrue(kdtree.root.right.left.left.y == 5);
    	assertTrue(kdtree.root.right.left.right.x == 6);
    	assertTrue(kdtree.root.right.left.right.y == 1);
    	
    	assertTrue(kdtree.root.right.right.left.x == 7);
    	assertTrue(kdtree.root.right.right.left.y == 9);
    	assertTrue(kdtree.root.right.right.right.x == 9);
    	assertTrue(kdtree.root.right.right.right.y == 6);
    	
    	assertTrue(kdtree.root.key == 4);
    	assertTrue(kdtree.root.left.key == 6);
    	assertTrue(kdtree.root.left.left.key == 1);
    	
    }

    public void testTree2() {        
    	
    	int[] x = new int[] {6, 5, 9, 3, 4, 4, 7, 2};
    	int[] y = new int[] {1, 5, 6, 6, 9, 0, 9, 9};
    	
    	/* sorted by x because depth=0 then divided by 2
    	 * 0  1  2  3    4  5  6  7
    	  {2, 3, 4, 4,   5, 6, 7, 9};
    	  {9, 6, 9, 0,   5, 1, 9, 6}
    	  
    	  
    	  left depth=1 sort by y
    	  {4, 3, 2, 4}
    	  {0, 6, 9, 9}
    	        
    	      left depth=2 sort by x
    	      {3, 4}
    	      {6, 0}
    	      
    	      right depth=2 sort by x
    	      {2, 4}
    	      {9, 9}
    	      
    	  
    	  right depth=1 sort by y
    	  {6, 5, 9, 7}
    	  {1, 5, 6, 9}
    	  
    	      left depth=2 sort by x
    	      {5, 6}
    	      {5, 1}
    	      
    	      right depth=2 sort by x
    	      {7, 9}
    	      {9, 6}

		 :LEFT 4:LEFT 6:LEFT 30(3,6)
		 :LEFT 4:LEFT 6:RIGHT 30(4,0)
		 :LEFT 4:RIGHT 6:LEFT 20(2,9)
		 :LEFT 4:RIGHT 6:RIGHT 20(4,9)
		 :RIGHT 4:LEFT 5:LEFT 50(5,5)
		 :RIGHT 4:LEFT 5:RIGHT 50(6,1)
		 :RIGHT 4:RIGHT 5:LEFT 70(7,9)
		 :RIGHT 4:RIGHT 5:RIGHT 70(9,6)
    	*/
    	KDTree kdtree = new KDTree(x, y);
    	kdtree.printTree();
    	
    	assertTrue(kdtree.root.left.left.left.x == 3);
    	assertTrue(kdtree.root.left.left.left.y == 6);
    	assertTrue(kdtree.root.left.left.right.x == 4);
    	assertTrue(kdtree.root.left.left.right.y == 0);
    	
    	assertTrue(kdtree.root.left.right.left.x == 2);
    	assertTrue(kdtree.root.left.right.left.y == 9);
    	assertTrue(kdtree.root.left.right.right.x == 4);
    	assertTrue(kdtree.root.left.right.right.y == 9);
    	
    	assertTrue(kdtree.root.right.left.left.x == 5);
    	assertTrue(kdtree.root.right.left.left.y == 5);
    	assertTrue(kdtree.root.right.left.right.x == 6);
    	assertTrue(kdtree.root.right.left.right.y == 1);
    	
    	assertTrue(kdtree.root.right.right.left.x == 7);
    	assertTrue(kdtree.root.right.right.left.y == 9);
    	assertTrue(kdtree.root.right.right.right.x == 9);
    	assertTrue(kdtree.root.right.right.right.y == 6);
    }

    public void testSearch() {
       
    	int[] x = new int[] {6, 5, 9, 3, 4, 4, 7, 2};
    	int[] y = new int[] {1, 5, 6, 6, 9, 0, 9, 9};
    	
    	KDTree kdtree = new KDTree(x, y);
    	KDTreeNode node = kdtree.findNearestNeighbor(6, 8);
    	assertNotNull(node);
    	assertEquals(7, node.x);
    	assertEquals(9, node.y);
    }
    	
    public void test3() {
        
	int[] x = new int[] {5, 5, 9, 13, 14, 21, 27, 32};
    	int[] y = new int[] {5, 5, 6,  6,  9,  1,  9,  9};
    	
        KDTree kdTree = new KDTree(x, y);
    	
        kdTree.printTree();
        
        KDTreeNode node = kdTree.findNearestNeighbor(12, 7);
    	assertNotNull(node);
    	assertEquals(13, node.x);
    	assertEquals(6, node.y);
    }
    
    public void test3_2() {
        
        int width = 35;
        int height = 10;
	int[] x = new int[] {5, 5, 9, 13, 14, 21, 27, 32};
    	int[] y = new int[] {5, 5, 6,  6,  9,  1,  9,  9};
    	
        PixelHelper ph = new PixelHelper();
        TLongSet pixIdxs = new TLongHashSet();
        for (int i = 0; i < x.length; ++i) {
            long pixIdx = ph.toPixelIndex(x[i], y[i], width);
            pixIdxs.add(pixIdx);
        }
        
        KDTree kdTree = new KDTree(pixIdxs, width, height);
    	
        kdTree.printTree();
        
        KDTreeNode node = kdTree.findNearestNeighbor(12, 7);
    	assertNotNull(node);
    	assertEquals(13, node.x);
    	assertEquals(6, node.y);
        assertTrue(node.xyAreSet());

        node = new KDTreeNode();
        assertFalse(node.xyAreSet());
    }
    
    public void testNN_ne() {
        
	int[] x = new int[] {5, 5, 9, 13, 14, 21, 27, 32};
    	int[] y = new int[] {5, 5, 6,  6,  9,  1,  9,  9};
    	
        KDTree kdTree = new KDTree(x, y);
    	
        kdTree.printTree();
        
        KDTreeNode node = kdTree.findNearestNeighbor(12, 7);
    	assertNotNull(node);
    	assertEquals(13, node.x);
    	assertEquals(6, node.y);
        
        node = kdTree.findNearestNeighborNotEquals(13, 6);
        // 9, 6 or 14, 9
    	assertNotNull(node);
        PairInt found = new PairInt(node.x, node.y);
        
    	Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(9, 6));
        expected.add(new PairInt(14, 9));
        
        assertTrue(expected.remove(found));
    	assertEquals(1, expected.size());
        
    }
    
    /**
     * Test suite
     @return static Test
    */
    public static Test suite(){
        System.out.println("Creating a TestSuite for"
            + " KDTree");

        return new TestSuite(KDTreeTest.class);
    }

    /**
     * Set up a Junit test runner
     @param args Not used.
    */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }

}
