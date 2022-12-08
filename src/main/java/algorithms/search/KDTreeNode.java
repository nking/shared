package algorithms.search;

import algorithms.util.ObjectSpaceEstimator;

/**
 *
 * @author nichole
 */
public class KDTreeNode {

	KDTreeNode right = null;
	KDTreeNode left = null;
	KDTreeNode parent = null;
	
    /**
     *
     */
    public static final int sentinel = Integer.MIN_VALUE;
	int x = sentinel;
	int y = sentinel;
	int key = sentinel; // median value
	int nChildren = 0;

    /**
     *
     @return
     */
    public boolean xyAreSet() {
        return (x != sentinel && y != sentinel);
    }
    
    /**
     @return the x
     */
    public int getX() {
        return x;
    }
    
    /**
     @param theX
     */
    public void setX(int theX) {
        x = theX;
    }

    /**
     @return the y
     */
    public int getY() {
        return y;
    }
    
    /**
     @param theY
     */
    public void setY(int theY) {
        y = theY;
    }

    /**
     @return the key
     */
    public int getKey() {
        return key;
    }
    
    /**
     *
     @return
     */
    public static long estimateSizeOnHeap() {
                
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNIntFields(4);
        est.setNObjRefsFields(3);
        
        return est.estimateSizeOnHeap();
    }
    
}
