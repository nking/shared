package algorithms.search;

import algorithms.sort.MiscSorter;
import algorithms.util.ObjectSpaceEstimator;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import gnu.trove.iterator.TLongIterator;
import gnu.trove.set.TLongSet;
import java.util.HashSet;
import java.util.Set;
import java.util.Arrays;

/**
 * k-dimension tree is a binary tree used to store coordinates for quick nearest
 * neighbor or range searches.
 *    -- values are always stored in leaves
 *    -- the meaning of an internal node depends upon the depth of the
 *       node within the tree.
 *       
 * Note that learning the true medium while 
 * constructing the tree takes more
 * time, but leads to a better balanced tree.
 * 
 * adapted from pseudocode from
 * http://ldots.org/kdtree which licenses the content as:
 * http://creativecommons.org/licenses/by-sa/2.0/
 * 
 * useful reading regarding best distances in nearest neighbor search:
 *    http://web.stanford.edu/class/cs106l/handouts/assignment-3-kdtree.pdf
 * 
 * construction is O(N*log_2(N))
 * queries are: O(log_2(N)) to O(sqrt(N) + k)?
 *              need to revisit code to estimate this...
 * 
 * first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
<pre>
   mem is in MB: 
     width=5000, height=7000 n=  35000000 mem=      1869  w=25  rt=  25
     width=5000, height=7000 n=   3500000 mem=       186  w=25  rt=  22
     width=1024, height=1024 n=   1048576 mem=        56  w=20  rt=  20
     width=1024, height=1024 n=    104858 mem=         5  w=20  rt=  17
     width= 512, height= 512 n=    262144 mem=        14  w=18  rt=  18
     width= 512, height= 512 n=     26214 mem=         1  w=18  rt=  15
     width= 256, height= 256 n=     65536 mem=         3  w=16  rt=  16
     width= 256, height= 256 n=      6554 mem=         0  w=16  rt=  13
     width= 128, height= 128 n=     16384 mem=         0  w=14  rt=  14
     width= 128, height= 128 n=      1638 mem=         0  w=14  rt=  11
     width=  64, height=  64 n=      4096 mem=         0  w=12  rt=  12
     width=  64, height=  64 n=       410 mem=         0  w=12  rt=   9
  </pre>
 * @author nichole
 */
public class KDTree {
	
    /**
     *
     */
    protected KDTreeNode root = null;
		
    /**
     *
     @param x
     @param y
     */
    public KDTree(int[] x, int[] y) {
		
		if (x == null) {
			throw new IllegalArgumentException("x cannot be null");
		}
		if (y == null) {
			throw new IllegalArgumentException("y cannot be null");
		}
		if (x.length < 2) {
			throw new IllegalArgumentException("x must be larger than 2");
		}
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y must have same number of points");
		}
	    
		int lastUsableIndex = reduceToUniqueWithMoveUp(x, y);

        if (lastUsableIndex < (x.length - 1)) {
            x = Arrays.copyOf(x, lastUsableIndex + 1);
            y = Arrays.copyOf(y, lastUsableIndex + 1);
        }

        this.root = buildTree(0, x, y, 0, lastUsableIndex);
    }

    /**
     *
     @param pixIdxs
     @param width
     @param height
     */
    public KDTree(TLongSet pixIdxs, int width, int height) {
		
		if (pixIdxs == null || pixIdxs.size() < 2) {
			throw new IllegalArgumentException(
                "pixIdxs cannot be null or less than 2 in size");
		}
        
        PixelHelper ph = new PixelHelper();
        int n = pixIdxs.size();
        int[] x = new int[n];
        int[] y = new int[n];
        int[] xy = new int[2];
        TLongIterator iter = pixIdxs.iterator();
        int count = 0;
        while (iter.hasNext()) {
            ph.toPixelCoords(iter.next(), width, xy);
            x[count] = xy[0];
            y[count] = xy[1];
            count++;
        }
        pixIdxs = null;
	    
		int lastUsableIndex = reduceToUniqueWithMoveUp(x, y);

        if (lastUsableIndex < (x.length - 1)) {
            x = Arrays.copyOf(x, lastUsableIndex + 1);
            y = Arrays.copyOf(y, lastUsableIndex + 1);
        }

        this.root = buildTree(0, x, y, 0, lastUsableIndex);
    }
    
    /**
     *
     @return
     */
    public KDTreeNode getRoot() {
		return root;
	}
	
	/**
	 * remove non-unique values by moving up items underneath them.  returns
	 * the last index which should be used in the modified arrays given as arguments.
	 @param x x coordinates
	 @param y y coordinates
	 @return the number of items to use in the modified x and y arrays which have only unique
	 * pairs in them now.
	 */
	static int reduceToUniqueWithMoveUp(int[] x, int[] y) {
		// reduce points to unique by moving them up if they already exist.
            
        Set<PairInt> added = new HashSet<PairInt>();
        
		int count = 0;
		boolean moveUp = false;
		for (int i = 0; i < x.length; i++) {
			if (moveUp) {
				x[count] = x[i];
				y[count] = y[i];
			}
            PairInt p = new PairInt(x[i], y[i]);
			if (!added.contains(p)) {
				added.add(p);
				count++;
			} else {
				moveUp = true;
			}
		}
		return (count-1);
	}
	
	/**
	 * sort by x and return index of median x value for which all points 
	 * below are less than median x and .geq. for all points above median x index,
	 * in other words 
	 * for indexes : 0 .lt. index have values x .lt. x[index] where x[index] is median x.
	 * for indexes : index .lte. n have values x .gte. x[index] where x[index] is median x.
	 @param x x coordinates
	 @param y y coordinates
     @param startSortRange
     @param stopSortRangeExclusive
	 @return index that divides the arrays as x .lt. median value and x .geq. median value
	 */
	int partitionByX(int[] x, int[] y, int startSortRange, int stopSortRangeExclusive) {
				 
        MiscSorter.sortBy1stArg2(x, y, startSortRange, stopSortRangeExclusive);
        
		int n = (stopSortRangeExclusive - startSortRange);
		int index = startSortRange + (n >> 1);
		int xMedian = x[index];
		// find last index with same value as xMedian in range
		int index0 = MiscBisectingSearch.ceiling(x, xMedian, index, stopSortRangeExclusive-1);
		return index;
	}
	
	/**
	 * sort by y and return index of median y value for which all points 
	 * below are less than median y and .geq. for all points above median y index,
	 * in other words 
	 * for indexes : 0 .lt. index have values y .lt. y[index] where y[index] is median y.
	 * for indexes : index .lte. n have values y .gte. y[index] where y[index] is median y.
	 @param x x coordinates
	 @param y y coordinates
     @param startSortRange
     @param stopSortRangeExclusive
	 @return index that divides the arrays as y .lt. median value and y .geq. median value
	 */
	int partitionByY(int[] x, int[] y, int startSortRange, int stopSortRangeExclusive) {
				
        MiscSorter.sortBy1stArg2(y, x, startSortRange, stopSortRangeExclusive);
		int n = (stopSortRangeExclusive - startSortRange);
		int index = startSortRange + (n >> 1); // rounds towards zero
		int yMedian = y[index];
		// find last index with same value as xMedian in range
		int index0 = MiscBisectingSearch.ceiling(y, yMedian, index, stopSortRangeExclusive-1);
		return index;
	}
        		
	/**
	 * build a tree at from given depth for x, y points using a depth-first algorithm.
	 * Exits from the recursion when npoints = 0, 1 or a leaf.
	 * 
	 @param depth depth from which to build tree
	 @param x x coordinates
	 @param y y coordinates
     @param startSortRange index start range to use in partition
     @param stopSortRangeExclusive index stop range to use in partition
	 * 
	 @return tree
	 */
	protected KDTreeNode buildTree(int depth, int[] x, int[] y, 
        int startSortRange, int stopSortRangeExclusive) {

		if (x == null || y == null || x.length == 0 || y.length == 0) {
			return null;
		}
				
		if (stopSortRangeExclusive == startSortRange) {
			// return a leaf of 2 points only
			KDTreeNode leaf = new KDTreeNode();
			leaf.x = x[startSortRange];
			leaf.y = y[startSortRange];			
			return leaf;
		}
		
		int medianIndex = -1;
		int median = 1;
		
	    // if depth of tree is even, partition the points by x, else y
		if (depth % 2 == 0) {
			// sort by x and return last index within range that has same value as median
			medianIndex = partitionByX(x, y, startSortRange, stopSortRangeExclusive);
			median = x[medianIndex];
		} else {
			//sort by y and return last index within range that has same value as median
			medianIndex = partitionByY(x, y, startSortRange, stopSortRangeExclusive);
			median = y[medianIndex];
		}
		
		depth++;
		
		// left points are  startSortRange through medianX
		KDTreeNode leftChildren = buildTree(depth, x, y, startSortRange, medianIndex);

		// right points are medianIndex    through stopSortRangeExclusive
		KDTreeNode rightChildren = buildTree(depth, x, y, medianIndex+1, stopSortRangeExclusive);
		
		KDTreeNode parent = new KDTreeNode();
		
		parent.key = median;
	    
	    if (leftChildren != null) {
			parent.left = leftChildren;
			leftChildren.parent = parent;
	    	parent.nChildren += 1 + parent.left.nChildren;
	    }
	    if (rightChildren != null) {
			parent.right = rightChildren;
			rightChildren.parent = parent;
	    	parent.nChildren += 1 + parent.right.nChildren;
	    }
	    
	    return parent;
	}
        
    /**
     *
     @param x
     @param y
     @return
     */
    public KDTreeNode findNearestNeighbor(int x, int y) {

        KDTreeNode[] best = new KDTreeNode[1];
        double[] bestDist = new double[]{Double.MAX_VALUE};

        nearestNeighborSearch(root, x, y, 0, best, bestDist, false);
        
        return best[0];
    }
    
    /**
     *
     @param x
     @param y
     @return
     */
    public KDTreeNode findNearestNeighborNotEquals(int x, int y) {

        KDTreeNode[] best = new KDTreeNode[1];
        double[] bestDist = new double[]{Double.MAX_VALUE};

        nearestNeighborSearch(root, x, y, 0, best, bestDist, true);
        
        return best[0];
    }
	
    /**
     *
     @param tree
     @param leftValue
     @param rightValue
     @param depth
     @param best
     @param bestDist
     @param excludeEquals
     @return
     */
    protected KDTreeNode nearestNeighborSearch(KDTreeNode tree, int leftValue, 
        int rightValue, int depth, KDTreeNode[] best, double[] bestDist,
        boolean excludeEquals) {

        if (tree.nChildren == 0 ) {
			return tree;
		}
		
		int medianValue = tree.getKey();
		
		float diffMedValSq;
		
		KDTreeNode subTree1, subTree2;
        int partition;

        if ((depth & 1) == 0) {
            partition = leftValue;
        } else {
            partition = rightValue;
        }
        
        diffMedValSq = medianValue - partition;
        if (partition <= medianValue) {
            subTree1 = tree.left;
            subTree2 = tree.right;
        } else {
            subTree1 = tree.right;
            subTree2 = tree.left;
        }
        diffMedValSq *= diffMedValSq;
	
        KDTreeNode retVal1 = nearestNeighborSearch(
            subTree1, leftValue, rightValue, depth + 1, best, bestDist, 
            excludeEquals);
		
        double dist1 = Double.MAX_VALUE;
        if (retVal1 != null) {
            dist1 = distanceSq(retVal1, leftValue, rightValue);
            if (excludeEquals) {
                if (dist1 > 0) {
                    if (dist1 < bestDist[0]) {
                        bestDist[0] = dist1;
                        best[0] = retVal1;
                    }
                }
            } else {
                if (dist1 < bestDist[0]) {
                    bestDist[0] = dist1;
                    best[0] = retVal1;
                    if (dist1 == 0) {
                        // this is the point
                        return retVal1;
                    }
                }
            }
        }
        
        //TODO: this may need to be revised for a radius.
        //   basically, if (leftValue, rightValue) is closer to
        //      the median than it is to retVal1,
        //      search the other tree too.
        
		if ((2*diffMedValSq) < dist1) {
            
			KDTreeNode retVal2 = nearestNeighborSearch(
                subTree2, leftValue, rightValue, depth + 1, best, bestDist,
                excludeEquals);
            
            double dist2 = Double.MAX_VALUE;
            
            if (retVal2 != null) {
                dist2 = distanceSq(retVal2, leftValue, rightValue);
                if (excludeEquals) {
                    if (dist2 > 0) {
                        if (dist2 < bestDist[0]) {
                            bestDist[0] = dist2;
                            best[0] = retVal2;
                        }
                    }
                } else {
                    if (dist2 < bestDist[0]) {
                        bestDist[0] = dist2;
                        best[0] = retVal2;
                        if (dist2 == 0) {
                            // this is the point
                            return retVal2;
                        }
                    }
                }
            }
        }
        
		return best[0];
    }
    
    private double distanceSq(KDTreeNode tree, 
        float leftValue, float rightValue) {

        float diffX = tree.getX() - leftValue;
        float diffY = tree.getY() - rightValue;
        
        return (diffX * diffX) + (diffY * diffY);
    }
    
    /**
     *
     */
    public void printTree() {
		printTree(root, " ");
	}
	private void printTree(KDTreeNode node, String preString) {
		if (node.left != null) {
		    printTree(node.left, preString + ":LEFT " + node.getKey());
		}
		if (node.right != null) {
		    printTree(node.right, preString + ":RIGHT " + node.getKey());
		}
		if (node.left == null && node.right == null) {
			System.out.println(preString + node.getKey() + "(" + node.getX() + "," + node.getY() + ")");
		}
	}
    
    /**
     *
     @param numberOfPoints
     @return
     */
    public static long estimateSizeOnHeap(int numberOfPoints) {
                
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNObjRefsFields(1);
        
        long total = est.estimateSizeOnHeap();
        
        total += numberOfPoints * KDTreeNode.estimateSizeOnHeap();
        
        return total;
    }
}
