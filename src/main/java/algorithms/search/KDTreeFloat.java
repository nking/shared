package algorithms.search;

import algorithms.sort.MiscSorter;
import algorithms.util.PairFloat;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashSet;
import java.util.Set;

/**
 * k-dimension tree is a binary tree used to store coordinates for quick nearest
 * neighbor or range searches.
 * This one is a 2D tree only.
 * 
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
 * first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

 * @author nichole
 */
public class KDTreeFloat {
	
    /**
     *
     */
    protected KDTreeNodeFloat root = null;
		
    /**
     *
     @param xPoints
     @param yPoints
     @param alreadySorted
     */
    public KDTreeFloat(float[] xPoints, float[] yPoints, boolean alreadySorted) {
		
		if (xPoints == null) {
			throw new IllegalArgumentException("xPoints cannot be null");
		}
		if (yPoints == null) {
			throw new IllegalArgumentException("yPoints cannot be null");
		}
		if (xPoints.length < 2) {
			throw new IllegalArgumentException("xPoints must be larger than 2");
		}
		if (xPoints.length != yPoints.length) {
			throw new IllegalArgumentException("xPoints and yPoints must have same number of points");
		}
        
        int lastUsableIndex = reduceToUniqueWithMoveUp(xPoints, yPoints);
	    
        if (lastUsableIndex < (xPoints.length - 1)) {
            xPoints = Arrays.copyOf(xPoints, lastUsableIndex + 1);
            yPoints = Arrays.copyOf(yPoints, lastUsableIndex + 1);
        }
        
        //TODO: to better handle space for large number of points,
        //  could change out the merge sorts for quick sorts
        
        if (!alreadySorted) {
            MiscSorter.sortBy1stArgThen2nd(xPoints, yPoints);
        }
                
        this.root = buildTree(0, xPoints, yPoints, 0, lastUsableIndex);
	}
    
	/**
	 * remove unique values by moving up items underneath them.  returns
	 * the last index which should be used in the modified arrays given as arguments.
	 @param x x coordinates
	 @param y y coordinates
	 @return the number of indexes usable in x and y after the arrays have been condensed to remove redundant pairs
	 */
	static int reduceToUniqueWithMoveUp(float[] x, float[] y) {
		// reduce points to unique by moving them up if they already exist.
            
        Set<PairFloat> added = new HashSet<PairFloat>();
        
		int count = 0;
		// [0]  0,0      count=0  0,0
		// [1]  1,1            1  1,1
		// [2]  2,2            2  2,2
		// [3]  2,2            *
		// [4]  3,3            3  3,3 
		// [5]  4,4            4  4,4
		// [6]  4,4            *
		// [7]  5,5            5  5,5
		boolean moveUp = false;
		for (int i = 0; i < x.length; i++) {
			if (moveUp) {
				x[count] = x[i];
				y[count] = y[i];
			}
            PairFloat p = new PairFloat(x[i], y[i]);
			if (!added.contains(p)) {
				added.add(p);
				count++;
			} else {
				moveUp = true;
			}
		}
        
		return added.size() - 1;
	}
	
    /**
	 * build a tree at from given depth for x, y points using a depth-first algorithm.
	 * Exits from the recursion when npoints = 0, 1 or a leaf.
	 * 
	 @param depth depth from which to build tree
	 @param x x coordinates
	 @param y y coodinates
     @param startSortRange
     @param stopSortRangeExclusive
	 * 
	 @return tree
	 */
	protected KDTreeNodeFloat buildTree(int depth, float[] x, float[] y, 
        int startSortRange, int stopSortRangeExclusive) {

		if (x == null || y == null || x.length == 0 || y.length == 0) {
			return null;
		}
				
		if (stopSortRangeExclusive == startSortRange) {
			KDTreeNodeFloat leaf = new KDTreeNodeFloat();
            leaf.depth = depth;
            leaf.key = Integer.MIN_VALUE;
			leaf.x = x[startSortRange];
			leaf.y = y[startSortRange];			
			return leaf;
		}
		
		int medianIndex = -1;
		float median = 1;
		
	    // if depth of tree is even, partition the points by x, else y
		if ((depth & 1) == 0) {
			medianIndex = partitionByX(
                x, y, startSortRange, stopSortRangeExclusive);
			median = x[medianIndex];
		} else {
			medianIndex = partitionByY(
                x, y, startSortRange, stopSortRangeExclusive);
			median = y[medianIndex];
		}
                
		depth++;
       
		// left points are  startSortRange through medianX
		KDTreeNodeFloat leftChildren = buildTree(depth, x, y, 
            startSortRange, medianIndex);
		
		// right points are medianIndex    through stopSortRangeExclusive
		KDTreeNodeFloat rightChildren = buildTree(depth, x, y, 
            medianIndex+1, stopSortRangeExclusive);
		
		KDTreeNodeFloat parent = new KDTreeNodeFloat();
		
        parent.depth = depth - 1;
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
	 * sort by x and return index of median x value for which all points 
	 * below are less than median x and .geq. for all points above median x index,
	 * in other words 
	 * for indexes : 0 .lt. index have values x .lt. x[index] where x[index] is median x.
	 * for indexes : index .leq. n have values x .geq. x[index] where x[index] is median x.
	 @param x x coordinates
	 @param y y coordinates
     @param startSortRange
     @param stopSortRangeExclusive
	 @return index that divides the arrays as x .lt. median value and x .geq. median value
	 */
	int partitionByX(float[] x, float[] y, int startSortRange, int stopSortRangeExclusive) {
				
        MiscSorter.sortBy1stArg(
            x, y, startSortRange, stopSortRangeExclusive);
        
		int n = (stopSortRangeExclusive - startSortRange);
		int index = startSortRange + (n >> 1);
		float xMedian = x[index];
		while ((index+1) < stopSortRangeExclusive) {
			if (x[index + 1] == xMedian) {
				index++;
			} else {
				break;
			}
		}
		return index;
	}
	
	/**
	 * sort by y and return index of median y value for which all points 
	 * below are less than median y and .geq. for all points above median y index,
	 * in other words 
	 * for indexes : 0 .lt. index have values y .lt. y[index] where y[index] is median y.
	 * for indexes : index .leq. n have values y .geq. y[index] where y[index] is median y.
	 @param x x coordinates
	 @param y y coordinates
     @param startSortRange
     @param stopSortRangeExclusive
	 @return index that divides the arrays as y .lt. median value and y .geq. median value
	 */
	int partitionByY(float[] x, float[] y, int startSortRange, int stopSortRangeExclusive) {
				
        MiscSorter.sortBy1stArg(y, x, startSortRange, stopSortRangeExclusive);
		int n = (stopSortRangeExclusive - startSortRange);
		int index = startSortRange + (n >> 1); // rounds towards zero
		float yMedian = y[index];
		while ((index+1) < stopSortRangeExclusive) {
			if (y[index + 1] == yMedian) {
				index++;
			} else {
				break;
			}
		}
		return index;
	}	
    
    private Set<KDTreeNodeFloat> bestNode = null;
    private double bestDist = Double.MAX_VALUE;

    /**
     *
     */
    protected Set<KDTreeNodeFloat> visited = null;
    
    /**
     * find the nearest neighbor, and if it is equidistant to
     * others, return those too.
     * Note that in the worse case, a brute force search over
     * all members would be faster than this which can follow
     * each branch if the equidistant points are left and
     * right of root and close to equidistant to the median split for each level.
     @param x x coordinate
     @param y y coordinate
     @return nearest neighbors
     */
	public Set<PairFloat> findNearestNeighbor(float x, float y) {
        
        bestNode = new HashSet<KDTreeNodeFloat>();
        bestDist = Double.MAX_VALUE;
        visited = new HashSet<KDTreeNodeFloat>();
        
        Set<KDTreeNodeFloat> nodes = 
            nearestNeighborSearch(root, x, y, 0);
        
        if (nodes == null || nodes.size() == 0) {
            return null;
        }
        
        Set<PairFloat> set = new HashSet<PairFloat>(nodes.size());
        for (KDTreeNodeFloat node : nodes) {
            PairFloat p = new PairFloat(node.getX(), node.getY());    
            set.add(p);
        }
        return set;
	}

    /**
     *
     @param tree
     @param leftValue
     @param rightValue
     @param depth
     @return
     */
    protected Set<KDTreeNodeFloat> nearestNeighborSearch(
        KDTreeNodeFloat tree, float leftValue, float rightValue,
        int depth) {
 
        //TODO: make this iterative for java...
        
		if (tree.nChildren == 0 ) {
            double dist = distanceSq(tree, leftValue, rightValue);
			if (dist == bestDist) {
                bestNode.add(tree);
            } else if (dist < bestDist) {
                bestNode.clear();
                bestDist = dist;
                bestNode.add(tree);
            }
            return bestNode;
		}
        
		float medianValue = tree.getKey();

        float diffMedValSq;
		
		KDTreeNodeFloat subTree1, subTree2;

        if ((depth & 1) == 0) {
            diffMedValSq = medianValue - leftValue;
            if (leftValue <= medianValue) {
                subTree1 = tree.left;
                subTree2 = tree.right;
            } else {
                subTree1 = tree.right;
                subTree2 = tree.left;
            }
        } else {
            diffMedValSq = medianValue - rightValue;
            if (rightValue <= medianValue) {
                subTree1 = tree.left;
                subTree2 = tree.right;
            } else {
                subTree1 = tree.right;
                subTree2 = tree.left;
            }
        }
        diffMedValSq *= diffMedValSq;
	 
        Set<KDTreeNodeFloat> retVal1 = null;
        if (!visited.contains(subTree1)) {
		    retVal1 = nearestNeighborSearch(
                subTree1, leftValue, rightValue, depth + 1);
		    visited.add(subTree1);
        }
        
        double dist1 = Double.MAX_VALUE;
        if (retVal1 != null && !retVal1.isEmpty()) {
            dist1 = distanceSq(retVal1.iterator().next(), 
                leftValue, rightValue);
            // TODO: consider a tolerance
            if (dist1 == 0) {
                // this is the point
                bestDist = dist1;
                bestNode = retVal1;
                return retVal1;
            }
        }
        
        //TODO: this may need to be revised for a radius.
        //   basically, if (leftValue, rightValue) is closer to
        //      the median than it is to retVal1,
        //      search subtree2 too.
        
        //System.out.println("dist1=" + dist1 + " (med-val)=" + diffMedValSq
        //+ " best=" + bestDist);        
        
		if (!visited.contains(subTree2) && diffMedValSq < dist1) {
			
            Set<KDTreeNodeFloat> retVal2 = nearestNeighborSearch(
                subTree2, leftValue, rightValue, depth + 1);
            
            visited.add(subTree2);
            
            double dist2 = Double.MAX_VALUE;
            if (retVal2 != null && !retVal2.isEmpty()) {
                dist2 = distanceSq(retVal2.iterator().next(), 
                    leftValue, rightValue);
                // TODO: consider a tolerance
                if (dist2 == 0) {
                    // this is the point
                    bestDist = dist2;
                    bestNode = retVal2;
                    return bestNode;
                }
                if (dist2 == dist1) {
                    if (dist1 == bestDist) {
                        bestNode.addAll(retVal2);
                    } else if (dist1 < bestDist) {
                        bestNode.clear();
                        bestDist = dist2;
                        bestNode.addAll(retVal2);
                    }
                } else if (dist2 < dist1) {
                    dist1 = dist2;
                    retVal1 = retVal2;
                }
            }
        }
        
        if (dist1 == bestDist && retVal1 != null) {
            bestNode.addAll(retVal1);
        } else if (dist1 < bestDist && retVal1 != null) {
            bestNode.clear();
            bestDist = dist1;
            bestNode.addAll(retVal1);
        }
        
		return bestNode;
	}
        		
    /**
     *
     */
    public void printTree() {
		printTree(root, " ");
	}
	private void printTree(KDTreeNodeFloat node, String preString) {
		if (node == null) {
            return;
        }
        
        Deque<KDTreeNodeFloat> q0 = new ArrayDeque<KDTreeNodeFloat>();
        Deque<KDTreeNodeFloat> q1 = new ArrayDeque<KDTreeNodeFloat>();
        q0.offer(node);

        int count = 0;
        boolean skip = true;
        while(!q0.isEmpty()) {
            while(!q0.isEmpty()) {
                node = q0.poll();
                System.out.println("level=" + node.depth + " (med=" + node.getKey() +
                    " x=" + node.getX() + " y=" + node.getY() + ")");
                if (node.left != null) {
                    q1.offer(node.left);
                }
                if (node.right != null) {
                    q1.offer(node.right);
                }
            }
            if (!skip) {
                count++;
            } else {
                skip = false;
            }
            q0.addAll(q1);
            q1.clear();
        }
	}
    
    /**
     *
     @return
     */
    public KDTreeNodeFloat getRoot() {
		return root;
	}

    private double distanceSq(KDTreeNodeFloat tree, 
        float leftValue, float rightValue) {

        float diffX = tree.getX() - leftValue;
        float diffY = tree.getY() - rightValue;
        
        return (diffX * diffX) + (diffY * diffY);
    }
	
}
