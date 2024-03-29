package thirdparty.edu.princeton.cs.algs4;

import java.util.ArrayList;
import java.util.List;

/******************************************************************************
 A specialization of RangeSearch for an Interval parameterized type.
 
 from RangeSearch in algs4.jar
   from the book "Algorithms" by Sedgewick and Wayne
 * http://algs4.cs.princeton.edu/
 * copyright for authors Robert Sedgewick and Kevin Wayne
 * is GPLV3, http://algs4.cs.princeton.edu/faq/
   (see bottom of this file)
 
   first adapted in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)
   then moved to this shared library project which has the same copyright
   and license.
 
 NOTE that the intervals cannot overlap.  If a put of an interval
 intersects with a key, the existing interval in the tree gets
 the value of the new interval, but the key range does not
 change.
  
 @param <T> the data type used in the Intervals
 @param <Value> the data type of the key associated with each
 * tree interval.
 ******************************************************************************/
public class IntervalRangeSearch<T extends Comparable<T>, Value> extends
    RangeSearch<Interval<T>, Value> {
    
    /**
     *
     @param interval
     @return
     */
    public Queue<Interval<T>> range0(Interval<T> interval) {
        Queue<Interval<T>> list = new Queue<Interval<T>>();
        //System.out.println("root=" + root);
        //System.out.println("srch interval=" + interval);
        range0(root, interval, list);
        return list;
    }
    
    /**
     *
     @return
     */
    public List<Interval<T>> getAllIntervals() {
        List<Interval<T>> list = new ArrayList<Interval<T>>();
        IntervalRangeSearch.this.getAllIntervals(root, list);
        return list;
    }
    
    /**
     *
     @param outputIntervals
     @param outputValues
     */
    public void getAllIntervals(List<Interval<T>> outputIntervals, 
        List<Value> outputValues) {
        getAllIntervals(root, outputIntervals, outputValues);
    }
    
    private void getAllIntervals(RangeSearchNode<Interval<T>, Value> x, 
        List<Interval<T>> outputIntervals, List<Value> outputValues) {
        
        if (x == null) {
            return;
        }
        outputIntervals.add(x.key);
        outputValues.add(x.val);
            
        if (x.left != null) {
            getAllIntervals(x.left, outputIntervals, outputValues);
        }
        
        if (x.right != null) {
            getAllIntervals(x.right, outputIntervals, outputValues);
        }
    }
    
    private void getAllIntervals(RangeSearchNode<Interval<T>, Value> x, 
        List<Interval<T>> list) {
        
        if (x == null) {
            return;
        }
        list.add(x.key);
            
        if (x.left != null) {
            getAllIntervals(x.left, list);
        }
        
        if (x.right != null) {
            getAllIntervals(x.right, list);
        }
    }
    
    private void range0(RangeSearchNode<Interval<T>, Value> x, 
        Interval<T> interval, Queue<Interval<T>> list) {
       
        if (x == null) return;
       
        boolean intersects = interval.intersects(x.key);
        if (intersects) {
            list.enqueue(x.key);
        }
        
        /*
        interval has min max for search.
        
        tree has left as the larger keys
        
                     x
              lft         rgt
           lft  rgt     lft  rgt
        
        
        or viewed by increasing values--->
               xmin--x--xmax smin
              rgt         lft
           rgt  lft     rgt  lft
        */
          
        // if x.max < interval.min
        //   search lft
        // if xmin > interval.max
        //   search rgt
        
        if ((x.left != null) && (intersects || 
            (x.key.max().compareTo(interval.min()) < 1)) ) {
            range0(x.left, interval, list);
        }
        
        if ((x.right != null) && (intersects || 
            (interval.max().compareTo(x.key.min()) < 1))) {
            range0(x.right, interval, list);
        }
    }

}
