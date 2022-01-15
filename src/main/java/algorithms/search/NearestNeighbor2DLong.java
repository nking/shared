package algorithms.search;

import algorithms.heapsAndPQs.YFastTrieLong;
import algorithms.util.ObjectSpaceEstimator;
import algorithms.util.PairInt;
import gnu.trove.iterator.TLongIterator;
import gnu.trove.map.TLongLongMap;
import gnu.trove.map.hash.TLongLongHashMap;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 * a nearest neighbor's algorithm using XFastTrie
 * for predecessor and successor queries
 * on spatially indexed numbers.
 * 
 * The algorithm performs better on dense data
 * (that is because the base of the prefix tree
 * is filled, leaving smaller number of nodes to
 * create in linear time).
 * The queries depend upon the maximum of x and
 * maximum of y to be entered.
 * 
 * At this time, all entries must be non-negative numbers and
 * so must inquiries.
 * 
 * A worst case query would be when column 0 is
 * filled with points and no others filled elsewhere, 
 * and the last point in the last row and last 
 * column is the query point.
 * In this worst case, the query time would scale
 * roughly as maxY * O(log_2(w))
 
  The algorithm starts with a predecessor and successor 
  call on the query point.  The minimum distance among
  those 2 becomes the goal to search to for completeness
  as rows above and below the query in the same column.
  The search continues to higher rows making predecessor and
  successor calls until the goal is reached.  The next higher
  row is one less than the predecessor result.
  The goal to the complete search is reduced by smaller distance
  answers.  The same is repeated for lower rows.
  <pre>
  an example would be:
  
     0   1   2   3   4

     5   6   7   8   9

     10  11 *12  13  14     q='18'.  pred='15', succ='null' --> goals(3, 23)
                                     13.pred='12'
    *15  16  17 *18  19              goals change to (13,23)
                                     13.succ='15', not closer than 12.
     20  21  22  23  24              23.pred and 23.succ not closer than 12
                             ans='12'.  queries: 3 pred, 3 succ queries.
                                     at O(log_2(w)) each
                                     complexity was 
                                           6 * O(log_2(w))
                                     for max index = 24, have w = 6 
 </ore>

   first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright

   TODO: consider for the y-axis, adding locality based hashing or similar
 
 * @author nichole
 */
public class NearestNeighbor2DLong {
    
    private final YFastTrieLong xbt;
            
    private final int width;
    
    private final int height;
    
    private final long maxIndex;
    
    private boolean useCache = true;
    
    //TODO: consider replacing these with weak reference key hash maps.
    //  even though the later take more space, they can evict keys no longer
    //  used elsewhere.  the google guava library has several structures 
    //  could use for a cache and there are a couple libraries such as cache 2k.
    private TLongLongMap pCache = new TLongLongHashMap();
    private TLongLongMap sCache = new TLongLongHashMap();
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * 
     * @param points non-negative coordinates
     * @param imgWidth maximum x value of any data point + 1 including
     *    those to be queries
     * @param imgHeight maximum y value of any data point + 1 including
     *    those to be queries
     */
    public NearestNeighbor2DLong(Set<PairInt> points,
        int imgWidth, int imgHeight) {
        
        this.width = imgWidth;
        this.height = imgHeight ;
                
        maxIndex = (long)width * height;
        
        int maxW = 1 + (int)Math.ceil(Math.log(maxIndex)/Math.log(2));
        
        xbt = new YFastTrieLong(maxW);
        
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
                        
            if (x > width || x < 0) {
                throw new IllegalArgumentException(
                    "x cannot be larger than "
                    + " maxX given in constructor " + width
                    + ". x=" + x);
            }

            if (y > height || y < 0) {
                throw new IllegalArgumentException(
                    "y cannot be larger than "
                    + " maxY given in constructor " + height + ". y=" + y);
            }
            
            long index = getInternalIndex(x, y);
            
            //System.out.format("add %d  (%d, %d)\n", index, x, y);
            
            xbt.add(index);
        }
    }
    
    /**
     * 
     * @param pointIdxs pixel indexes formed from relationship
     *   pixIdx = (row * width) + col
     * @param imgWidth maximum x value of any data point  + 1including
     *    those to be queries
     * @param imgHeight maximum y value of any data point  + 1including
     *    those to be queries
     */
    public NearestNeighbor2DLong(TLongSet pointIdxs, 
        int imgWidth, int imgHeight) {
        
        this.width = imgWidth;
        this.height = imgHeight;
                
        maxIndex = (long)width * height;
        
        int maxW = 1 + (int)Math.ceil(Math.log(maxIndex)/Math.log(2));
        
        xbt = new YFastTrieLong(maxW);
        
        TLongIterator iter = pointIdxs.iterator();
        while (iter.hasNext()) {
            
            long pixIdx = iter.next();
            
            int x = getCol(pixIdx);
            int y = getRow(pixIdx);
            
            if (x > width || x < 0) {
                throw new IllegalArgumentException(
                    "x cannot be larger than "
                    + " maxX given in constructor " + width
                    + ". x=" + x);
            }

            if (y > height || y < 0) {
                throw new IllegalArgumentException(
                    "y cannot be larger than "
                    + " maxY given in constructor " + height + ". y=" + y);
            }
            
            xbt.add(pixIdx);
        }
    }
    
    public void doNotUseCache() {
        useCache = false;
    }
    
    protected long getInternalIndex(int col, int row) {
        long t = ((long)width * row) + col;
        
        return t;
    }
    
    protected int getRow(long internalIndex) {
        int row = (int)(internalIndex/width);        
        return row;
    }
    
    protected int getCol(long internalIndex) {
        //int row = (int)(internalIndex/width);
        //int col = (int)(internalIndex - ((long)row * width));
        return (int)(internalIndex % width);
    }
    
    /**
    <pre>
      runtime complexity is
         best case: 
            Note that caching leads to an O(1) term
            over time.
            
         worst case: 
         
         Note, worst case is: first column
         filled with points and all else is empty and
         the number of rows is same or larger than 
         number of columns and the
         query is for the point in the last column and
         last row... a predecessor call is necessary for
         each row in the worst case.
          
     Note: maxW = 1 + Math.ceil(Math.log(maxX * maxY)/Math.log(2));            
     </ore>
    
     * @param x non-negative x coord to query for
     * @param y non-negative y coord to query for
     */
    public Set<PairInt> findClosest(final int x, final int y) {
        
        return findClosestWithinTolerance(x, y, 0);
    }
    
    /**
     * NOTE: NOT READY FOR USE
     * method to return only the nearest point and any
     * that are at the same distance within a tolerance.
     * This is meant to be a nearest neighbor method
     * with ability to return more than one at same distance within a tolerance
     * of that distance.  
     * TODO: calculate the runtime complexity bounds....
     * @param x non-negative x coord to query for
     * @param y non-negative y coord to query for
     * @param tolerance
     * @return 
     */
    public Set<PairInt> findClosestWithinTolerance(int x, int y,
        double tolerance) {
        
        if (xbt.size() == 0) {
            log.warning("xbt is empty");
            return new HashSet<PairInt>();
        }
        
        if (x >= width || x < 0) {
            //throw new IllegalArgumentException(
            log.fine(
            "x cannot be larger than "
                + " maxX given in constructor " + width
                + ". x=" + x);
            return null;
        }
        
        if (y >= height || y < 0) {
            //throw new IllegalArgumentException(
            log.fine(
                "y cannot be larger than "
                + " maxY given in constructor " + height + ". y=" + y);
            return null;
        }
        
        double closestDist = Double.MAX_VALUE;
        double closestDistPlusTol = Double.MAX_VALUE;
        
        TLongSet closestIndexes = new TLongHashSet();
        
        long idx = getInternalIndex(x, y);
        
        //System.out.format("find  %d  (=%d, %d)\n", idx, x, y);
        
        //O(1)
        long q = xbt.find(idx);
        if (q != -1) {
            // have found nearest, but still need to search
            // within tolerance distance for others.
            closestDist = 0;
            closestDistPlusTol = tolerance;
            closestIndexes.add(idx);
        }
                
        long predecessor = -1;
        long successor = -1;
        
        if (useCache && pCache.containsKey(idx)) {
            predecessor = pCache.get(idx);
        } else {
            //O(log_2(maxW))
            predecessor = xbt.predecessor(idx);
            if (useCache && predecessor != -1) {
                pCache.put(idx, predecessor);
            }
        }
        if (useCache && sCache.containsKey(idx)) {
            successor = sCache.get(idx);
        } else {
            //O(log_2(maxW))
            successor = xbt.successor(idx);
            if (useCache && successor != -1){
                sCache.put(idx, successor);
            }
        }
        
        double dp2 = dist(x, y, predecessor);
        double ds2 = dist(x, y, successor);
        double dMin = Math.min(dp2, ds2);
        
        //System.out.println("p=" + predecessor + " s=" + successor);
 
        /*
        if smallest is smaller than closest 
           if the new closest diff with current is greater 
               than tol, clear the indexes and reset closest 
               vars and add smallest to indexes
               also add the other if within tolerance
           else if closer is within tolerance,
              update closest vars and add whichever or both 
              s2 and p2 to indexes (delaying detailed checks 
              of indexes until end of method)
        else if smallest is <= closestPlusTol
            add s2 and/or p2 to indexes
        */
        
        if (dMin <= closestDist) {
            if (Math.abs(closestDist - dMin) > tolerance) {
                closestIndexes.clear();
                closestDist = dMin;
                closestDistPlusTol = closestDist + tolerance;
            }
            if ((predecessor != -1) && 
                (dp2 <= closestDistPlusTol)) {
                closestIndexes.add(predecessor);
            }
            if ((successor != -1) &&
                (ds2 <= closestDistPlusTol)) {
                closestIndexes.add(successor);
            }
        } else if (dMin <= closestDistPlusTol) {
            if ((predecessor != -1) && 
                (dp2 <= closestDistPlusTol)) {
                closestIndexes.add(predecessor);
            }
            if ((successor != -1) &&
                (ds2 <= closestDistPlusTol)) {
                closestIndexes.add(successor);
            }
        }
        
        //add tolerance to goal
        int goal = (closestDist != Double.MAX_VALUE) ?
            (int)Math.ceil(closestDistPlusTol) : 0;
        
        int yLow = estimateLowBound(y, goal);
       
        int yCurrent;
        if (predecessor == -1) {
            yCurrent = Integer.MIN_VALUE;
        } else {
            int pRow = getRow(predecessor);
            if (pRow < y) {
                yCurrent = pRow;
            } else {
                yCurrent = pRow - 1;
            }
        }
        
        //System.out.println("yCurrent=" + yCurrent + " yLow=" + yLow);
        
        // predecessor searches until reach yLow, adjusting goal by
        //   min distances
        long p2 = -1; 
        long s2 = -1;
        while (yCurrent >= yLow) {
            long cIdx = getInternalIndex(x, yCurrent);

            //O(1)
            q = xbt.find(cIdx);
            if (q != -1) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                if (useCache && pCache.containsKey(cIdx)) {
                    p2 = pCache.get(cIdx);
                } else {
                    p2 = xbt.predecessor(cIdx);
                    if (useCache && p2 != -1) {
                        pCache.put(cIdx, p2);
                    }
                }
                if (useCache && sCache.containsKey(cIdx)) {
                    s2 = sCache.get(cIdx);
                } else {
                    //O(log_2(maxW))
                    s2 = xbt.successor(cIdx);
                    if (useCache && s2 != -1) {
                        sCache.put(cIdx, s2);
                    }
                }
                
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
        
            dMin = Math.min(dp2, ds2);
            if (dMin <= closestDist) {
                if (Math.abs(closestDist - dMin) > tolerance) {
                    closestIndexes.clear();
                    closestDist = dMin;
                    closestDistPlusTol = closestDist + tolerance;
                    goal = (int)Math.ceil(closestDistPlusTol);
                    yLow = estimateLowBound(y, goal); 
                }
                if ((p2 != -1) && dp2 <= closestDistPlusTol) {
                    closestIndexes.add(p2);
                }
                if ((s2 != -1) && ds2 <= closestDistPlusTol) {
                    closestIndexes.add(s2);
                }
            } else if (dMin <= closestDistPlusTol) {
                if ((p2 != -1) && dp2 <= closestDistPlusTol) {
                    closestIndexes.add(p2);
                }
                if ((s2 != -1) && ds2 <= closestDistPlusTol) {
                    closestIndexes.add(s2);
                }
            }    
            
            if (p2 != -1) {
                long expectedNext = getInternalIndex(x, yCurrent - 1);
                if (p2 > expectedNext) {
                    yCurrent -= 1;
                } else {
                    yCurrent = getRow(p2) - 1;
                }
            } else {
                yCurrent = Integer.MIN_VALUE;
            }
        }
       
        //System.out.println("yCurrent=" + yCurrent + " yLow=" + yLow);
        //System.out.println("p=" + p2 + " s=" + s2);
        
        // successor searches to higher bounds
        if (successor == -1) {
            yCurrent = Integer.MAX_VALUE;
        } else {
            int sr = getRow(successor);
            if (sr > y) {
                yCurrent = sr;
            } else {
                yCurrent = sr + 1;
            }
        }
        int yHigh = estimateHighBound(y, goal);
        
        //System.out.println("yCurrent=" + yCurrent + " yHigh=" + yHigh);
        
        while (yCurrent <= yHigh) {
            long cIdx = getInternalIndex(x, yCurrent);
            
            //O(1)
            q = xbt.find(cIdx);
            if (q != -1) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                if (useCache && pCache.containsKey(cIdx)) {
                    p2 = pCache.get(cIdx);
                } else {
                    //O(log_2(maxW))
                    p2 = xbt.predecessor(cIdx);
                    if (useCache && p2 != -1) {
                        pCache.put(cIdx, p2);
                    }
                }
                if (useCache && sCache.containsKey(cIdx)) {
                    s2 = sCache.get(cIdx);
                } else {
                    //O(log_2(maxW))
                    s2 = xbt.successor(cIdx);
                    if (useCache && s2 != -1) {
                        sCache.put(cIdx, s2);
                    }
                }
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
            
            dMin = Math.min(dp2, ds2);
            if (dMin <= closestDist) {
                if (Math.abs(closestDist - dMin) > tolerance) {
                    closestIndexes.clear();
                    closestDist = dMin;
                    closestDistPlusTol = closestDist + tolerance;
                    goal = (int)Math.ceil(closestDistPlusTol);
                    yHigh = estimateHighBound(y, goal); 
                }
                if ((p2 != -1) && dp2 <= closestDistPlusTol) {
                    closestIndexes.add(p2);
                }
                if ((s2 != -1) && ds2 <= closestDistPlusTol) {
                    closestIndexes.add(s2);
                }
            } else if (dMin <= closestDistPlusTol) {
                if ((p2 != -1) && dp2 <= closestDistPlusTol) {
                    closestIndexes.add(p2);
                }
                if ((s2 != -1) && ds2 <= closestDistPlusTol) {
                    closestIndexes.add(s2);
                }
            }    
            
            if (s2 != -1) {
                long expectedNext = getInternalIndex(x, yCurrent + 1);
                if (s2 < expectedNext) {
                    yCurrent += 1;
                } else {
                    yCurrent = getRow(s2) + 1;
                }
            } else {
                yCurrent = Integer.MAX_VALUE;
            }
            
            //System.out.println("yCurrent=" + yCurrent + " yHigh=" + yHigh);
        }
        
        //filter results for closest and tolerance
        Set<PairInt> results = new HashSet<PairInt>();
        TLongIterator iter = closestIndexes.iterator();
        while (iter.hasNext()) {
            long index2 = iter.next();
            if (dist(x, y, index2) <= closestDistPlusTol) {
                int x2 = getCol(index2);
                int y2 = getRow(index2);
                PairInt p3 = new PairInt(x2, y2);
                results.add(p3);
            }
        }
 
        return results;
    }
       
    /**
    <pre>
      runtime complexity is
         best case: 2 * O(log_2(maxW)).
            Note that caching leads to an O(1) term
            over time instead of the logarithmic term.
            
         worst case: nRows * 2 * O(log_2(maxW))
         
         Note, worst case is: first column
         filled with points and all else is empty and
         the number of rows is same or larger than 
         number of columns and the
         query is for the point in the last column and
         last row... a predecessor call is necessary for
         each row in the worst case.
          
     Note: maxW = 1 + Math.ceil(Math.log(maxX * maxY)/Math.log(2));            
     </ore>
    
     * @param x non-negative x coord to query for
     * @param y non-negative y coord to query for
     */
    public Set<PairInt> findClosestNotEqual(final int x, final int y) {
        
        return findClosest(x, y, Integer.MAX_VALUE, false);
    }
    
    /**
    <pre>
      runtime complexity is
         best case: 
            Note that caching leads to an O(1) term
            over time.
            
         worst case: 
         
      Note: maxW = 1 + Math.ceil(Math.log(maxX * maxY)/Math.log(2));
     </ore>
    
     * @param x
     * @param y
     * @param dMax
     * @return a set of points within dMax that are the 
     * closest points, else returns an empty set
     */
    public Set<PairInt> findClosest(int x, int y, int dMax) {
        return findClosest(x, y, dMax, true);
    }
    
    /**
    <pre>
      runtime complexity is
         best case: 
            Note that caching leads to an O(1) term
            over time.
            
         worst case: 
         
      Note: maxW = 1 + Math.ceil(Math.log(maxX * maxY)/Math.log(2));
     </ore>
    
     * @param x
     * @param y
     * @param dMax
     * @return a set of points within dMax that are the 
     * closest points, else returns an empty set
     */
    private Set<PairInt> findClosest(int x, int y, int dMax, boolean includeEquals) {
        
        if (x >= width || x < 0) {
            //throw new IllegalArgumentException(
            log.fine(
            "x cannot be larger than "
                + " maxX given in constructor " + width
                + ". x=" + x);
            return null;
        }
        
        if (y >= height || y < 0) {
            //throw new IllegalArgumentException(
            log.fine(
                "y cannot be larger than "
                + " maxY given in constructor " + height + ". y=" + y);
            return null;
        }
        
        long idx = getInternalIndex(x, y);
        
        if (includeEquals) {
            long q = xbt.find(idx);
            if (q != -1) {
                Set<PairInt> results = new HashSet<PairInt>();
                results.add(new PairInt(x, y));
                return results;
            }
        }
                
        TLongSet closestIndexes = new TLongHashSet();
        
        double closestDist = Double.MAX_VALUE;
        
        long predecessor = -1;
        long successor = -1;
        
        if (useCache && pCache.containsKey(idx)) {
            predecessor = pCache.get(idx);
        } else {
            //O(log_2(maxW))
            predecessor = xbt.predecessor(idx);
            if (useCache && predecessor != -1) {
                pCache.put(idx, predecessor);
            }
        }
        if (useCache && sCache.containsKey(idx)) {
            successor = sCache.get(idx);
        } else {
            //O(log_2(maxW))
            successor = xbt.successor(idx);
            if (useCache && successor != -1) {
                sCache.put(idx, successor);
            }
        }
        
        double dp2 = dist(x, y, predecessor);
        double ds2 = dist(x, y, successor);
        if (!includeEquals) {
            if (dp2 == 0) {
                ds2 = Double.MAX_VALUE;
            } else if (ds2 == 0) {
                ds2 = Double.MAX_VALUE;
            }
        }
        if (dp2 <= ds2 && (dp2 <= dMax)) {
            closestDist = dp2;
            closestIndexes.add(predecessor);
            if (dp2 == ds2) {
                closestIndexes.add(successor);
            }
        } else if (ds2 < dp2 && (ds2 <= dMax)) {
            closestDist = ds2;
            closestIndexes.add(successor);
        }
        
        int goal = (closestDist != Double.MAX_VALUE) ?
            (int)Math.ceil(closestDist) : dMax;
        
        if (goal > dMax) {
            goal = dMax;
        }
        
        int yLow = estimateLowBound(y, goal);
       
        int yCurrent;
        if (predecessor == -1) {
            yCurrent = Integer.MIN_VALUE;
        } else {
            int pRow = getRow(predecessor);
            if (pRow < y) {
                yCurrent = pRow;
            } else {
                yCurrent = pRow - 1;
            }
        }
        
        // predecessor searches until reach yLow, adjusting goal by
        //   min distances
        long p2 = -1; 
        long s2 = -1;
        while (yCurrent >= yLow) {
            long cIdx = getInternalIndex(x, yCurrent);
            
            long q = xbt.find(cIdx);
            if (q != -1) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                if (useCache && pCache.containsKey(cIdx)) {
                    p2 = pCache.get(cIdx);
                } else {
                    //O(log_2(maxW))
                    p2 = xbt.predecessor(cIdx);
                    if (useCache && p2 != -1) {
                        pCache.put(cIdx, p2);
                    }
                }
                if (useCache && sCache.containsKey(cIdx)) {
                    s2 = sCache.get(cIdx);
                } else {
                    //O(log_2(maxW))
                    s2 = xbt.successor(cIdx);
                    if (useCache && s2 != -1) {
                        sCache.put(cIdx, s2);
                    }
                }
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
            if (!includeEquals) {
                if (s2 != -1 && s2 == idx) {
                    ds2 = Double.MAX_VALUE;
                }
            }
            if ((dp2 < ds2) && (dp2 < closestDist) && (dp2 <= dMax)) {
                closestIndexes.clear();
                closestDist = dp2;
                closestIndexes.add(p2);
                goal = (int)Math.ceil(closestDist);
                if (goal > dMax) {
                    goal = dMax;
                }
                yLow = estimateLowBound(y, goal);                
            } else if ((ds2 < dp2) && (ds2 < closestDist) && (ds2 <= dMax)) {
                           
                closestIndexes.clear();
                closestDist = ds2;
                closestIndexes.add(s2);
                goal = (int)Math.ceil(closestDist);
                if (goal > dMax) {
                    goal = dMax;
                }
                yLow = estimateLowBound(y, goal);
            } else if (dp2 == closestDist && (dp2 != Double.MAX_VALUE)
                && (dp2 <= dMax)) {               
                closestIndexes.add(p2);
                if (dp2 == ds2) {
                    closestIndexes.add(s2);
                }
            } else if (ds2 == closestDist && (ds2 != Double.MAX_VALUE)
                && (ds2 <= dMax)) {                
                closestIndexes.add(s2);
            }
            
            if (p2 != -1) {
                long expectedNext = getInternalIndex(x, yCurrent - 1);
                if (p2 > expectedNext) {
                    yCurrent -= 1;
                } else {
                    yCurrent = getRow(p2) - 1;
                }
            } else {
                yCurrent = Integer.MIN_VALUE;
            }
        }
        
        // successor searches to higher bounds
        if (successor == -1) {
            yCurrent = Integer.MAX_VALUE;
        } else {
            int sr = getRow(successor);
            if (sr > y) {
                yCurrent = sr;
            } else {
                yCurrent = sr + 1;
            }
        }
        int yHigh = estimateHighBound(y, goal);
        
        while (yCurrent <= yHigh) {
            long cIdx = getInternalIndex(x, yCurrent);
            long q = xbt.find(cIdx);
            if (q != -1) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                if (useCache && pCache.containsKey(cIdx)) {
                    p2 = pCache.get(cIdx);
                } else {
                    //O(log_2(maxW))
                    p2 = xbt.predecessor(cIdx);
                    if (useCache && p2 != -1) {
                        pCache.put(cIdx, p2);
                    }
                }
                if (useCache && sCache.containsKey(cIdx)) {
                    s2 = sCache.get(cIdx);
                } else {
                    //O(log_2(maxW))
                    s2 = xbt.successor(cIdx);
                    if (useCache && s2 != -1) {
                        sCache.put(cIdx, s2);
                    }
                }
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
            if (!includeEquals) {
                if (p2 != -1 && p2 == idx) {
                    dp2 = Double.MAX_VALUE;
                } else if (s2 != -1 && s2 == idx) {
                    ds2 = Double.MAX_VALUE;
                }
            }
            if ((dp2 < ds2) && (dp2 < closestDist) && (dp2 <= dMax)) {
                closestIndexes.clear();
                closestDist = dp2;                
                closestIndexes.add(p2);
                goal = (int)Math.ceil(closestDist);
                if (goal > dMax) {
                    goal = dMax;
                }
                yHigh = estimateHighBound(y, goal);
            } else if ((ds2 < dp2) && (ds2 < closestDist) && (ds2 <= dMax)) {                
                closestIndexes.clear();
                closestDist = ds2;
                closestIndexes.add(s2);
                goal = (int)Math.ceil(closestDist);
                if (goal > dMax) {
                    goal = dMax;
                }
                yHigh = estimateHighBound(y, goal);
            } else if (dp2 == closestDist && (dp2 != Double.MAX_VALUE)
                && (dp2 <= dMax)) {                
                closestIndexes.add(p2);
                if (dp2 == ds2) {                    
                    closestIndexes.add(s2);
                }
            } else if (ds2 == closestDist && (ds2 != Double.MAX_VALUE)
                && (ds2 <= dMax)) {                
                closestIndexes.add(s2);
            } 
            
            if (s2 != -1) {
                long expectedNext = getInternalIndex(x, yCurrent + 1);
                if (s2 < expectedNext) {
                    yCurrent += 1;
                } else {
                    yCurrent = getRow(s2) + 1;
                }
            } else {
                yCurrent = Integer.MAX_VALUE;
            }
        }
        
        Set<PairInt> results = new HashSet<PairInt>();
        TLongIterator iter = closestIndexes.iterator();
        while (iter.hasNext()) {
            long index2 = iter.next();
            int x2 = getCol(index2);
            int y2 = getRow(index2);
            results.add(new PairInt(x2, y2));
        }
        
        return results;
    }

    private double dist(int x, int y, long p2) {
        
        if (p2 == -1) {
            return Double.MAX_VALUE;
        }
        
        int x2 = getCol(p2);
        int y2 = getRow(p2);
        
        int diffX = x2 - x;
        int diffY = y2 - y;
        double dist = Math.sqrt(diffX * diffX + diffY * diffY);
        
        return dist;
    }

    private int estimateLowBound(int y, int goal) {

        int low = y - goal;
        
        if (low < 0) {
            low = 0;
        }
        
        return low;
    }
    
    private int estimateHighBound(int y, int goal) {

        int high = y + goal;
        
        if (high > height) {
            high = height;
        }
        
        return high;
    }
 
    public static long estimateSizeOnHeap(int numberOfPoints,
        int maxBitLength) {
         
        long[] yTotals = YFastTrieLong.estimateSizeOnHeap(numberOfPoints, 
            maxBitLength);
        
        //System.out.println("yft estimates=" + Arrays.toString(yTotals));
        
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNObjRefsFields(4);
        est.setNIntFields(2);
        est.setNLongFields(1);
        est.setNBooleanFields(1);
        
        long total = est.estimateSizeOnHeap();
       
        total += yTotals[1];
                
        return total;
    }
}
