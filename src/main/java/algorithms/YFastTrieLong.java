package algorithms;

import gnu.trove.map.TLongLongMap;
import gnu.trove.map.TLongObjectMap;
import gnu.trove.map.hash.TLongLongHashMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import thirdparty.edu.princeton.cs.algs4.RedBlackBSTLongInt;
import thirdparty.ods.Longizer;
import thirdparty.ods.XFastTrieLong;
import thirdparty.ods.XFastTrieNodeLong;

/**
 * NOTE: has not been tuned for best runtime complexity yet.  
 * 
 * from wikipedia
 *     https://en.wikipedia.org/wiki/Y-fast_trie
 *  
 * a y-fast trie is a data structure for storing 
 * integers from a bounded domain. It supports exact and predecessor 
 * or successor queries in time O(log log M), using O(n) space, 
 * where n is the number of stored values and M is the maximum 
 * value in the domain. 
 * The structure was proposed by Dan Willard in 1982[1] to decrease 
 * the O(n log M) space used by an x-fast trie.
   
   The Y-Fast trie has the ordered associative array operations + successor and
   predecessor.
   
   NOTE that the runtime complexities listed are not yet achieved.
   The current operations depend upon settings, but should be 
   constant runtime complexity .lte. O(10).

   Find(k): find the value associated with the given key.
       runtime complexity is O(log log(M))
   Successor(k): find the key/value pair with the smallest key larger than or 
       equal to the given key.
       runtime complexity is O(log log(M))
   Predecessor(k): find the key/value pair with the largest key less than or 
       equal to the given key.
       runtime complexity is O(log log(M))
   Insert(k, v): insert the given key/value pair.
       runtime complexity is O(log log(M))
   Delete(k): remove the key/value pair with the given key.
       runtime complexity is O(log log(M))
 
 * Note, have not read the Willard paper yet, just a few online
 * lecture notes to implement this.  
 *  
    first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright

 * @author nichole
 */
public class YFastTrieLong {

    /*    
    designing from browsing a few different lecture notes
    online. the yfast trie uses same w and maxC as
    the XFastTrie.
      - creates w red black trees to hold inserted heap nodes.
        the w trees each have range size of maxC/w and
        start from 0 extending to last one holding maxC.
      - each tree has a representative if it has any nodes
        and those are stored
        in the XFastTrie of this YFastTrie.
      - because the XFastTrie only holds w xft values,
        the space complexity is reduced.
-------------------------------
YFastTrie

   - w bits set by maximum expected value to be added.
   - one XFastTrie to hold the representives (at most w in number)
   - w red black trees to keep ordered points.
     - because some of the items added may have more than
       one with same key value, the values in the red black tree
       will be linked lists.

    NOTE: topics to consider for improvements:
          the distribution of rb trees, that is their partitions,
          could be improved dynamically.
          For example, if maxC were value 127, but the majority
          of nodes at some point in time were in bin 0 at values
          near 4, one would prefer to divide that tree 
          into more than one tree to speed up searches.
          This begins to look like a good reason to
          compare to multi-level-buckets.  The only implementation
          I could find was the Andrew Goldberg MLB offered 
          under a license that is freely available for non-commercial
          use, else need to contact for permission... not wanting
          to include mixed license restrictions for now...
          (so I didn't download and read the code.  am reading
          his 2 papers on the subject, but they depend upon other
          papers too, so gathering all the specs for the MLB
          algorithms is not complete...)
          -- one possible work around without making dynamic
          partitions in the YFastTrie would be to know or 
          estimate the population of data ahead of time and 
          then make separate YFastTrie's for manually partitioned 
          data (changing zero-points, a.k.a. bias levels as needed
          before and after use of more than one YFastTrie)
    */
    
    private int n = 0;
    
    private final int w;
    
    private final long maxC;
    
    private final long binSz;
    
    private long nBins;

    private final XFastTrieLong<XFastTrieNodeLong<Long>, Long> xft;
    
    // key = bin index (which is node/binSz), value = repr value.
    // each repr value is the minimum stored in the bin.
    private final TLongLongMap xftReps = new TLongLongHashMap();
    
    // there are w items in rbs
    // each list item is a sorted binary search tree of numbers in that bin.
    //    the value is the tree holds the number of times that number is present.
    // hashkey= node/binSz;     treemap key=node, value=multiplicity
    private final TLongObjectMap<RedBlackBSTLongInt> rbs;

    private boolean chooseByN = true;

    public YFastTrieLong(int wBits) {
        
        if (wBits < 63 && wBits > 1) {
            this.w = wBits;
        } else {
            throw new IllegalStateException("wBits "
                + " shoulw be greater than 1 and less than 32");
        }
        maxC = (1L << w) - 1;
        
        /*
        LG binsize = (int)Math.ceil((float)maxC/(float)w);
        MID binsize = 10
        
        w,     binsz,     n,             rt
        31,    69273666,  31,            26.   LG
        31,    5.0,       429496729.6,   2.32  MID
        
        10,    102,       10,            6.67  LG
        10,    4.0,       256.0,         2.0   MID
        
        aiming for a runtime of about O(10) or better without increasing n too
        much.
        
        alternatively, could keep n as large as possible within good
        performance range, hence the binsz will be small, hence the
        runtime will be small.
        
        n            TreeMaps
        n * binSz    objects in TreeMaps
        
        heapsize: 100's or more MB
        */
        
        if (chooseByN) {
            binSz = chooseBinSizeByN();
        } else {
            long tmpLg = (long)Math.ceil((double) maxC / (double) w);
            if ((Math.log(tmpLg) / Math.log(2)) < 10) {
                binSz = tmpLg;
            } else {
                binSz = 1024;
            }
        }
        
        nBins = (long)Math.ceil((double)maxC/(double)binSz);
    
        System.out.println("nBins=" + nBins + "  rt of ops=" +
            (Math.log(binSz)/Math.log(2)));
        
        rbs = new TLongObjectHashMap<RedBlackBSTLongInt>();
         
        XFastTrieNodeLong<Long> clsNode = new XFastTrieNodeLong<Long>();
        Longizer<Long> it = new Longizer<Long>() {
            @Override
            public long longValue(Long x) {
                return x;
            }
        };
        
        xft = new XFastTrieLong<XFastTrieNodeLong<Long>, Long>(clsNode, it, w);
    }
    
    public YFastTrieLong() {
        
        this.w = 62;
        
        maxC = (1L << w) - 1;
        
        if (chooseByN) {
            binSz = chooseBinSizeByN();
        } else {
            long tmpLg = (long) Math.ceil((double) maxC / (double) w);
            if ((Math.log(tmpLg) / Math.log(2)) < 10) {
                binSz = tmpLg;
            } else {
                binSz = 1024;
            }
        }
        
        nBins = (long)Math.ceil((double)maxC/(double)binSz);
        
        System.out.println("nBins=" + nBins + "  rt of ops=" +
            (Math.log(binSz)/Math.log(2)));
        
        rbs = new TLongObjectHashMap<RedBlackBSTLongInt>();
        
        XFastTrieNodeLong<Long> clsNode = new XFastTrieNodeLong<Long>();
        Longizer<Long> it = new Longizer<Long>() {
            @Override
            public long longValue(Long x) {
                return x;
            }
        };
        
        xft = new XFastTrieLong<XFastTrieNodeLong<Long>, Long>(clsNode, it, w);
    }
    
    protected RedBlackBSTLongInt getTreeMap(long index) {
        RedBlackBSTLongInt tree = rbs.get(index);
        if (tree == null) {
            tree = new RedBlackBSTLongInt();
            rbs.put(index, tree);
        }
        return tree;
    }

    /**
     * 
     * @param node
     * @param index 
     */
    private void addToRBTree(long node, long index) {
        
        RedBlackBSTLongInt tree = getTreeMap(index);
        
        assert(tree != null);
                
        int[] output = new int[2];
        
        tree.get(node, output);
    
        int multiplicity;
        if (output[0] == -1) {
            multiplicity = 1;
        } else {
            multiplicity = 1 + output[1];
        }
        
        tree.put(node, multiplicity);        
    }
    
    /**
     * 
     * @param node
     * @param index 
     */
    private boolean deleteFromRBTree(long node, long index) {
                
        RedBlackBSTLongInt tree = getTreeMap(index);
        
        assert(tree != null);
        
        int[] output = new int[2];
        
        tree.get(node, output);
    
        if (output[0] == -1) {
            return false;
        }
        
        int multiplicity = output[1];
        if (multiplicity > 0) {
            multiplicity = output[1] - 1;
            if (multiplicity > 0) {
                tree.put(node, multiplicity);
            }
        }
        if (multiplicity == 0) {
            tree.delete(node);
        }
        
        return true;
    }
    
    /**
     * runtime complexity is roughly .lte. O(10) and may be better than this
     * for some datasets.  The runtime is dependent on the bit length of the
     * largest number to hold or query, and the balance between the number
     * of bins and number of items per bin.
     * 
     * @param node a number >= 0 and having bit length 
     * less than or equal to w.
     * @return 
     */
    public boolean add(long node) {

        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node.key must "
                + "be less than " + maxC + " node=" + node);
        }
        
        long index = node/binSz;
        
        long existingRepr = xftReps.get(index);
                
        if (!xftReps.containsKey(index)) {
            // insert is O(log_2(w)) + O(l-w)
            xft.add(Long.valueOf(node));
            xftReps.put(index, node);
        } else if (node < existingRepr) {
            // delete is O(log_2(w)) + O(l-w)
            // insert is O(log_2(w)) + O(l-w)
            xft.remove(Long.valueOf(existingRepr));
            xft.add(Long.valueOf(node));
            xftReps.put(index, node);
        }
                
        addToRBTree(node, index);
        
        n++;
        
        return true;
    }

    /**
     * runtime complexity is roughly .lte. O(10) and may be better than this
     * for some datasets.  The runtime is dependent on the bit length of the
     * largest number to hold or query, and the balance between the number
     * of bins and number of items per bin.
     * 
     * @param node
     * @return 
     */
    public boolean remove(long node) {
        
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC);
        }
        
        long index = node/binSz;
                
        boolean removed = deleteFromRBTree(node, index);
                
        if (!removed) {
            return false;
        }
        
        if (!xftReps.containsKey(index)) {
            return false;
        }
        
        RedBlackBSTLongInt tree = getTreeMap(index);
      
        long existingRepr = xftReps.get(index);
      
        if (tree.isEmpty()) {
            // just deleted the last item so remove from rbs
            // delete is O(log_2(w)) + O(w-l)
            if (xftReps.containsKey(index)) {
                xft.remove(Long.valueOf(existingRepr));
                xftReps.remove(index);
            }
        } else if (node == existingRepr) {
            
            //existingRepr is maintained as the minimum in the bin,
            //   so if a node w/ this value is removed and the multiplicity
            //      was 1, need to assign a new repr
            
            int[] output = new int[2];
            tree.get(node, output);
            
            int multiplicity = output[1];
    
            if (output[0] == -1) {
                // remove the current repr and assign a new one
                // delete is O(log_2(w)) + O(w-l)
                xft.remove(Long.valueOf(existingRepr));
                xftReps.remove(index);
            
                // O(log_2(N/w))
                long minKey = tree.min();
                xft.add(minKey);
                xftReps.put(index, minKey); 
            }            
        }
        
        n--;
        
        return true;
    }

    /**
     * runtime complexity is roughly .lte. O(10) and may be better than this
     * for some datasets.  The runtime is dependent on the bit length of the
     * largest number to hold or query, and the balance between the number
     * of bins and number of items per bin.
     * 
     * @param node
     * @return 
     */
    public long find(long node) {
                
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC + ". node=" + node);
        }
                
        long index = node/binSz;
                
        RedBlackBSTLongInt tree = getTreeMap(index);
        
        int[] output = new int[2];
        tree.get(node, output);
        
        if (output[0] == -1) {
            return -1;
        }
        
        return node;
    }

    /**
     * runtime complexity is roughly .lte. O(10) and may be better than this
     * for some datasets.  The runtime is dependent on the bit length of the
     * largest number to hold or query, and the balance between the number
     * of bins and number of items per bin.
     * 
     * @param node
     * @return value preceding node, else -1 if there is not one
     */
    public long predecessor(long node) {
    
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC);
        }
                
        long nodeIndex = node/binSz;
        
        // the repr is stored in xft and it is always the minium for the bin
        boolean isAMinimum = xft.find(Long.valueOf(node)) != null;
        
        /*
        if the node is not a minima, the answer is in
           the node's map if its size is larger > 1
        */
        
        RedBlackBSTLongInt tree = getTreeMap(nodeIndex);
        
        if (!isAMinimum && (tree.size() > 1)) {
            long[] output = new long[2];
            tree.lower(node, output);
            if (output[0] != -1) {
                return output[1];
            }
        }
       
        // else, predeccessor is in the closest bin < nodeIndex that has
        //    items in it.
                
        Long prev = xft.predecessor(Long.valueOf(node));
        if (prev == null) {
            return -1;
        }
        
        long prev0Index = prev.longValue()/binSz;
            
        tree = getTreeMap(prev0Index);
        
        if (tree.isEmpty()) {
            return -1;
        }
        
        long lastKey = tree.max();
        
        return lastKey;
    }
    
    /**
     * runtime complexity is roughly .lte. O(10) and may be better than this
     * for some datasets.  The runtime is dependent on the bit length of the
     * largest number to hold or query, and the balance between the number
     * of bins and number of items per bin.
     * 
     * @param node
     * @return 
     */
    public long successor(long node) {
                
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC);
        }
        
        Long nodeKey = Long.valueOf(node);
        
        long nodeIndex = node/binSz;
        
        boolean isAMinimum = xft.find(nodeKey) != null;
        
        RedBlackBSTLongInt tree = getTreeMap(nodeIndex);
        
        if (isAMinimum) {
            // if tree size > 1, the next key is the successor
            // else, the xft sucessor to nodeIndex is the successor
            
            if (tree.size() > 1) {
                long[] output = new long[2];
                tree.higher(node, output);
                assert(output[0] != -1);
                return output[1];
            }
            
            Long successorRepr = xft.successor(nodeKey);
            if (successorRepr == null) {
                return -1;
            }
            
            // the successor representative is then the next value
            return successorRepr;
        }
        
        // else, the node is not a repr
        //   if there is a tree successor to the node, that is the successor
        //   else, the xft successor to nodeIndex is the successor
                    
        long[] output = new long[2];
        tree.higher(node, output);
        
        if (output[0] != -1) {
            return output[1];
        }
        
        Long successorRepr = xft.successor(nodeKey);
        if (successorRepr == null) {
            return -1;
        }

        // the successor representative is then the next value
        return successorRepr;
    }

    /**
     * runtime complexity is roughly .lte. O(10) and may be better than this
     * for some datasets.  The runtime is dependent on the bit length of the
     * largest number to hold or query, and the balance between the number
     * of bins and number of items per bin.
     * 
     * runtime complexity is O(log_2(w)) 
     * @return minimum, else -1 if empty
     */
    public long minimum() {
        
        if (xft.size() == 0) {
            return -1;
        }
        
        Long repr = xft.minimum();
        
        // cannot be null if size > 0
        assert(repr != null);
       
        return repr.longValue();
    }

    /**
     * runtime complexity is roughly .lte. O(10) and may be better than this
     * for some datasets.  The runtime is dependent on the bit length of the
     * largest number to hold or query, and the balance between the number
     * of bins and number of items per bin.
     * 
     * @return maximum, else -1 if empty
     */
    public long maximum() {
        
        if (xft.size() == 0) {
            return -1;
        }
        
        Long maxRepr = xft.maximum();
        
        assert(maxRepr != null);
        
        long index = maxRepr.longValue()/binSz;
        
        RedBlackBSTLongInt tree = getTreeMap(index);
        
        assert(tree != null);
        assert(!tree.isEmpty());
        
        long lastKey = tree.max();
                
        return lastKey;
    }
    
    /**
     * runtime complexity is roughly .lte. O(10) and may be better than this
     * for some datasets.  The runtime is dependent on the bit length of the
     * largest number to hold or query, and the balance between the number
     * of bins and number of items per bin.
     * 
     * TODO: calc runtime complexity again
     * 
     * @return minimum, else -1 if empty
     */
    public long extractMinimum() {
        
        //O(log_2(w))
        long min = minimum();

        if (min == -1) {
            assert(xft.size() == 0);
            return -1;
        }
                
        remove(min);
        
        return min;
    }
    
    /**
     * runtime complexity is roughly .lte. O(10) and may be better than this
     * for some datasets.  The runtime is dependent on the bit length of the
     * largest number to hold or query, and the balance between the number
     * of bins and number of items per bin.
     * 
     * TODO: calc runtime complexity again
     * 
     * @return maximum, else -1 if empty
     */
    public long extractMaximum() {
        
        long max = maximum();

        if (max == -1) {
            assert(xft.size() == 0);
            return -1;
        }
                
        remove(max);
        
        return max;
    }
    
    public int size() {
        return n;
    }

    private long chooseBinSizeByN() {
        
        long totalMemory = Runtime.getRuntime().totalMemory();
        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        long heapUsage = mbean.getHeapMemoryUsage().getUsed();
        long avail = totalMemory - heapUsage;

        long n = avail/32;
        
        // n = maxC/binsz
        long bs = maxC/n;
        
        double rt = Math.log(bs)/Math.log(2);
        
        if (rt < 10) {
            if (bs > 10) {
                //this is the number of items, that is capacity, of each map
                return bs;
            }
        }
        
        // 2^10 = binSz/2 -->binSz = 2^11
        //return 2048;
        
        // else, fall back to using the default for rt = O(10)
        return (long)Math.ceil((double)maxC/(double)w);
        
        /*
        LG binsize = (int)Math.ceil((float)maxC/(float)w);
        MID binsize = 10
        
        w,     binsz,     n,             rt
        31,    69273666,  31,            26.   LG
        31,    5.0,       429496729.6,   2.32  MID
        
        10,    102,       10,            6.67  LG
        10,    4.0,       256.0,         2.0   MID
        
        aiming for a runtime of about O(10) or better without increasing n too
        much.
        
        alternatively, could keep n as large as possible within good
        performance range, hence the binsz will be small, hence the
        runtime will be small.
        
        n            TreeMaps
        n * binSz    objects in TreeMaps
        
        heapsize: 100's or more MB
        */
        
    }
    
    protected long getBinSz() {
        return binSz;
    }
}
