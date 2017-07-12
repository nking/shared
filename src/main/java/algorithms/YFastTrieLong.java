package algorithms;

import algorithms.util.ObjectSpaceEstimator;
import gnu.trove.map.TLongLongMap;
import gnu.trove.map.TLongObjectMap;
import gnu.trove.map.hash.TLongLongHashMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.util.Arrays;
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
   The current operations depend upon settings to compromise between
   number of secondary hashmaps and the number of objects within those, 
   but the runtime should be 
   constant runtime complexity usually .lte. O(10).

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
   
    add the notes below about the runtime complexities and memory usage here
   
    */
    
    /**
    <pre>
    3 choices:
     AUTOMATIC = code looks at memory available to heap and if there is 
         enough, chooses the fast runtime model of operations O(log_2(w)).
          This is the default.
     FAST_RT = user has specified choice should be the fastest runtime 
          model.  O(log_2(w))
     SPACE_CONSERTVING = user has specified choice should be the smallest 
          memory use model w/ runtimes that are O(log_2(N_entries/w)) which 
          ends up being approx 
          O(20) for nEntries=35e6 and O(12) for nEntries=35e3.
    </pre>
     */
    public static enum BinSizeModel {
        AUTOMATIC, FAST_RT, SPACE_CONSERTVING
    }
    private final BinSizeModel binSzModel;
    
    private int n = 0;
    
    private final int w;
    
    private final long maxC;
    
    private final long binSz;
    
    private long nBins;
    
    /**
     * the minimum of each bin range, if it is populated, is the representative
     * node, and that is held in 2 data structures:
     *    xft holds the number, allowing fast repr prev and next lookups.
     *    xftReps holds the repr as the value, found by key = binNumber.
     */
    private final XFastTrieLong<XFastTrieNodeLong<Long>, Long> xft;
    
    // key = bin index (which is node/binSz), value = repr value.
    // each repr value is the minimum stored in the bin.
    // the max number of map entries will be nBins.
    private final TLongLongMap xftReps = new TLongLongHashMap();
    
    // all inserts of this class are held in w trees in list rbs.
    // each list item is a sorted binary search tree of numbers in that bin.
    //    the value is the tree holds the number of times that number is present.
    // hashkey= node/binSz;     treemap key=node, value=multiplicity
    private final TLongObjectMap<RedBlackBSTLongInt> rbs;

    /**
     * constructor specifying the maximum number of bits of any future add or
     * find, etc, and is by default choosing the model for the fast runtime
     * which may be expensive in space requirements.
     * 
     * @param wBits 
     */
    public YFastTrieLong(int wBits) {
        
        if (wBits < 63 && wBits > 1) {
            this.w = wBits;
        } else {
            throw new IllegalStateException("wBits "
                + " shoulw be greater than 1 and less than 32");
        }
        maxC = (1L << w) - 1;
        
        binSzModel = BinSizeModel.FAST_RT;
                
        binSz = w;
        
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
    
    /**
     * 
     * @param wBits
     * @param numberOfEntries this number of entries to be inserted is used in
     * memory estimates to decide the bin size model.
     */
    public YFastTrieLong(int wBits, int numberOfEntries) {
        
        if (wBits < 63 && wBits > 1) {
            this.w = wBits;
        } else {
            throw new IllegalStateException("wBits "
                + " shoulw be greater than 1 and less than 32");
        }
        maxC = (1L << w) - 1;
        
        binSzModel = BinSizeModel.AUTOMATIC;
               
        binSz = calculateBinSize(numberOfEntries, wBits, binSzModel);
        
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
    
    /**
     * 
     * @param wBits
     * @param binSizeModel 
        <pre>
        3 choices:
         AUTOMATIC = code looks at memory available to heap and if there is 
             enough, chooses the fast runtime model of operations O(log_2(w)).
              This is the default.
         FAST_RT = user has specified choice should be the fastest runtime 
              model.  O(log_2(w))
         SPACE_CONSERTVING = user has specified choice should be the smallest 
              memory use model w/ runtimes that are O(log_2(N_entries/w)) which 
              ends up being approx 
              O(20) for nEntries=35e6 and O(12) for nEntries=35e3.
        </pre>
     */
    public YFastTrieLong(int numberOfEntries, int wBits, 
        BinSizeModel binSizeModel) {
        
        if (wBits < 63 && wBits > 1) {
            this.w = wBits;
        } else {
            throw new IllegalStateException("wBits "
                + " shoulw be greater than 1 and less than 32");
        }
        maxC = (1L << w) - 1;
        
        this.binSzModel = binSizeModel;
               
        binSz = calculateBinSize(numberOfEntries, w, binSizeModel);
        
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
    
    /**
     * constructor using the the maximum number of bits of 62
     * and the default model for the fast runtime
     * which may be expensive in space requirements.
     * 
     */
    public YFastTrieLong() {
        
        this.w = 62;
        
        maxC = (1L << w) - 1;
        
        binSzModel = BinSizeModel.FAST_RT;
                
        binSz = w;
        
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
    
    /**
     * 
     * @param numberOfEntries
     * @param binSizeModel 
       <pre>
        3 choices:
         AUTOMATIC = code looks at memory available to heap and if there is 
             enough, chooses the fast runtime model of operations O(log_2(w)).
              This is the default.
         FAST_RT = user has specified choice should be the fastest runtime 
              model.  O(log_2(w))
         SPACE_CONSERTVING = user has specified choice should be the smallest 
              memory use model w/ runtimes that are O(log_2(N_entries/w)) which 
              ends up being approx 
              O(20) for nEntries=35e6 and O(12) for nEntries=35e3.
        </pre>
     */
    public YFastTrieLong(int numberOfEntries, BinSizeModel binSizeModel) {
        
        this.w = 62;
        
        maxC = (1L << w) - 1;
        
        this.binSzModel = binSizeModel;
                
        binSz = calculateBinSize(numberOfEntries, w, binSizeModel);
        
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
     * 
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
     * 
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

    private static long calculateBinSize(int numberOfEntries, 
        int maxNumberOfBits, BinSizeModel theBinSizeModel) {

        /*        
        avail = max memory to use
              = (xFt+xftReps holding at most nBins inserts each) 
                + nBins * (nRBs trees) + nEntries
        
             largest scaling component is nEntries nodes
        
        wanting small nBins, but with nEntries/nBins not too large nor too small
        
        let w = max bit length specified by user, else is default 62 bits.
        
        a few examples here with data in mind for image pixel indexes, in which
        the maximum pixel index = width * height:
        
        if nBins = w:
            nEntries        w  nBins   binsz,   runtime      space
           
            5000*7000=35e6  62  62     564516     19       in progress
            5000*7000=35e6  26  26   1346153      20 
         .1*5000*7000=35e5  62  62     56451      16
         .1*5000*7000=35e5  22  22     56451      18
        
            500*700=35e4    62  62      5645      13
            500*700=35e4    19  19     18421      14
         .1*500*700=35e3    62  62       564      10
         .1*500*700=35e3    16  16       564      12
        
        can see that can improve the runtime by decreasing the binSz at
        expense of increasing nBins (which leads to creating more items
        on the space hog xfasttrielong).
        
        If instead, one uses binSz = w to get the suggested yfasttrie runtime,
        one would have nEntries/binSz number of representatives.
        This becomes a larger consumer of resources, though still 
        better than an xfasttrie alone by a factor of w.
        
        The calculations above for the suggested yfasttrie model:
        
        if binsz = w:
            nEntries        w  nBins   binsz,   runtime      space
           
            5000*7000=35e6  62  564516  62       6       in progress
            5000*7000=35e6  26  1346153 26       5   
         .1*5000*7000=35e5  62  56451   62       6
         .1*5000*7000=35e5  22  56451   22       5
        
            500*700=35e4    62  5645    62       6   
            500*700=35e4    19  18421   19       5    
         .1*500*700=35e3    62  564     62       6    
         .1*500*700=35e3    16  564     16       4    

        Recalculating the memory requirements for both models
        and if have enough memory, prefering the one with 
        better runtime unless a user setting requests smaller space use.        
        */
        
        if (theBinSizeModel.equals(BinSizeModel.FAST_RT)) {
            
            // binSize = w bits
            return maxNumberOfBits;
            
        } else if (theBinSizeModel.equals(BinSizeModel.SPACE_CONSERTVING)) {
            
            // nBins = w bits
            // binSize = number of entries / w
            
            return numberOfEntries/maxNumberOfBits;
        }
        
        long maxNumber = (1L << maxNumberOfBits) - 1;

        long totalMemory = Runtime.getRuntime().totalMemory();
        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        long heapUsage = mbean.getHeapMemoryUsage().getUsed();
        long avail = totalMemory - heapUsage;

        //double rt = Math.log(binSz)/Math.log(2);
        
        long[] est = estimateSizeOnHeap(numberOfEntries, 
            maxNumberOfBits, BinSizeModel.FAST_RT);
        
        System.out.println("Memory size estimates for fastest: " + 
            Arrays.toString(est) + " avail=" + avail);

        long fastMemory = est[0];
        
        if (fastMemory < avail) {
            return maxNumberOfBits;
        } else if (est[1] < avail) {
            long remaining = avail - est[1];
            if (((double)remaining/(double)est[1]) > 0.25) {
                return maxNumberOfBits;
            }
        }
        
        return numberOfEntries/maxNumberOfBits;
    }
    
    protected long getBinSz() {
        return binSz;
    }
    
    /**
     * estimate the size that an instance of YFastTrieLong with
     * n added entries, maxNumberOfBits, and 
     * use binSzModel
     * would occupy in heap space in Bytes.
     * 
     * NOTE: there are some varying components to estimating the memory that
     * depend upon the properties of the numbers inserted.
     * For example:
     * <pre>
     *     -- xft is an XFastTrie instantiated with maxNumberOfBits.
     *        It will have at most, 
     *        * nBins number of entries, where nBins
     *        is determined by the BinSizeModel.
     *        In addition to the number of inserted items (which is only one
     *        per bin of numberOfEntries), there will be some undetermined
     *        number of prefix nodes created in the process.
     *        A factor of 5 more is assumed here to over estimate the total 
     *        number of trie nodes that includes the internal prefix nodes.
     *        -- THE LOGIC is still in progress to determine
     *           an upper and lower limit to estimate the number of populated 
     *           nBins w/o knowing properties of the numbers, such as whether 
     *           they are sequential, or have large gaps, etc.
     *    -- xftReps is a hashMap with same number of inserts as xft,
     *       so has the same need for an upper and lower estimate.
     * </pre>
     * 
     * @param numberOfEntries amount of space for this object's instance
     * with n entries in Bytes on the heap.
     * @param maxNumberOfBits all entries must have bit lengths .lte. this
     * @param binSzModel must not be AUTOMATIC
     * 
     * @return array with 2 estimates, (1) estimate using all bins and a
     * factor of 5 for creating trie prefix nodes,
       (2) estimate from using 1/4 of the bins and a factor of 3 for creating
       the trie prefix nodes.
     */
    public static long[] estimateSizeOnHeap(int numberOfEntries, int
        maxNumberOfBits, BinSizeModel binSzModel) {
        
        if (binSzModel.equals(BinSizeModel.AUTOMATIC)) {
            throw new IllegalArgumentException("binSizeModel cannot be "
                + " AUTOMATIC.  That decision is based upon sing this method");
        }
    
        long ww = maxNumberOfBits;
        
        long maxNumber = (1L << ww) - 1;
        
        long binSz;
        if (binSzModel.equals(BinSizeModel.FAST_RT)) {
            binSz = maxNumberOfBits;
        } else {
            binSz = numberOfEntries/maxNumberOfBits;
        }
        
        int nBins = (int)Math.ceil((double)maxNumber/(double)binSz);
        
        long total = 0;
        
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNIntFields(2);
        est.setNLongFields(3);
        est.setNBooleanFields(1);
        //objects: xft, xftReps, rbs
        est.setNObjRefsFields(3);
       
        total += est.estimateSizeOnHeap();
       
        // --------- include contents of the objects --------------
        
        /*
         the minimum of each bin range, if it is populated, is the representative
         node, and that is held in 2 data structures:
            xft holds the number, allowing fast repr prev and next lookups.
            xftReps holds the repr as the value, found by key = binNumber.
         * the max number of trie entries will be nBins
             but there will be prefix trie nodes too
        private final XFastTrieLong<XFastTrieNodeLong<Long>, Long> xft;
        
        // key = bin index (which is node/binSz), value = repr value.
        // each repr value is the minimum stored in the bin.
        // * the max number of map entries will be nBins.
        private final TLongLongMap xftReps = new TLongLongHashMap();
   
        // all inserts of this class are held in 
        //    * at most nBins number of trees which each 
        //      hold at most binSz number of entries.
        // each list item is a sorted binary search tree of numbers in that bin.
        //    the value in the tree holds multiplicity of the number.
        // each list index can be found by node/binSz
        // each sorted tree has 
        //    key = node (inserted number), w/ value=
        //        the number of times that number is present (multiplicity).
        TLongObjectMap<RedBlackBSTLongInt> rbs;
        */
    
        // returning 2 estimates
        // (1) estimate using all bins w/ factor 5 for tries
        // (2) estimate from using 1/4 of the bins w/ factor 3 for tries
                
        // using factor of 5 for total w/ prefix nodes
        long total2_1 = numberOfEntries * 5 *
            XFastTrieNodeLong.estimateSizeOnHeap();
        
        // all nBins are filled w/ a repr
        total2_1 += XFastTrieLong.estimateSizeOnHeap(nBins, maxNumberOfBits);
        
        long total2_2 = numberOfEntries * 3 *
            XFastTrieNodeLong.estimateSizeOnHeap();
        
        // 1/4 of nBins are filled w/ a repr
        total2_2 += XFastTrieLong.estimateSizeOnHeap(nBins/4, maxNumberOfBits);
        
        //TLongLongMap
        total2_1 += ObjectSpaceEstimator.estimateTLongLongHashMap();
        
        //nBins number of repr entries in map
        total2_1 += (2 * nBins * ObjectSpaceEstimator.estimateLongSize());
        
        
        //TLongLongMap
        total2_2 += ObjectSpaceEstimator.estimateTLongLongHashMap();
        
        //nBins/4 number of repr entries in map
        total2_1 += ((nBins/2) * ObjectSpaceEstimator.estimateLongSize());
        
        // 1 TLongObjectMap<RedBlackBSTLongInt> rbs;
        total2_1 += ObjectSpaceEstimator.estimateTLongObjectHashMap();
        
        total2_2 += ObjectSpaceEstimator.estimateTLongObjectHashMap();
        
        // nBins number of RedBlackBSTLongInt
        long rbtree = RedBlackBSTLongInt.estimateSizeOnHeap(0);
        
        total2_1 += (nBins * rbtree);
        
        total2_2 += ((nBins/4) * rbtree);
        
        // nEntries number of long, int nodes
        total2_1 += numberOfEntries * RedBlackBSTLongInt.estimateNodeSizeOnHeap();
        
        total2_2 += numberOfEntries * RedBlackBSTLongInt.estimateNodeSizeOnHeap();
            
        return new long[]{total2_1, total2_2};
    }
    
    public void debugPrint() {
        
        long[] binIndexes = rbs.keys();
        Arrays.sort(binIndexes);
        for (long binIdx : binIndexes) {
            System.out.println("binNumber=" + binIdx);
            RedBlackBSTLongInt rbt = rbs.get(binIdx);
            rbt.printPreOrderTraversal();
        }
    }
}
