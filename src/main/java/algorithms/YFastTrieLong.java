package algorithms;

import algorithms.util.ObjectSpaceEstimator;
import gnu.trove.map.TLongLongMap;
import gnu.trove.map.TLongObjectMap;
import gnu.trove.map.hash.TLongLongHashMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import java.util.Arrays;
import thirdparty.edu.princeton.cs.algs4.RedBlackBSTLongInt;
import thirdparty.ods.Longizer;
import thirdparty.ods.XFastTrieLong;
import thirdparty.ods.XFastTrieNodeLong;

/**
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
 
   NOTE that one may need to use the object heap size estimator before using
   this on a large number of objects and compare the result to the
   available heap memory.
        long totalMemory = Runtime.getRuntime().totalMemory();
        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        long heapUsage = mbean.getHeapMemoryUsage().getUsed();
        long avail = totalMemory - heapUsage;

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
      
YFastTrie
   
    let w = max bit length specified by user, else is default 62 bits.

    a few examples here with data in mind for image pixel indexes, in which
    the maximum pixel index = width * height:

    If instead, one uses binSz = w to get the suggested yfasttrie runtime,
    one would have nEntries/binSz number of representatives.
    This is fewer entries into xfasttrie than if using xfasttrie alone 
    so is conserving space.

    The calculations above for the suggested yfasttrie model:

    binsz = w

        nEntries        w  nBins   binsz,   runtime    space(MB)

        5000*7000=35e6  62  564516  62       6       505405, 61507
        5000*7000=35e6  26  1346153 26       5       14515, 9141
     .1*5000*7000=35e5  62  56451   62       6       492908, 53336
     .1*5000*7000=35e5  22  56451   22       5       1435, 912

        500*700=35e4    62  5645    62       6       491658, 52519
        500*700=35e4    19  18421   19       5       145, 91
     .1*500*700=35e3    62  564     62       6       491533, 52437
     .1*500*700=35e3    16  564     16       4       14, 9
   
    */
    
    private int n = 0;
    
    private final int w;
    
    private final long maxC;
    
    private final long binSz;
    
    private long nBins;
    
    /**
       the minimum of each bin range, if it is populated, is the representative
         node, and that is held in 2 data structures:
            xft holds the number, allowing fast repr prev and next lookups.
            xftReps holds the repr as the value, found by key = binNumber.
         * the max number of trie entries will be nBins
             but there will be prefix trie nodes too
     */
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
     * constructor using the the maximum number of bits of 62
     * and the default model for the fast runtime
     * which may be expensive in space requirements.
     * 
     */
    public YFastTrieLong() {
        
        this.w = 62;
        
        maxC = (1L << w) - 1;
                        
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
     * add node to the data structure.
     * 
     * runtime complexity should usually be O(log_2(w)) but
     * filling the prefix tree in the xfasttrie can sometimes add a small amount
     * + O(l-w)
     * where w is the maximum bit length and l is a level in the prefix tree.
     * 
     * @param node a number >= 0 and having bit length 
     * less than or equal to w.
     * @return true if successfully added node.
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
     * remove node from the data structure.
     * 
     * runtime complexity should usually be O(log_2(w)) but
     * filling the prefix tree in the xfasttrie can sometimes add a small amount
     * + O(l-w)
     * where w is the maximum bit length and l is a level in the prefix tree.
     * 
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
     * find node in the data structure and return it, else return -1.
     * 
     * runtime complexity is at most O(log_2(w)) but since the map might not
     * be completely populated, the complexity might be smaller.
     * 
     * @param node
     * @return the value of the node if present, else -1
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
     * find the largest node smaller in value than node in the datastructure.
     * 
     * runtime complexity should usually be O(log_2(w)) but
     * filling the prefix tree in the xfasttrie can sometimes add a small amount
     * + O(l-w)
     * where w is the maximum bit length and l is a level in the prefix tree.
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
     * find the smallest value larger than node in the data structure.
     * 
     * runtime complexity should usually be O(log_2(w)) but
     * filling the prefix tree in the xfasttrie can sometimes add a small amount
     * + O(l-w)
     * where w is the maximum bit length and l is a level in the prefix tree.
     * 
     * @param node
     * @return the next node in the ordered data strucure, else -1 if no such
     * node exists.
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
        
        if (tree.size() > 1) {
            
            // if tree size > 1, the next key is the successor
            // else, the xft sucessor to nodeIndex is the successor

            long[] output = new long[2];
            tree.higher(node, output);
                
            if (isAMinimum) {
                assert(output[0] != -1);
                return output[1];
            }
        
            // else, the node is not a repr
            //   if there is a tree successor to the node, that is the successor
            //   else, the xft successor to nodeIndex is the successor
        
            if (output[0] != -1) {
                return output[1];
            }
        }
                
        Long successorRepr = xft.successor(nodeKey);
        if (successorRepr == null) {
            return -1;
        }

        // the successor representative is then the next value
        return successorRepr;
    }

    /**
     * find the smallest node in the datastructure.
     * 
     * runtime complexity should usually be O(log_2(w)) but
     * filling the prefix tree in the xfasttrie can sometimes add a small amount
     * + O(l-w)
     * where w is the maximum bit length and l is a level in the prefix tree.
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
     * find the largest node in the data structure.
     * runtime complexity should usually be O(log_2(w)) but
     * filling the prefix tree in the xfasttrie can sometimes add a small amount
     * + O(l-w)
     * where w is the maximum bit length and l is a level in the prefix tree.
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
     * find and remove the smallest value in the data structure.
     * 
     * runtime complexity should usually be O(log_2(w)) but
     * filling the prefix tree in the xfasttrie can sometimes add a small amount
     * + O(l-w)
     * where w is the maximum bit length and l is a level in the prefix tree.
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
     * find and remove the largest value in the data structure.
     * runtime complexity should usually be O(log_2(w)) but
     * filling the prefix tree in the xfasttrie can sometimes add a small amount
     * + O(l-w)
     * where w is the maximum bit length and l is a level in the prefix tree.
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
    
    /**
     * get the current number of values stored in the data structure.
     * @return 
     */
    public int size() {
        return n;
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
     * 
     * @return array with 2 estimates, (1) estimate using all bins and a
     * factor of 5 for creating trie prefix nodes,
       (2) estimate from using 1/4 of the bins and a factor of 3 for creating
       the trie prefix nodes.
     */
    public static long[] estimateSizeOnHeap(int numberOfEntries, int
        maxNumberOfBits) {
        
        long ww = maxNumberOfBits;
        
        long maxNumber = (1L << ww) - 1;
        
        long binSz = maxNumberOfBits;
        
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
        // (2) estimate from using nBinsSparse of the nBins w/ factor 3 for tries
        
        int nBinsSparse = nBins/10;
        if (nBinsSparse < 1) {
            nBinsSparse = 1;
        }
        
        // using factor of 5 for total w/ prefix nodes
        long total2_1 = numberOfEntries * 5 *
            XFastTrieNodeLong.estimateSizeOnHeap();
        
        // all nBins are filled w/ a repr
        total2_1 += XFastTrieLong.estimateSizeOnHeap(nBins, maxNumberOfBits);
        
        long total2_2 = numberOfEntries * 3 *
            XFastTrieNodeLong.estimateSizeOnHeap();
        
        // nBinsSparse of nBins are filled w/ a repr
        total2_2 += XFastTrieLong.estimateSizeOnHeap(nBinsSparse, maxNumberOfBits);
        
        
        //TLongLongMap
        total2_1 += ObjectSpaceEstimator.estimateTLongLongHashMap();
        
        //nBins number of repr entries in map
        total2_1 += (2 * nBins * ObjectSpaceEstimator.estimateLongSize());
        
        
        //TLongLongMap
        total2_2 += ObjectSpaceEstimator.estimateTLongLongHashMap();
        
        //nBins/4 number of repr entries in map
        total2_2 += (2 * nBinsSparse * ObjectSpaceEstimator.estimateLongSize());
        
        
        // 1 TLongObjectMap<RedBlackBSTLongInt> rbs;
        total2_1 += ObjectSpaceEstimator.estimateTLongObjectHashMap();
        
        total2_2 += ObjectSpaceEstimator.estimateTLongObjectHashMap();
        
        
        // nBins number of RedBlackBSTLongInt
        long rbtree = RedBlackBSTLongInt.estimateSizeOnHeap(0);
        
        total2_1 += (nBins * rbtree);
        
        total2_2 += (nBinsSparse * rbtree);
        
        
        // nEntries number of long, int nodes
        total2_1 += numberOfEntries * RedBlackBSTLongInt.estimateNodeSizeOnHeap();
        
        total2_2 += numberOfEntries * RedBlackBSTLongInt.estimateNodeSizeOnHeap();
            
        return new long[]{total2_1 + total, total2_2 + total};
    }
    
    /**
     * print all entries in data structure to standard out.
     */
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
