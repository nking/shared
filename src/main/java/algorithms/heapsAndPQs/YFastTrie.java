package algorithms.heapsAndPQs;

import algorithms.util.ObjectSpaceEstimator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import thirdparty.ods.Integerizer;
import thirdparty.ods.XFastTrie;
import thirdparty.ods.XFastTrieNode;

/** 
 * 
 * from wikipedia
 *     https://en.wikipedia.org/wiki/Y-fast_trie
 *  
 * a y-fast trie is a data structure for storing 
 * integers from a bounded domain. It supports exact and predecessor 
 * or successor queries in time O(log log M), using O(n) space, 
 * where n is the number of stored values and M is the number of bits
 * of the maximum value in the domain. 
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
public class YFastTrie {

    /*    
    designing from browsing a few different lecture notes
    online. the yfast trie uses same w and maxC as
    the XFastTrie where maxC is the maximum value that the
    trie will hold and w is the number of bits needed to
    represent maxC.
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
    */
    
    protected int n = 0;
    
    private final int w;
    
    private final int maxC;
    
    private final int binSz;
    
    private int nBins;

    private final XFastTrie<XFastTrieNode<Integer>, Integer> xft;
    
    // key = bin index, value = repr value.
    // each repr value is the minimum stored in the bin.
    private final TIntIntMap xftReps = new TIntIntHashMap();
    
    // there are w items in rbs
    // each list item is a sorted binary search tree of numbers in that bin.
    //    the value is the tree holds the number of times that number is present.
    private final TIntObjectMap<TreeMap<Integer, Integer>> rbs;

    public YFastTrie(int wBits) {
        
        if (wBits < 31 && wBits > 1) {
            this.w = wBits;
        } else {
            throw new IllegalStateException("wBits "
                + " should be greater than 1 and less than 32"
                + " wBits=" + wBits);
        }
        maxC = (1 << w) - 1;

        // for the runtime complexity of operations to be 
        // math.log( math.log(maxC)/math.log(2))/math.log(2)
        // each balanced binary search tree is size
        //   n = binSz = math.log(maxC)/math.log(2)
        //
        // therefore, splitting maxC numbers so that each bin size is n
        // requires nBins = maxC / n

        // e.g. given w = 31 bits, then maxC = 2.15E9, n = 31(==w), nBins=7E7
        // e.g. given w = 16 bits, then maxC = 65535, n = 16(==w), nBins=4.1E3

        binSz = w;

        nBins = (int)Math.ceil((double)maxC/(double)binSz);
        
        //System.out.println("nBins=" + nBins + "  rt of ops=" +
        //    (Math.log(binSz)/Math.log(2)));
        
        rbs = new TIntObjectHashMap<TreeMap<Integer, Integer>>();
         
        XFastTrieNode<Integer> clsNode = new XFastTrieNode<Integer>();
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        // xft operations have runtime complexity O(log_2(w)) 
        //      w/ large space complexity of O(n * w), hence reducing n here: 
        xft = new XFastTrie<XFastTrieNode<Integer>, Integer>(clsNode, it, w);
    }
    
    /** construct with default maximum size of numbers to store being Integer.MAX_VALUE.
    The operations for this instance will have runtime complexity 
    O(log_2(31)) = O(5).
    */
    public YFastTrie() {
        
        this.w = 30;
        
        maxC = (1 << w) - 1;
                                
        binSz = w;
        
        nBins = (int)Math.ceil((double)maxC/(double)binSz);
                
        //System.out.println("nBins=" + nBins + "  rt of ops=" +
        //    (Math.log(binSz)/Math.log(2)));
        
        rbs = new TIntObjectHashMap<TreeMap<Integer, Integer>>();
        
        XFastTrieNode<Integer> clsNode = new XFastTrieNode<Integer>();
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        xft = new XFastTrie<XFastTrieNode<Integer>, Integer>(clsNode, it, w);
    }
    
    // runtime complexity is O(1)
    protected TreeMap<Integer, Integer> getTreeMap(int index) {
        Integer key = Integer.valueOf(index);
        TreeMap<Integer, Integer> map = rbs.get(key);
        if (map == null) {
            map = new TreeMap<Integer, Integer>();
            rbs.put(key, map);
        }
        return map;
    }

    /**
     * runtime complexity is O(log_2(wBits))
     * @param node
     * @param index 
     */
    private void addToRBTree(int node, int index) {
        
        // runtime complexity is O(1)
        TreeMap<Integer, Integer> map = getTreeMap(index);
        
        assert(map != null);
        
        Integer key = Integer.valueOf(node);
        
        // runtime complexity is O(log_2(w))
        Integer multiplicity = map.get(key);
    
        if (multiplicity == null) {
            multiplicity = Integer.valueOf(1);
        } else {
            multiplicity = Integer.valueOf(1 + multiplicity.intValue());
        }
        
        map.put(key, multiplicity);        
    }
    
    /**
     * runtime complexity is O(log_2(wBits))
     * @param node
     * @param index 
     */
    private boolean deleteFromRBTree(int node, int index) {
                
        TreeMap<Integer, Integer> map = getTreeMap(index);
        
        assert(map != null);
        
        Integer key = Integer.valueOf(node);

        Integer multiplicity = map.get(key);
    
        if (multiplicity == null) {
            return false;
        }
        
        if (multiplicity.intValue() > 0) {
            multiplicity = Integer.valueOf(multiplicity.intValue() - 1);
            if (multiplicity.intValue() > 0) {
                map.put(key, multiplicity);
            }
        }
        if (multiplicity.intValue() == 0) {
            map.remove(key);
        }
        
        return true;
    }
    
    /**
     * runtime complexity is roughly O(log_2(wBits)).
     * @param node a number >= 0 and having bit length 
     * less than or equal to w.
     * @return 
     */
    public boolean add(int node) {

        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node.key must "
                + "be less than " + maxC + " node=" + node);
        }
        
        int index = node/binSz;
        
        //O(1) to find the repr
        int existingRepr = xftReps.get(index);
                
        // worse case runtime here is O(log_2(w)) + O(l-w), else is 0
        if (!xftReps.containsKey(index)) {
            // insert is O(log_2(w)) + O(l-w)
            xft.add(Integer.valueOf(node));
            xftReps.put(index, node);
        } else if (node < existingRepr) {
            // delete is O(log_2(w)) + O(l-w)
            // insert is O(log_2(w)) + O(l-w)
            xft.remove(Integer.valueOf(existingRepr));
            xft.add(Integer.valueOf(node));
            xftReps.put(index, node);
        }
                
        // runtime complexity here is
        addToRBTree(node, index);
        
        n++;
        
        return true;
    }

    /**
     * runtime complexity is roughly O(log_2(wBits)).
     * 
     * @param node
     * @return 
     */
    public boolean remove(int node) {
        
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC);
        }
        
        int index = node/binSz;
                
        boolean removed = deleteFromRBTree(node, index);
                
        if (!removed) {
            return false;
        }
        
        if (!xftReps.containsKey(index)) {
            return false;
        }
        
        TreeMap<Integer, Integer> map = getTreeMap(index);
      
        int existingRepr = xftReps.get(index);
      
        if (map.isEmpty()) {
            // just deleted the last item so remove from rbs
            // delete is O(log_2(w)) + O(w-l)
            if (xftReps.containsKey(index)) {
                xft.remove(Integer.valueOf(existingRepr));
                xftReps.remove(index);
            }
        } else if (node == existingRepr) {
            
            //existingRepr is maintained as the minimum in the bin,
            //   so if a node w/ this value is removed and the multiplicity
            //      was 1, need to assign a new repr
            
            // O(log_2(N/w))
            Integer multiplicity = map.get(Integer.valueOf(node));
    
            if (multiplicity == null) {
                // remove the current repr and assign a new one
                // delete is O(log_2(w)) + O(w-l)
                xft.remove(Integer.valueOf(existingRepr));
                xftReps.remove(index);
            
                // O(log_2(N/w))
                Entry<Integer, Integer> minEntry = map.firstEntry(); 
                xft.add(minEntry.getKey());
                xftReps.put(index, minEntry.getKey()); 
            }            
        }
        
        n--;
        
        return true;
    }

    /**
     * runtime complexity is roughly O(log_2(wBits)).
     * 
     * @param node
     * @return returns node if found, else -1
     */
    public int find(int node) {
                
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC + ". node=" + node);
        }
                
        int index = node/binSz;
                
        TreeMap<Integer, Integer> map = getTreeMap(index);
        
        Integer multiplicity = map.get(Integer.valueOf(node));
        if (multiplicity == null) {
            return -1;
        }
        
        return node;
    }

    /**
     * runtime complexity is roughly O(log_2(wBits)).
     * 
     * @param node
     * @return value preceding node, else -1 if there is not one
     */
    public int predecessor(int node) {
    
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC);
        }
        
        Integer nodeKey = Integer.valueOf(node);
        
        int nodeIndex = node/binSz;
        
        // the repr is stored in xft and it is always the minium for the bin
        boolean isAMinimum = xft.find(nodeKey) != null;
        
        /*
        if the node is not a minima, the answer is in
           the node's map if its size is larger > 1
        */
        
        TreeMap<Integer, Integer> map = getTreeMap(nodeIndex);
        
        if (!isAMinimum && (map.size() > 1)) {
            Entry<Integer, Integer> pred = map.lowerEntry(nodeKey);
            if (pred != null) {
                return pred.getKey().intValue();
            }
        }
       
        // else, predeccessor is in the closest bin < nodeIndex that has
        //    items in it.
                
        Integer prev = xft.predecessor(nodeKey);
        if (prev == null) {
            return -1;
        }
        
        int prev0Index = prev.intValue()/binSz;
            
        map = getTreeMap(prev0Index);
                
        Entry<Integer, Integer> lastItem = map.lastEntry();
               
        if (lastItem == null) {
            return -1;
        }
        
        return lastItem.getKey();
    }
    
    /**
     * runtime complexity is roughly O(log_2(wBits)).
     * 
     * @param node
     * @return 
     */
    public int successor(int node) {
                
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC);
        }
        
        Integer nodeKey = Integer.valueOf(node);
        
        int nodeIndex = node/binSz;
        
        boolean isAMinimum = xft.find(nodeKey) != null;
        
        TreeMap<Integer, Integer> nodeMap = getTreeMap(nodeIndex);
        
        if (isAMinimum) {
            // if tree size > 1, the next key is the successor
            // else, the xft sucessor to nodeIndex is the successor
            
            if (nodeMap.size() > 1) {
                Entry<Integer, Integer> successor = nodeMap.higherEntry(nodeKey);
                assert(successor != null);
                return successor.getKey();
            }
            
            Integer successorRepr = xft.successor(nodeKey);
            if (successorRepr == null) {
                return -1;
            }
            
            // the successor representative is then the next value
            return successorRepr;
        }
        
        // else, the node is not a repr
        //   if there is a tree successor to the node, that is the successor
        //   else, the xft successor to nodeIndex is the successor
                    
        Entry<Integer, Integer> sEntry = nodeMap.higherEntry(nodeKey);
        
        if (sEntry != null) {
            return sEntry.getKey();
        }
        
        Integer successorRepr = xft.successor(nodeKey);
        if (successorRepr == null) {
            return -1;
        }

        // the successor representative is then the next value
        return successorRepr;
    }

    /**
     * runtime complexity is roughly O(log_2(wBits)).
     * @return minimum, else -1 if empty
     */
    public int minimum() {
        
        if (xft.size() == 0) {
            return -1;
        }
        
        Integer repr = xft.minimum();
        
        assert(repr != null);
       
        return repr.intValue();
    }

    /**
     * runtime complexity is roughly O(log_2(wBits)).
     * 
     * @return maximum, else -1 if empty
     */
    public int maximum() {
        
        if (xft.size() == 0) {
            return -1;
        }
        
        Integer maxRepr = xft.maximum();
        
        assert(maxRepr != null);
        
        int index = maxRepr.intValue()/binSz;
        
        TreeMap<Integer, Integer> map = getTreeMap(index);
        
        assert(map != null);
        
        Entry<Integer, Integer> lastItem = map.lastEntry();
        
        assert(lastItem != null);
        
        return lastItem.getKey();
    }
    
    /**
     * runtime complexity is roughly O(log_2(wBits)).
     * 
     * @return minumum, else -1 if empty
     */
    public int extractMinimum() {
        
        //O(log_2(w))
        int min = minimum();

        if (min == -1) {
            assert(xft.size() == 0);
            return -1;
        }
                
        remove(min);
        
        return min;
    }
    
    /**
     * runtime complexity is roughly O(log_2(wBits)).
     * @return maximum, else -1 if empty
     */
    public int extractMaximum() {
        
        int max = maximum();

        if (max == -1) {
            assert(xft.size() == 0);
            return -1;
        }
                
        remove(max);
        
        return max;
    }
    
    /**
    get the number of items stored in the trie.  runtime complexity is O(1).
    @return number of items in the trie.
    */
    public int size() {
        return n;
    }

    protected int getBinSz() {
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
        TLongObjectMap<RedBlackBSTLongInt2> rbs;
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
            XFastTrieNode.estimateSizeOnHeap();
        
        // all nBins are filled w/ a repr
        total2_1 += XFastTrie.estimateSizeOnHeap(numberOfEntries);
        
        long total2_2 = numberOfEntries * 3 *
            XFastTrieNode.estimateSizeOnHeap();
        
        // nBinsSparse of nBins are filled w/ a repr
        total2_2 += XFastTrie.estimateSizeOnHeap(numberOfEntries);
        
        
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
        
        // nBins number of TreeMap<Integer, Integer> 
        ObjectSpaceEstimator est2 = new ObjectSpaceEstimator();
        est2.setNBooleanFields(1);
        est2.setNObjRefsFields(2);
        long totalEntry = est2.estimateSizeOnHeap();
        totalEntry += 3. * totalEntry;
        est2 = new ObjectSpaceEstimator();
        est2.setNIntFields(2);
        long rbtree = est2.estimateSizeOnHeap() + totalEntry;
        long rbtreeNodes = numberOfEntries * totalEntry;
        
        total2_1 += (nBins * rbtree);
        
        total2_2 += (nBinsSparse * rbtree);
        
        
        // nEntries number of long, int nodes
        
        total2_1 += rbtreeNodes;
        
        total2_2 += rbtreeNodes;
           
        return new long[]{total2_1 + total, total2_2 + total};
    }
    
}
