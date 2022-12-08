package algorithms.bipartite;

import algorithms.heapsAndPQs.Heap;
import algorithms.heapsAndPQs.HeapNode;

/**
 *  NOTE: each heap has O(1) decreaseKey, but if multiple bins are used,
 *     the method takes O(log_2(n)) + O(1).
 * 
 * @author nichole
 */
public class FibonacciHeapWrapper {

    private final int nApprox;
    
    private final Heap[] heaps;
    
    private long lastKnownMinKey = 0;
    private int lastKnownMinKeyIdx = 0;
    private long lastKnownMaxKey = -1;
    
    private final int binSz;
    
    private int n = 0;
    
    /**
     *
     @param nEstimate
     @param maxC
     */
    public FibonacciHeapWrapper(int nEstimate, int maxC) {
        
        int pow2 = (int)Math.ceil(Math.log(nEstimate)/Math.log(2));
    
        nApprox = nEstimate;
        
        if (pow2 > 6 && (maxC > 1)) {
            
            // make multiple maps separated by sequential partitions of values
            int nBins = pow2;
            
            if (maxC < pow2) {
                //pow2 range is 5 to 31
                nBins = maxC;
            }
            
            binSz = (int)Math.ceil((float)maxC/(float)nBins);
            
            heaps = new Heap[nBins];
            
        } else {
            heaps = new Heap[1];
            binSz = Integer.MAX_VALUE;
        }
    }
    
    /**
     * runtime is O(1).  makes no attempt to consolidate tree.
     @param node 
     */
    public void insert(HeapNode node) {
        
        int key = (int)node.getKey();
        
        int binIdx;
        if (heaps.length == 1) {
            binIdx = 0;
        } else {
            binIdx = key/binSz;
        }
        
        if (heaps[binIdx] == null) {
            heaps[binIdx] = new Heap();
        }
        
        heaps[binIdx].insert(node);
        
        n++;
        
        if (key < lastKnownMinKey) {
            lastKnownMinKey = key;
            lastKnownMinKeyIdx = binIdx;
        }
        if (key > lastKnownMaxKey) {
            lastKnownMaxKey = key;
        }
    }

    /**
     * note, key1 and key2 are in same heap, the runtime is 
     * O(1) for decrease key, else the runtime is a 
     * delete and insert as O(log_2(n)).
     * 
     @param node
     @param key2 
     */
    public void decreaseKey(HeapNode node, long key2) {
        
        int key = (int)node.getKey();
        
        int binIdx1, binIdx2;
        if (heaps.length == 1) {
            binIdx1 = 0;
            binIdx2 = 0;
        } else {
            binIdx1 = key/binSz;
            binIdx2 = (int)key2/binSz;
        }
                
        if (binIdx1 == binIdx2) {
            
            heaps[binIdx1].decreaseKey(node, key2);
            
        } else {
        
            //runtime is that of extractMin, O(log_2(n))
            heaps[binIdx1].remove(node);
                
            node.setKey(key2);
     
            if (heaps[binIdx2] == null) {
               heaps[binIdx2] = new Heap();
            }
            
            //runtime is O(1).  makes no attempt to consolidate tree.
            heaps[binIdx2].insert(node);
        }
        
        if (key2 < lastKnownMinKey) {
            lastKnownMinKey = key2;
            lastKnownMinKeyIdx = binIdx2;
        }
    }
    
    /**
     * runtime is O(log_2 N) or better.
     @return 
     */
    public HeapNode extractMin() {
        
        if (n == 0) {
            return null;
        }
     
        for (int i = (int)lastKnownMinKey; i < heaps.length; ++i) {
            
            if ((heaps[i] == null) || (heaps[i].getNumberOfNodes() == 0)) {
                continue;
            }
            
            HeapNode node = heaps[i].extractMin();
            lastKnownMinKey = i;
            n--;
            
            return node;
        }
        
        return null;
    }
    
    /**
     * runtime complexity is O(log_2(n))
     @param node 
     */
    public void remove(HeapNode node) {
        
        int key = (int)node.getKey();
        
        int binIdx;
        if (heaps.length == 1) {
            binIdx = 0;
        } else {
            binIdx = key/binSz;
        }
        //runtime is that of extractMin, O(log_2(n))
        heaps[binIdx].remove(node);
        
        n--;
    }
    
    /**
     *
     @return
     */
    public long getNumberOfNodes() {
        return n;
    }
}
