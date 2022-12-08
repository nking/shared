package algorithms.bipartite;

import algorithms.heapsAndPQs.YFastTrie;
import algorithms.heapsAndPQs.HeapNode;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.util.logging.Logger;

/**
 * A min heap for the MinCostUnbalancedAssignment.
 * It uses the "YFastTrie min priority queue algorithm" 
 * pattern by default if the VM has enough memory, else
 * uses a Fibonacci Heap.
 * All operations for the "YFastTrie" are constant time.

    <pre>
     prefer to use YFastTrie as all operations are essentially
     O(log_2 log_2(w)) where w is the largest integer to store in trie.
     in contrast Min-Heap/PriorityQueue and Fibonacci have these trade-offs
                                 FibHeap             MinHeap/PQ
      extractMin()        O(log_2(n))                  O(1)
      insert()                   O(1)                  O(log_2(n))
      decreaseKey()       O(log_2(n))                  O(log_2(n))
    </pre>
 * 
 * The YFastTrie has O(log log(M)) operations including successor and
 * predecessor where M is the maximum value in the domain.
 * 
 * @author nichole
 */
public class MinHeapForRT2012 {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    // 1 = Fibonacci, 2 = YFastTrie
    private final int algorithm;
        
    private int lastKnownMinKey0 = 0;
    
    // for use in tuning the capacity
    private int lastKnownMaxKey0 = 0;
    
    private final FibonacciHeapWrapper heap1;

    private final YFastTrieWrapper heap2;
    
    void printLastKnownMinMax() {
        log.fine("min=" + lastKnownMinKey0
            + " max=" + lastKnownMaxKey0);
    }
    
    /**
     * 
     @param maxValue estimate of maximum value to store.
     @param approxN approximate number of nodes expected
     * to be in the heap as a rough maximum at a given time.
     * (it's used to help determine which algorithm to use
     * internally).
     * 
     * IA min heap for the MinCostUnbalancedAssignment.
     * It uses the "YFastTrie min priority queue algorithm" 
     * pattern by default if the VM has enough memory, else
     * uses a Fibonacci Heap.
     * All operations for the "YFastTrie" are constant time.
     * 
     * The Fibonacci Heap has O(1) operations excepting
     * extractMin which is O(lg_2(N_nodes)).
     @param maxNumberOfBits
     * 
     */
    public MinHeapForRT2012(int maxValue, int approxN, int maxNumberOfBits) {

        //use yfasttrie if theres enough memory        
        long totalMemory = Runtime.getRuntime().totalMemory();
        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        long heapUsage = mbean.getHeapMemoryUsage().getUsed();
        long avail = totalMemory - heapUsage;

        long[] yftEstimate = YFastTrie.estimateSizeOnHeap(maxValue, 
                maxNumberOfBits);
        
        log.fine("avail=" + avail + " yftEst=" + yftEstimate[1] + " < " +
            (yftEstimate[1] < avail));
        
        if (yftEstimate[1] < avail) {
            // wanting the base of the prefix tree to be filled
            // to improve performance.   for larger N
            
            algorithm = 2;
            
            heap2 = new YFastTrieWrapper(maxValue);
            
            heap1 = null;
            
        } else {
            
            algorithm = 1;
        
            heap1 = new FibonacciHeapWrapper(approxN, maxValue);
        
            heap2 = null;
        }
    }
    
    /**
     * for fib. heap runtime is runtime is O(1).
     * for yft runtime complexity for best case is O(1) when there
     * is already a similar key in the YFastTrie, else is O(log_2(w)) + O(w-l)
     * where w is the number of bits set in the constructor and l is the prefix 
     * tree already filled leading up to the value node.  The performance of 
     * the YFastTrie increases when more nodes are in it (can see that in the
     * l term).
     @param node 
     */
    public void insert(HeapNode node) {
        
        if (node.getKey() < 0) {
            throw new IllegalArgumentException("key must be >= 0");
        }
        
        switch(algorithm) {
            case 1:
                // Fibonacci
                insert1(node);
                break;
            default:
                // YFastTrie
                insert2(node);
                break;
        }
    }
    
    /**
     if using a YFastTrie the runtime complexity is O(log log(M))
        where M is the number of bits of the maximum value the trie
        was initialized with.
     If using a FibonacciHeap, the runtime complexity is O(log_2 N) or better
       where N is the number of entries in the heap.
        
     @return extracted node which is the minimum in the queue
     */
    public HeapNode extractMin() {
        
        switch(algorithm) {
            case 1:
                return extractMin1();
            default:
                return extractMin2();
        }
        
    }
    
    /**
     * extract min-key from Fibonacci heap.
     * runtime is O(log_2 N) or better.
     @return 
     */
    private HeapNode extractMin1() {
        // runtime is O(log_2 N) or better
        HeapNode node = heap1.extractMin();
        if (node != null) {
            return node;
        } else {
            return null;
        }
    }
    
    /**
     * extract min-key from YFastTrie min priority queue.
     * runtime complexity is O(log log(M)) where M is the number of bits 
     * of the maximum value that the trie can hold (set during construction). 
     @return 
     */
    private HeapNode extractMin2() { 
        // runtime is runtime complexity is 
        return heap2.extractMin();
    }
    
    /**
     * insert into Fibonacci heap. r.t. complexity is O(1)
     @param node 
     */
    private void insert1(HeapNode node) {
         
        int key = (int)node.getKey();
        
        heap1.insert(node);
        
        log.fine("insert into fib minHeap at key =" + key);        
    }
    
    /**
     * insert node into YFastTrie min priority queue.
     * runtime complexity for best case is O(1) when there
     * is already a similar key in the YFastTrie, else is
     * is O(log_2(w)) + O(w-l)
     * where w is the number of bits set in the constructor
     * and l is the prefix tree already filled leading
     * up to the value x.  The performance of the YFastTrie
     * increases when more nodes are in it (can see that in the
     * l term).
     @param node
     */
    private void insert2(HeapNode node) {
         
        int key = (int)node.getKey();
        
        heap2.insert(node);
        
        log.fine("insert into yft minHeap at key =" + key);        
    }
    
    /**
     if using a YFastTrie the runtime complexity is O(log log(M))
        where M is the number of bits of the maximum value the trie
        was initialized with.
     If using a FibonacciHeap, the runtime complexity is O(1) if the node
     and the new key reside in same heap in the group of heaps, else
     the remove and insert makes the runtime complexity O(log_2(n)).
     @param node
     @param key2
     */
    public void decreaseKey(HeapNode node, long key2) {
    
        switch(algorithm) {
            case 1:
                decreaseKey1(node, key2);
                break;
            default:
                decreaseKey2(node, key2);
                break;
        }
    }
     
    private void decreaseKey1(HeapNode node, long key2) {

        log.fine("decreaseKey in fibHeap from key=" + 
            node.getKey() + " to key=" + key2);
        
        heap1.decreaseKey(node, key2);
    }
    
    private void decreaseKey2(HeapNode node, long key2) {

        log.fine("decreaseKey in yft from key=" + 
            node.getKey() + " to key=" + key2);
        
        heap2.decreaseKey(node, key2);
    }
    
     /**
     remove the node from the min-priority queue.
     If using a FibonacciHeap, the runtime complexity is O(log_2(n)).
     If using YFT, runtime complexity is O(log_2(w)) where w is the bitlength
     of the maximum value storable in the trie set at construction time.
     @param node
     */
    public void remove(HeapNode node) {
    
        switch(algorithm) {
            case 1:
                remove1(node);
                break;
            default:
                remove2(node);
                break;
        }
    }
     
    /**
     * remove from Fibonacci heap.
     * runtime complexity is O(lg_2(n))
     @param node 
     */
    private void remove1(HeapNode node) {

        log.fine("remove in fibHeap key=" + 
            node.getKey());
        
        heap1.remove(node);
    }
    
    /**
     * remove node from YFastTrie min priority queue.
     * runtime complexity is O(log log(M)) where M is the number of bits 
     * of the maximum value that the trie can hold (set during construction). 
     @param node
     */
    private void remove2(HeapNode node) {

        log.fine("remove in yft key=" + 
            node.getKey());
        
        heap2.remove(node);
    }
    
    /**
     *
     @return
     */
    public long getNumberOfNodes() {
        
        switch(algorithm) {
            case 1:
                return heap1.getNumberOfNodes();
            default:
                return heap2.getNumberOfNodes();
        }
    }

    @Override
    public String toString() {
    
        StringBuilder sb = new StringBuilder();
        sb.append("min heap type = ");
        switch(algorithm) {
            case 1:
                sb.append("Fibonacci Heap"); break;
            default:
                sb.append("YFastTrie min priority queue"); break;
        }
        sb.append(". size=").append(getNumberOfNodes());
        
        return sb.toString();
    }
    
    
}
