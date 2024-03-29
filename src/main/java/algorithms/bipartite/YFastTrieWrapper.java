package algorithms.bipartite;

import algorithms.heapsAndPQs.YFastTrie;
import algorithms.heapsAndPQs.HeapNode;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * a wrapper for the YFastTrie to provide methods of a minimum
 * heap that can handle more than one node of the same key.
 * @author nichole
 */
public class YFastTrieWrapper {

    private final int w;
    
    private final int maxC;
    
    private final YFastTrie yft;
    
    private final TIntObjectHashMap<Set<HeapNode>> map =
        new TIntObjectHashMap<Set<HeapNode>>();
    
    private long lastKnownMinKey = 0;
    private long lastKnownMaxKey = -1;
        
    private int n = 0;
    
    /**
     *
     @param maxC
     */
    public YFastTrieWrapper(int maxC) {
            
        w = 1 + (int)Math.ceil(Math.log(maxC)/Math.log(2));
        
        this.maxC = maxC;
        
        yft = new YFastTrie(w);
        
    }
    
    /**
     *
     @return
     */
    public int getW() {
        return w;
    }
    
    /**
     * runtime complexity for best case is O(1) when there
     * is already a similar key in the underlying XFastTrie, else is
     * is O(log_2(w)) + O(w-l)
     * where w is the number of
     * bits set in the constructor
     * and l is the prefix tree already filled leading
     * up to the value x.  The performance of the XFastTrie
     * increases when more nodes are in it (can see that in the
     * l term).
     @param node 
     */
    public void insert(HeapNode node) {
        
        int keyIdx = (int)node.getKey();
        
        Set<HeapNode> set = map.get(keyIdx);
        
        if (set == null) {
            Integer key = Integer.valueOf(keyIdx);
            // O(log_2(w)) + O(w-l)
            boolean added = yft.add(key);
            assert(added);
            
            set = new HashSet<HeapNode>();
            map.put(keyIdx, set);
        }
        
        //O(1)
        set.add(node);
        
        n++;
        
        if (keyIdx < lastKnownMinKey) {
            lastKnownMinKey = keyIdx;
        }
        if (keyIdx > lastKnownMaxKey) {
            lastKnownMaxKey = keyIdx;
        }
    }
    
    /**
     * runtime complexity for best case is O(1) when there
     * is already a similar key in the XFastTrie, else is
     * is O(log_2(w)) + O(w-l)
     * where w is the number of
     * bits set in the constructor
     * and l is the prefix tree already filled leading
     * up to the value x.  The performance of the XFastTrie
     * increases when more nodes are in it (can see that in the
     * l term).
     @param node 
     @param key2 
     */
    public void decreaseKey(HeapNode node, long key2) {

        int keyIdx = (int)node.getKey();
                
        Set<HeapNode> set0 = map.get(keyIdx);
        
        assert(set0 != null);
        
        set0.remove(node);
        
        if (set0.isEmpty()) {
            boolean removed = yft.remove(Integer.valueOf(keyIdx));
            assert(removed);
            map.remove(keyIdx);
        }
                        
        node.setKey(key2);
        
        Integer index2 = Integer.valueOf((int)key2);
        
        Set<HeapNode> set2 = map.get((int)key2);
         
        if (set2 == null) {
            boolean added = yft.add(index2);
            assert(added);
            
            set2 = new HashSet<HeapNode>();
            map.put((int)key2, set2);
        }
        
        set2.add(node);
     
        if (key2 < lastKnownMinKey) {
            lastKnownMinKey = key2;
        }
    }
    
    /**
     * runtime complexity is roughly O(log_2(w)) where w is the bitlength of the
       maximum value stored in the trie, set at construction time.
     @param node 
     */
    public void remove(HeapNode node) {

        int keyIdx = (int)node.getKey();
                
        Set<HeapNode> set0 = map.get(keyIdx);
        
        assert(set0 != null);
        
        set0.remove(node);
        
        if (set0.isEmpty()) {
            boolean removed = yft.remove(Integer.valueOf(keyIdx));
            assert(removed);
            map.remove(keyIdx);
        }
        
        n--;
    }
    
    /**
     * runtime complexity for best case is O(1) when there
     * is already a similar key in the XFastTrie, else is
     * is O(log_2(w)) + O(w-l)
     * where w is the number of
     * bits set in the constructor
     * and l is the prefix tree already filled leading
     * up to the value x.  The performance of the XFastTrie
     * increases when more nodes are in it (can see that in the
     * l term).
     @return node 
     */
    public HeapNode extractMin() {
        
        if (n == 0) {
            return null;
        }
     
        Integer key = yft.minimum();
        
        Set<HeapNode> set = map.get(key.intValue());
        
        HeapNode node = set.iterator().next();
        set.remove(node);
        
        if (set.isEmpty()) {
            map.remove(key.intValue());
            yft.remove(key);
        }
       
        lastKnownMinKey = key.intValue();
        n--;
        
        return node;
    }
    
    /**
     *
     @return
     */
    public long getNumberOfNodes() {
        return n;
    }
}
