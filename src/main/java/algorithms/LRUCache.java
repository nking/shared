package algorithms;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * a simple wrapping of java's LinkedHashMap in order to implement a "least
 * recently accessed cache".
 * 
 * @author nichole
 @param <K>
 @param <V>
 */
public class LRUCache<K, V> extends LinkedHashMap<K, V> {

    private static final long serialVersionUID = 1800104040800209030L;

    private static final boolean accessOrder = true;
    
    private final int initialCapacity;
    
    /**
     *
     */
    public LRUCache() {
        super(3, 0.75f, accessOrder);
        this.initialCapacity = 3;
    }
    
    /**
     *
     @param initialCapacity
     */
    public LRUCache(int initialCapacity) {
        super(initialCapacity, 0.75f, accessOrder);
        this.initialCapacity = initialCapacity;
    }
    
    /**
     *
     @param initialCapacity
     @param loadFactor
     */
    public LRUCache(int initialCapacity, float loadFactor) {
        super(initialCapacity, loadFactor, accessOrder);
        this.initialCapacity = initialCapacity;
    }
    
    /**
     *
     @param initialCapacity
     @param loadFactor
     @param accessOrder
     */
    public LRUCache(int initialCapacity, float loadFactor, boolean accessOrder) {
        super(initialCapacity, loadFactor, accessOrder);
        if (!accessOrder) {
            throw new IllegalArgumentException("for the Least Recently Accessed Cache, need accessOrder=true");
        }
        this.initialCapacity = initialCapacity;
    }
    
    /**
     *
     @param m
     */
    public LRUCache(Map<? extends K, ? extends V> m) {
        super(3, 0.75f, accessOrder);
        this.initialCapacity = 3;
        putAll(m);
    }
    
    @Override
    protected boolean removeEldestEntry(Map.Entry<K, V> eldest) {
        return size() > initialCapacity;
    }

}
