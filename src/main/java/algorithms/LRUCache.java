package algorithms;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * a simple wrapping of java's LinkedHashMap in order to implement a "least
 * recently accessed cache".
 * 
 * @author nichole
 */
public class LRUCache<K, V> extends LinkedHashMap<K, V> {
        
    private static final boolean accessOrder = true;
    
    private final int initialCapacity;
    
    public LRUCache() {
        super(3, 0.75f, accessOrder);
        this.initialCapacity = 3;
    }
    
    public LRUCache(int initialCapacity) {
        super(initialCapacity, 0.75f, accessOrder);
        this.initialCapacity = initialCapacity;
    }
    
    public LRUCache(int initialCapacity, float loadFactor) {
        super(initialCapacity, loadFactor, accessOrder);
        this.initialCapacity = initialCapacity;
    }
    
    public LRUCache(int initialCapacity, float loadFactor, boolean accessOrder) {
        super(initialCapacity, loadFactor, accessOrder);
        if (!accessOrder) {
            throw new IllegalArgumentException("for the Least Recently Accessed Cache, need accessOrder=true");
        }
        this.initialCapacity = initialCapacity;
    }
    
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
