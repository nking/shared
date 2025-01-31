package algorithms.cache;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * a simple wrapping of java's LinkedHashMap in order to implement a "least
 * recently accessed cache".
 * 
 * @author nichole
 @param <K> type for key, e.g. String
 @param <V> type for value, e.g. Integer
 */
public class LRUCache<K, V> extends LinkedHashMap<K, V> {

    private static final long serialVersionUID = 1800104040800209030L;

    private static final boolean accessOrder = true;
    
    private final int initialCapacity;
    
    /**
     * default size is 3
     */
    public LRUCache() {
        super(3, (3 + 2.f)/3.f, accessOrder);
        this.initialCapacity = 3;
    }
    
    /**
     * create a fixed size LRU cache with loadFactor = (initialCapacity + 2.f)/initialCapacity.
     @param initialCapacity
     */
    public LRUCache(int initialCapacity) {
        super(initialCapacity, (initialCapacity + 2.f)/initialCapacity, accessOrder);
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
     * given the map, create a fixed size LRU cache with loadFactor = (m.size() + 2.f)/m.size().
     * The items are inserted without order using map iterator.
     @param m
     */
    public LRUCache(Map<? extends K, ? extends V> m) {
        super(m.size(), (m.size() + 2.f)/ m.size(), accessOrder);
        this.initialCapacity = m.size();
        putAll(m);
    }
    
    @Override
    protected boolean removeEldestEntry(Map.Entry<K, V> eldest) {
        return size() > initialCapacity;
    }

}
