package ctci;

import java.lang.reflect.Array;

/**
 *
 * @author nichole
 @param <K>
 @param <V>
 */
public class HashTable<K,V> {
  
    /**
     *
     */
    protected final LinkedList<K, V>[] t;

    /**
     *
     @param capacity
     */
    public HashTable(int capacity) {

      t = (LinkedList<K, V>[]) Array.newInstance(ctci.LinkedList.class, capacity);

      for (int i = 0; i < capacity; ++i) {
        t[i] = new LinkedList<K, V>();
      }
    }

    /**
     *
     @param key
     @param value
     */
    public void addItem(K key, V value) {
    int i = hash(key);
    t[i].add(key, value);
  }

    /**
     *
     @param key
     @return
     */
    public boolean containsKey(K key) {
    int i = hash(key);
    return t[i].containsKey(key);
  }

    /**
     *
     @param key
     @return
     */
    public V get(K key) {
    int i = hash(key);
    Node<K,V> node = t[i].get(key);
    return node.getData();
  }

    /**
     *
     @param key
     @return
     */
    public V remove(K key) {
    int i = hash(key);
    Node<K,V> node = t[i].remove(key);
    return node.getData();
  }

    /**
     *
     @param key
     @return
     */
    protected int hash(K key) {
    // this can be improved
    int hc = key.hashCode();
    int i = Math.floorMod(hc, t.length);
    return i;
  }
}
