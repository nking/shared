package ctci;

/**
 *
 * @author nichole
 @param <K>
 @param <V>
 */
public class Node<K,V> {

    /**
     *
     */
    protected final K key;

    /**
     *
     */
    protected final V data;

    /**
     *
     */
    protected Node<K,V> prev = null;

    /**
     *
     */
    protected Node<K,V> next = null;

    /**
     *
     @param key
     @param value
     */
    public Node(K key, V value) {
    this.key = key;
    this.data = value;
  }

    /**
     *
     @param node
     */
    public void setNext(Node<K,V> node) {
    this.next = node;
  }

    /**
     *
     @param node
     */
    public void setPrev(Node<K,V> node) {
    this.prev = node;
  }

    /**
     *
     @return
     */
    public K getKey() {
    return key;
  }

    /**
     *
     @return
     */
    public V getData() {
    return data;
  }
  
    /**
     *
     @return
     */
    public Node<K,V> getNext() {
    return next;
  }

    /**
     *
     @return
     */
    public Node<K,V> getPrev() {
    return prev;
  }
}