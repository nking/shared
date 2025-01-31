package ctci;

/**
 *
 * @author nichole
 @param <K> parameter type of key
 @param <V> parameter type of value
 */
public class LinkedList<K, V> {

    /**
     *
     */
    protected Node<K, V> head = null;

    /**
     *
     */
    protected Node<K, V> tail = null;

    /**
     *
     */
    protected int n = 0;

    /**
     *
     */
    public LinkedList() {
  }

    /**
     *
     @param key
     @param value
     */
    public LinkedList(K key, V value) {
    add(key, value);
  }

  /**
  add a node to the end of the linked list
     @param key
     @param value
     @return 
  */
  public Node<K,V> add(K key, V value) {
    Node<K,V> node = new Node(key, value);
    if (head == null) {
      assert(tail == null);
      head = node;
      tail = node;
      // node next and prev are null instead of pointing to self
    } else {
      // add to end of nodes
      node.prev = tail;
      tail.next = node;
      tail = node;
    }
    ++n;
    return node;
  }

    /**
     *
     @return
     */
    public Node<K,V> removeFirst() {
    //cases:
    //  (1) head == null and tail == null
    //  (2) there's only item in list so head == tail
    //  (3) there's more than one item in list

    if (head == null) {
      assert(n == 0);
      return null;
    }
    assert(tail != null);

    Node<K,V> node = head;

    // case where head == tail, that is, only 1 item was in the list
    if (head.equals(tail)) {
       assert(n == 1);
       head = null;
       tail = null;
       --n;
       // can remove next link of node to improve integrity of remaining LinkedList
       node.next = null;
       node.prev = null;
       return node;
    }

    assert(head.next != null);
    assert(tail.prev != null);
    assert(n > 1);

    head = head.next;
    head.prev = null;
    --n;

    // can remove next link of node to improve integrity of remaining LinkedList
    node.next = null;
    node.prev = null;

    return node;
  }

    /**
     *
     @return
     */
    public Node<K,V> removeLast() {
    //cases:
    //(1) head and tail are null
    //(2) head == tail, that is, only 1 item in list
    //(3) more than 1 item in list

    if (tail == null) {
      assert(head == null);
      assert(n == 0);
      return null;
    }

    Node<K,V> node = tail;

    if (head.equals(tail)) {
      assert(n == 1);
      head = null;
      tail = null;
      --n;
      // can remove next link of node to improve integrity of remaining LinkedList
      node.next = null;
      node.prev = null;
      return node;
    }

    tail = tail.prev;
    tail.next = null;
    --n;

    // can remove next link of node to improve integrity of remaining LinkedList
    node.next = null;
    node.prev = null;
    return node;
  }

    /**
     *
     @param key
     @return
     */
    public boolean containsKey(K key) {
    Node<K, V> node = get(key);
    return (node != null);
  }

    /**
     *
     @param key
     @return
     */
    public Node<K,V> get(K key) {
    // search from head to null
    Node<K, V> node = head;
    while (node != null) {
      if (node.getKey().equals(key)) {
        return node;
      }
      node = node.getNext();
    }
    return node;
  }

    /**
     *
     @param key
     @return
     */
    public Node<K,V> remove(K key) {
    Node<K, V> node = get(key);
    if (node == null) {
      return null;
    }
    // cases:
    // (1) key == head, can use removeFirst()
    if (node.equals(head)) {
      return removeFirst();
    }
    // (2) key == tail, can use removeLast()
    if (node.equals(tail)) {
      return removeLast();
    }
    // (3) key is not head nor tail
    // prev
    // node
    // tail
    node.prev.next = node.next;
    node.next.prev = node.prev;
    --n;

    // can remove next link of node to improve integrity of remaining LinkedList
    node.next = null;
    node.prev = null;
    return node;
  }

    /**
     *
     @return
     */
    public Node<K,V> first() {
    return head;
  }

    /**
     *
     @return
     */
    public Node<K,V> last() {
    return tail;
  }

    /**
     *
     @return
     */
    public int size() {
    return n;
  }

    /**
     *
     @return
     */
    public boolean isEmpty() {
    return (n == 0);
  }
}
