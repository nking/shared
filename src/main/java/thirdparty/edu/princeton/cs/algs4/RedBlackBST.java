package thirdparty.edu.princeton.cs.algs4;

/******************************************************************************
 a red black tree is a balanced binary tree which uses 2 colors in the
 nodes and maintains the following properties upon modification (insert and delete):
   (1) all leaf nodes are black;
   (2) a red node does not have a red child; and
   (3) every path from a node to all of its descendant leaves has the same 
       number of black nodes.
   (Cormen, Leiserson, Rivest, and Stein additionally use a black root node, but others do not.)
 The balanced binary trees have search time O(log_2(n)) and h = log_2(n).
 
   from RedBlackBST.java from the 
   book "Algorithms" by Sedgewick and Wayne
   http://algs4.cs.princeton.edu/33balanced/RedBlackBST.java
   copyright for authors Robert Sedgewick and Kevin Wayne
   is GPLV3, http://algs4.cs.princeton.edu/faq/

 *  Compilation:  javac RedBlackBST.java
 *  Execution:    java RedBlackBST left-pipe input.txt
 *  Dependencies: StdIn.java StdOut.java  
 *  Data files:   http://algs4.cs.princeton.edu/33balanced/tinyST.txt  
 *    
 *  A symbol table implemented using a left-leaning red-black BST.
 *  This is the 2-3 version.
 *
 *  Note: commented out assertions because DrJava now enables assertions
 *        by default.
 *
 ******************************************************************************/

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *  The {@code BST} class represents an ordered symbol table of generic
 *  key-value pairs.
 *  It supports the usual <em>put</em>, <em>get</em>, <em>contains</em>,
 *  <em>delete</em>, <em>size</em>, and <em>is-empty</em> methods.
 *  It also provides ordered methods for finding the <em>minimum</em>,
 *  <em>maximum</em>, <em>floor</em>, and <em>ceiling</em>.
 *  It also provides a <em>keys</em> method for iterating over all of the keys.
 *  A symbol table implements the <em>associative array</em> abstraction:
 *  when associating a value with a key that is already in the symbol table,
 *  the convention is to replace the old value with the new value.
 *  Unlike {@link java.util.Map}, this class uses the convention that
 *  values cannot be {@code null}—setting the
 *  value associated with a key to {@code null} is equivalent to deleting the key
 *  from the symbol table.
 *  <p>
 *  This implementation uses a left-leaning red-black BST. It requires that
 *  the key type implements the {@code Comparable} interface and calls the
 *  {@code compareTo()} and method to compare two keys. It does not call either
 *  {@code equals()} or {@code hashCode()}.
 *  The <em>put</em>, <em>contains</em>, <em>remove</em>, <em>minimum</em>,
 *  <em>maximum</em>, <em>ceiling</em>, and <em>floor</em> operations each take
 *  logarithmic time in the worst case, if the tree becomes unbalanced.
 *  The <em>size</em>, and <em>is-empty</em> operations take constant time.
 *  Construction takes constant time.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/33balanced">Section 3.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *  For other implementations of the same API, see {@link ST}, {@link BinarySearchST},
 *  {@link SequentialSearchST}, {@link BST},
 *  {@link SeparateChainingHashST}, {@link LinearProbingHashST}, and {@link AVLTreeST}.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 @param <Key> paramter type of key
 @param <Value> paramter type of value
 */

public class RedBlackBST<Key extends Comparable<Key>, Value> {

    private static final boolean RED   = true;
    private static final boolean BLACK = false;

    private RBNode root;     // root of the BST

    // BST helper node data type
    private class RBNode {
        private Key key;           // key
        private Value val;         // associated data
        private RBNode left, right;  // links to left and right subtrees
        private boolean color;     // color of parent link
        private int size;          // subtree count

        public RBNode(Key key, Value val, boolean color, int size) {
            this.key = key;
            this.val = val;
            this.color = color;
            this.size = size;
        }
        
         @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("key=").append(key);
            sb.append(" val=").append(val);
            sb.append(" color=").append(color);
            sb.append(" size=").append(size);
            sb.append(" l=");
            if (left != null) {
                sb.append(left.key);
            }
            sb.append(" r=");
            if (right != null) {
                sb.append(right.key);
            }
            return sb.toString();
        }
        
        public String toString(RBNode p) {
            StringBuilder sb = new StringBuilder();
            sb.append("key=").append(key);
            sb.append(" val=").append(val);
            sb.append(" color=").append(color);
            sb.append(" size=").append(size);
            sb.append(" p=");
            if (p != null) {
                sb.append(p.key);
            }
            sb.append(" l=");
            if (left != null) {
                sb.append(left.key);
            }
            sb.append(" r=");
            if (right != null) {
                sb.append(right.key);
            }
            return sb.toString();
        }
    }

    /**
     * Initializes an empty symbol table.
     */
    public RedBlackBST() {
    }

   /***************************************************************************
    *  RBNode helper methods.
     @param x
     @return 
    ***************************************************************************/
    // is node x red; false if x is null ?
    private boolean isRed(RBNode x) {
        if (x == null) return false;
        return x.color == RED;
    }

    // number of node in subtree rooted at x; 0 if x is null
    private int size(RBNode x) {
        if (x == null) return 0;
        return x.size;
    } 


    /**
     * Returns the number of key-value pairs in this symbol table.
     @return the number of key-value pairs in this symbol table
     */
    public int size() {
        return size(root);
    }

   /**
     * Is this symbol table empty?
     @return {@code true} if this symbol table is empty and {@code false} otherwise
     */
    public boolean isEmpty() {
        return root == null;
    }


   /***************************************************************************
    *  Standard BST search.
    ***************************************************************************/

    /**
     * Returns the value associated with the given key.
     @param key the key
     @return the value associated with the given key if the key is in the symbol table
     *     and {@code null} if the key is not in the symbol table
     * @throws IllegalArgumentException if {@code key} is {@code null}
     */
    public Value get(Key key) {
        if (key == null) throw new IllegalArgumentException("argument to get() is null");
        return get(root, key);
    }

    // value associated with the given key in subtree rooted at x; null if no such key
    private Value get(RBNode x, Key key) {
        while (x != null) {
            int cmp = key.compareTo(x.key);
            if      (cmp < 0) x = x.left;
            else if (cmp > 0) x = x.right;
            else              return x.val;
        }
        return null;
    }

    /**
     * Does this symbol table contain the given key?
     @param key the key
     @return {@code true} if this symbol table contains {@code key} and
     *     {@code false} otherwise
     * @throws IllegalArgumentException if {@code key} is {@code null}
     */
    public boolean contains(Key key) {
        return get(key) != null;
    }

   /***************************************************************************
    *  Red-black tree insertion.
    ***************************************************************************/

    /**
     * Inserts the specified key-value pair into the symbol table, overwriting the old 
     * value with the new value if the symbol table already contains the specified key.
     * Deletes the specified key (and its associated value) from this symbol table
     * if the specified value is {@code null}.
     *
     @param key the key
     @param val the value
     * @throws IllegalArgumentException if {@code key} is {@code null}
     */
    public void put(Key key, Value val) {
        if (key == null) throw new IllegalArgumentException("first argument to put() is null");
        if (val == null) {
            delete(key);
            return;
        }

        root = put(root, key, val);
        root.color = BLACK;
        // assert check();
    }

    // insert the key-value pair in the subtree rooted at h
    private RBNode put(RBNode h, Key key, Value val) { 
        if (h == null) return new RBNode(key, val, RED, 1);

        int cmp = key.compareTo(h.key);
        if (cmp < 0) {
            h.left  = put(h.left,  key, val);
        } else if (cmp > 0) {
            h.right = put(h.right, key, val);
        } else {
            h.val   = val;
        }

        // fix-up any right-leaning links
        if (isRed(h.right) && !isRed(h.left))      h = rotateLeft(h);
        if (isRed(h.left)  &&  isRed(h.left.left)) h = rotateRight(h);
        if (isRed(h.left)  &&  isRed(h.right))     flipColors(h);
        h.size = size(h.left) + size(h.right) + 1;

        return h;
    }

   /***************************************************************************
    *  Red-black tree deletion.
    ***************************************************************************/

    /**
     * Removes the smallest key and associated value from the symbol table.
     * @throws NoSuchElementException if the symbol table is empty
     */
    public void deleteMin() {
        if (isEmpty()) throw new NoSuchElementException("BST underflow");

        // if both children of root are black, set root to red
        if (!isRed(root.left) && !isRed(root.right))
            root.color = RED;

        root = deleteMin(root);
        if (!isEmpty()) root.color = BLACK;
        // assert check();
    }

    // delete the key-value pair with the minimum key rooted at h
    private RBNode deleteMin(RBNode h) { 
        if (h.left == null)
            return null;

        if (!isRed(h.left) && !isRed(h.left.left))
            h = moveRedLeft(h);

        h.left = deleteMin(h.left);
        return balance(h);
    }


    /**
     * Removes the largest key and associated value from the symbol table.
     * @throws NoSuchElementException if the symbol table is empty
     */
    public void deleteMax() {
        if (isEmpty()) throw new NoSuchElementException("BST underflow");

        // if both children of root are black, set root to red
        if (!isRed(root.left) && !isRed(root.right))
            root.color = RED;

        root = deleteMax(root);
        if (!isEmpty()) root.color = BLACK;
        // assert check();
    }

    // delete the key-value pair with the maximum key rooted at h
    private RBNode deleteMax(RBNode h) { 
        if (isRed(h.left))
            h = rotateRight(h);

        if (h.right == null)
            return null;

        if (!isRed(h.right) && !isRed(h.right.left))
            h = moveRedRight(h);

        h.right = deleteMax(h.right);

        return balance(h);
    }

    /**
     * Removes the specified key and its associated value from this symbol table     
     * (if the key is in this symbol table).    
     *
     @param  key the key
     * @throws IllegalArgumentException if {@code key} is {@code null}
     */
    public void delete(Key key) { 
        if (key == null) throw new IllegalArgumentException("argument to delete() is null");
        if (!contains(key)) return;

        // if both children of root are black, set root to red
        if (!isRed(root.left) && !isRed(root.right))
            root.color = RED;

        root = delete(root, key);
        if (!isEmpty()) root.color = BLACK;
        // assert check();
        
    }

    // delete the key-value pair with the given key rooted at h
    private RBNode delete(RBNode h, Key key) { 
        // assert get(h, key) != null;

        if (key.compareTo(h.key) < 0)  {
            if (!isRed(h.left) && !isRed(h.left.left))
                h = moveRedLeft(h);
            h.left = delete(h.left, key);
        }
        else {
            if (isRed(h.left))
                h = rotateRight(h);
            if (key.compareTo(h.key) == 0 && (h.right == null))
                return null;
            if (!isRed(h.right) && !isRed(h.right.left))
                h = moveRedRight(h);
            if (key.compareTo(h.key) == 0) {
                RBNode x = min(h.right);
                
                h.key = x.key;
                h.val = x.val;
                
                // h.val = get(h.right, min(h.right).key);
                // h.key = min(h.right).key;
                h.right = deleteMin(h.right);
            }
            else h.right = delete(h.right, key);
        }
        return balance(h);
    }

   /***************************************************************************
    *  Red-black tree helper functions.
     @param h
     @return 
    ***************************************************************************/

    // make a left-leaning link lean to the right
    private RBNode rotateRight(RBNode h) {
        // assert (h != null) && isRed(h.left);
        RBNode x = h.left;
        h.left = x.right;
        x.right = h;
        x.color = x.right.color;
        x.right.color = RED;
        x.size = h.size;
        h.size = size(h.left) + size(h.right) + 1;
        return x;
    }

    // make a right-leaning link lean to the left
    private RBNode rotateLeft(RBNode h) {
        // assert (h != null) && isRed(h.right);
        RBNode x = h.right;
        h.right = x.left;
        x.left = h;
        x.color = x.left.color;
        x.left.color = RED;
        x.size = h.size;
        h.size = size(h.left) + size(h.right) + 1;
        return x;
    }

    // flip the colors of a node and its two children
    private void flipColors(RBNode h) {
        // h must have opposite color of its two children
        // assert (h != null) && (h.left != null) && (h.right != null);
        // assert (!isRed(h) &&  isRed(h.left) &&  isRed(h.right))
        //    || (isRed(h)  && !isRed(h.left) && !isRed(h.right));
        h.color = !h.color;
        h.left.color = !h.left.color;
        h.right.color = !h.right.color;
    }

    // Assuming that h is red and both h.left and h.left.left
    // are black, make h.left or one of its children red.
    private RBNode moveRedLeft(RBNode h) {
        // assert (h != null);
        // assert isRed(h) && !isRed(h.left) && !isRed(h.left.left);

        flipColors(h);
        if (isRed(h.right.left)) { 
            h.right = rotateRight(h.right);
            h = rotateLeft(h);
            flipColors(h);
        }
        return h;
    }

    // Assuming that h is red and both h.right and h.right.left
    // are black, make h.right or one of its children red.
    private RBNode moveRedRight(RBNode h) {
        // assert (h != null);
        // assert isRed(h) && !isRed(h.right) && !isRed(h.right.left);
        flipColors(h);
        if (isRed(h.left.left)) { 
            h = rotateRight(h);
            flipColors(h);
        }
        return h;
    }

    // restore red-black tree invariant
    private RBNode balance(RBNode h) {
        // assert (h != null);

        if (isRed(h.right))                      h = rotateLeft(h);
        if (isRed(h.left) && isRed(h.left.left)) h = rotateRight(h);
        if (isRed(h.left) && isRed(h.right))     flipColors(h);

        h.size = size(h.left) + size(h.right) + 1;
        return h;
    }


   /***************************************************************************
    *  Utility functions.
    ***************************************************************************/

    /**
     * Returns the height of the BST (for debugging).
     @return the height of the BST (a 1-node tree has height 0)
     */
    public int height() {
        return height(root);
    }
    private int height(RBNode x) {
        if (x == null) return -1;
        return 1 + Math.max(height(x.left), height(x.right));
    }

   /***************************************************************************
    *  Ordered symbol table methods.
    ***************************************************************************/

    /**
     * Returns the smallest key in the symbol table.
     @return the smallest key in the symbol table
     * @throws NoSuchElementException if the symbol table is empty
     */
    public Key min() {
        if (isEmpty()) throw new NoSuchElementException("called min() with empty symbol table");
        return min(root).key;
    } 

    // the smallest key in subtree rooted at x; null if no such key
    private RBNode min(RBNode x) { 
        // assert x != null;
        if (x.left == null) return x; 
        else                return min(x.left); 
    } 

    /**
     * Returns the largest key in the symbol table.
     @return the largest key in the symbol table
     * @throws NoSuchElementException if the symbol table is empty
     */
    public Key max() {
        if (isEmpty()) throw new NoSuchElementException("called max() with empty symbol table");
        return max(root).key;
    } 

    // the largest key in the subtree rooted at x; null if no such key
    private RBNode max(RBNode x) { 
        // assert x != null;
        if (x.right == null) return x; 
        else                 return max(x.right); 
    } 


    /**
     * Returns the largest key in the symbol table less than or equal to {@code key}.
     @param key the key
     @return the largest key in the symbol table less than or equal to {@code key}
     * @throws NoSuchElementException if there is no such key
     * @throws IllegalArgumentException if {@code key} is {@code null}
     */
    public Key floor(Key key) {
        if (key == null) throw new IllegalArgumentException("argument to floor() is null");
        if (isEmpty()) throw new NoSuchElementException("called floor() with empty symbol table");
        RBNode x = floor(root, key);
        if (x == null) return null;
        else           return x.key;
    }    

    // the largest key in the subtree rooted at x less than or equal to the given key
    private RBNode floor(RBNode x, Key key) {
        if (x == null) return null;
        int cmp = key.compareTo(x.key);
        if (cmp == 0) return x;
        if (cmp < 0)  return floor(x.left, key);
        RBNode t = floor(x.right, key);
        if (t != null) return t; 
        else           return x;
    }

    /**
     * Returns the smallest key in the symbol table greater than or equal to {@code key}.
     @param key the key
     @return the smallest key in the symbol table greater than or equal to {@code key}
     * @throws NoSuchElementException if there is no such key
     * @throws IllegalArgumentException if {@code key} is {@code null}
     */
    public Key ceiling(Key key) {
        if (key == null) throw new IllegalArgumentException("argument to ceiling() is null");
        if (isEmpty()) throw new NoSuchElementException("called ceiling() with empty symbol table");
        RBNode x = ceiling(root, key);
        if (x == null) return null;
        else           return x.key;  
    }

    // the smallest key in the subtree rooted at x greater than or equal to the given key
    private RBNode ceiling(RBNode x, Key key) {  
        if (x == null) return null;
        int cmp = key.compareTo(x.key);
        if (cmp == 0) return x;
        if (cmp > 0)  return ceiling(x.right, key);
        RBNode t = ceiling(x.left, key);
        if (t != null) return t; 
        else           return x;
    }

    /**
     * Return the kth smallest key in the symbol table.
     @param k the order statistic
     @return the {@code k}th smallest key in the symbol table
     * @throws IllegalArgumentException unless {@code k} is between 0 and
     *     <em>n</em>–1
     */
    public Key select(int k) {
        if (k < 0 || k >= size()) {
            throw new IllegalArgumentException("called select() with invalid argument: " + k);
        }
        RBNode x = select(root, k);
        return x.key;
    }

    // the key of rank k in the subtree rooted at x
    private RBNode select(RBNode x, int k) {
        // assert x != null;
        // assert k >= 0 && k < size(x);
        int t = size(x.left); 
        if      (t > k) return select(x.left,  k); 
        else if (t < k) return select(x.right, k-t-1); 
        else            return x; 
    } 

    /**
     * Return the number of keys in the symbol table strictly less than {@code key}.
     @param key the key
     @return the number of keys in the symbol table strictly less than {@code key}
     * @throws IllegalArgumentException if {@code key} is {@code null}
     */
    public int rank(Key key) {
        if (key == null) throw new IllegalArgumentException("argument to rank() is null");
        return rank(key, root);
    } 

    // number of keys less than key in the subtree rooted at x
    private int rank(Key key, RBNode x) {
        if (x == null) return 0; 
        int cmp = key.compareTo(x.key); 
        if      (cmp < 0) return rank(key, x.left); 
        else if (cmp > 0) return 1 + size(x.left) + rank(key, x.right); 
        else              return size(x.left); 
    } 

   /***************************************************************************
    *  Range count and range search.
    ***************************************************************************/

    /**
     * Returns all keys in the symbol table as an {@code Iterable}.
     * To iterate over all of the keys in the symbol table named {@code st},
     * use the foreach notation: {@code for (Key key : st.keys())}.
     @return all keys in the symbol table as an {@code Iterable}
     */
    public Iterable<Key> keys() {
        if (isEmpty()) return new Queue<Key>();
        return keys(min(), max());
    }

    /**
     * Returns all keys in the symbol table in the given range,
     * as an {@code Iterable}.
     *
     @param  lo minimum endpoint
     @param  hi maximum endpoint
     @return all keys in the sybol table between {@code lo} 
     *    (inclusive) and {@code hi} (inclusive) as an {@code Iterable}
     * @throws IllegalArgumentException if either {@code lo} or {@code hi}
     *    is {@code null}
     */
    public Iterable<Key> keys(Key lo, Key hi) {
        if (lo == null) throw new IllegalArgumentException("first argument to keys() is null");
        if (hi == null) throw new IllegalArgumentException("second argument to keys() is null");

        Queue<Key> queue = new Queue<Key>();
        // if (isEmpty() || lo.compareTo(hi) > 0) return queue;
        keys(root, queue, lo, hi);
        return queue;
    } 

    // add the keys between lo and hi in the subtree rooted at x
    // to the queue
    private void keys(RBNode x, Queue<Key> queue, Key lo, Key hi) { 
        if (x == null) return; 
        int cmplo = lo.compareTo(x.key); 
        int cmphi = hi.compareTo(x.key); 
        if (cmplo < 0) keys(x.left, queue, lo, hi); 
        if (cmplo <= 0 && cmphi >= 0) queue.enqueue(x.key); 
        if (cmphi > 0) keys(x.right, queue, lo, hi); 
    } 

    /**
     * Returns the number of keys in the symbol table in the given range.
     *
     @param  lo minimum endpoint
     @param  hi maximum endpoint
     @return the number of keys in the sybol table between {@code lo} 
     *    (inclusive) and {@code hi} (inclusive)
     * @throws IllegalArgumentException if either {@code lo} or {@code hi}
     *    is {@code null}
     */
    public int size(Key lo, Key hi) {
        if (lo == null) throw new IllegalArgumentException("first argument to size() is null");
        if (hi == null) throw new IllegalArgumentException("second argument to size() is null");

        if (lo.compareTo(hi) > 0) return 0;
        if (contains(hi)) return rank(hi) - rank(lo) + 1;
        else              return rank(hi) - rank(lo);
    }


   /***************************************************************************
    *  Check integrity of red-black tree data structure.
     @return 
    ***************************************************************************/
    protected boolean check() {
        if (!isBST())            System.out.println("Not in symmetric order");
        if (!isSizeConsistent()) System.out.println("Subtree counts not consistent");
        if (!isRankConsistent()) System.out.println("Ranks not consistent");
        if (!is23())             System.out.println("Not a 2-3 tree");
        if (!isBalanced())       System.out.println("Not balanced");
        return isBST() && isSizeConsistent() && isRankConsistent() && is23() && isBalanced();
    }

    // does this binary tree satisfy symmetric order?
    // Note: this test also ensures that data structure is a binary tree since order is strict
    private boolean isBST() {
        return isBST(root, null, null);
    }

    // is the tree rooted at x a BST with all keys strictly between min and max
    // (if min or max is null, treat as empty constraint)
    // Credit: Bob Dondero's elegant solution
    private boolean isBST(RBNode x, Key min, Key max) {
        if (x == null) return true;
        if (min != null && x.key.compareTo(min) <= 0) return false;
        if (max != null && x.key.compareTo(max) >= 0) return false;
        return isBST(x.left, min, x.key) && isBST(x.right, x.key, max);
    } 

    // are the size fields correct?
    private boolean isSizeConsistent() { return isSizeConsistent(root); }
    private boolean isSizeConsistent(RBNode x) {
        if (x == null) return true;
        if (x.size != size(x.left) + size(x.right) + 1) return false;
        return isSizeConsistent(x.left) && isSizeConsistent(x.right);
    } 

    // check that ranks are consistent
    private boolean isRankConsistent() {
        for (int i = 0; i < size(); i++)
            if (i != rank(select(i))) return false;
        for (Key key : keys())
            if (key.compareTo(select(rank(key))) != 0) return false;
        return true;
    }

    // Does the tree have no red right links, and at most one (left)
    // red links in a row on any path?
    private boolean is23() { return is23(root); }
    private boolean is23(RBNode x) {
        if (x == null) return true;
        if (isRed(x.right)) return false;
        if (x != root && isRed(x) && isRed(x.left))
            return false;
        return is23(x.left) && is23(x.right);
    } 

    // do all paths from root to leaf have same number of black edges?
    private boolean isBalanced() { 
        int black = 0;     // number of black links on path from root to min
        RBNode x = root;
        while (x != null) {
            if (!isRed(x)) black++;
            x = x.left;
        }
        return isBalanced(root, black);
    }

    // does every path from the root to a leaf have the given number of black links?
    private boolean isBalanced(RBNode x, int black) {
        if (x == null) return black == 0;
        if (!isRed(x)) black--;
        return isBalanced(x.left, black) && isBalanced(x.right, black);
    } 


    /**
     * Unit tests the {@code RedBlackBST} data type.
     *
     @param args the command-line arguments
     */
    public static void main(String[] args) { 
        RedBlackBST<String, Integer> st = new RedBlackBST<String, Integer>();
        int rd = -1;
        StringBuilder sb = new StringBuilder();
        int i = 0; 
        try {
            do {
                rd = System.in.read();
                if (rd > -1) {
                    sb.append(Byte.toString((byte)rd));
                    ++i;
                    st.put(sb.toString(), Integer.valueOf(i)); 
                }
            } while (rd != -1);
        } catch (IOException ex) {
            Logger.getLogger(RedBlackBST.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    /**
     * left subtree, root, right subtree
     */
    public void printInOrderTraversal() {
        List<RBNode> nodes = getInOrderTraversalIterative(root);
        for (RBNode node : nodes) {
            System.out.println("node=" + node.toString());
        }
    }
    
    /**
     * root, left subtree, right subtree
     */
    public void printPreOrderTraversal() {
        List<RBNode> nodes = getPreOrderTraversalIterative(root, 0);
        for (RBNode node : nodes) {
            RBNode p = findParent(node, nodes);
            System.out.println("node=" + node.toString(p));
        }
    }

    public void printLevelOrderTraversal() {
        List<RBNode> nodes = getLevelOrderTraversalIterative(root);
        for (RBNode node : nodes) {
            RBNode p = findParent(node, nodes);
            System.out.println("node=" + node.toString(p));
        }
    }

    /**
     *
     @param topNode
     */
    public void printPreOrderTraversal2(RBNode topNode) {
        List<RBNode> nodes = getPreOrderTraversalIterative(topNode, 0);
        for (RBNode node : nodes) {
            RBNode p = findParent(node, nodes);
            System.out.println("node=" + node.toString(p));
        }
    }
    private void printPreOrderTraversal(int addExtraToSize) {
        List<RBNode> nodes = getPreOrderTraversalIterative(root, addExtraToSize);
        for (RBNode node : nodes) {
            if (node == null) {
                continue;
            }
            RBNode p = findParent(node, nodes);
            System.out.println("node=" + node.toString(p));
        }
    }
    
    private RBNode findParent(RBNode h) {
        List<RBNode> nodes = getPreOrderTraversalIterative(root, 1);
        return findParent(h, nodes);
    }
    private RBNode findParent(RBNode h, List<RBNode> nodes) {
        for (RBNode node : nodes) {
            if (node == null) {
                continue;
            }
            if (node.left != null && node.left.key == h.key) {
                return node;
            }
            if (node.right != null && node.right.key == h.key) {
                return node;
            }
        }
        return null;
    }
    
    /**
     * left subtree, right subtree, root subtree
     */
    public void printPostOrderTraversal() {
        List<RBNode> nodes = getPostOrderTraversalIterative(root);
        for (RBNode node : nodes) {
            System.out.println("node=" + node.toString());
        }
    }

    /**
    visit each node using pattern left subtree, root, right subtree
    in an iterative manner rather than invoking the method recursively.
    * <pre>
    * example:
                   10
              5           20
                        30
        
            returns [5,10,30,20]
    </pre>
    <pre>
     from "Data Structures and Algorithms", Garnett, and Del Tongo 2008.
     The recursive version is:
         if root ̸= null {
             inorder(root.Left) 
             yield root.Value 
             inorder(root.Right)
         }
     </pre>
     @param node
     @return 
    */
    protected List<RBNode> getInOrderTraversalIterative(RBNode node) {
       
        if (isEmpty()) {
            return new ArrayList<RBNode>();
        }
        
        int sz = size();
       
        List<RBNode> out = new ArrayList<RBNode>();;
        int count = 0;
        
        Stack<RBNode> stack = new Stack<>();
               
        while (!stack.isEmpty() || (node != null)) {
            if (node != null) {
                 
                stack.push(node);
                
                node = node.left;
            
            } else if (count < sz) {
                
                node = stack.pop();
                
                out.add(node);
                count++;
                
                //System.out.println(node.key);
                
                node = node.right;
            }
        }
        return out;
    }
    
    /**
    visit each node using pattern root node, then all direct children of root node (=level 2),
    then all direct children of those children (=level 3), etc
    in an iterative manner.
     @param node
     @return 
    */
    protected List<RBNode> getLevelOrderTraversalIterative(RBNode node) {

        List<RBNode> array = new ArrayList<RBNode>();
        int count = 0;
        
        if (node == null) {
            return array;
        }
        
        java.util.Queue<RBNode> q = new ArrayDeque<RBNode>();
        q.add(node);
        while (!q.isEmpty()) {
            node = q.poll();
            array.add(node);
            if (node.left != null) {
                q.add(node.left);
            }
            if (node.right != null) {
                q.add(node.right);
            }
        }

        return array;
    }
    
    /**
     * visit each node using pattern: root, left subtree, right subtree
     * in an iterative manner rather than invoking the method recursively.
     <pre>
     from "Data Structures and Algorithms", Garnett, and Del Tongo 2008.
     The recursive version is:
         if root ̸= null {
             yield root.Value 
             preorder(root.Left) 
             preorder(root.Right)
         }
     </pre>
     @param node
     @param addExtraToSize
     @return 
     */
    protected List<RBNode> getPreOrderTraversalIterative(RBNode node, int addExtraToSize) {
       
        //NOTE: added additional integer and conditions 
        //   for size because may be printing tree
        //   in the middle of a put where the node size is not yet updated.
        // The count conditionals below are otherwise, not needed.
        
        List<RBNode> out = new ArrayList<RBNode>();
        
        if (isEmpty()) {
            return out;
        }
                
        int sz = size() + addExtraToSize;
        
        int count = 0;
        
        Stack<RBNode> stack = new Stack<>();
        
        while (count < sz && (!stack.isEmpty() || node != null)) {
            if (node != null) {
              
                out.add(node);
                
                count++;
                
                if (count < sz) {
                    stack.push(node);
                }
                
                node = node.left;
            
            } else {
                
                node = stack.pop();
                
                node = node.right;
            }
        }
        
        return out;
    }

    /**
     * visit each node using pattern: left subtree, right subtree, root subtree
     * in an iterative manner rather than invoking the method recursively.
     <pre>
       example:
                 10
              5           20
                        30
        
        gives out=[5, 30, 20, 10]
     </pre>
     <pre>
     from "Data Structures and Algorithms", Garnett, and Del Tongo 2008.
     The recursive version is:
         if root ̸= null {
             postorder(root.Left) 
             postorder(root.Right)
             yield root.Value
         }
     </pre>
     @param node
     @return 
     */
    protected List<RBNode> getPostOrderTraversalIterative(RBNode node) {
        
        if (isEmpty()) {
            return new ArrayList<RBNode>();
        }
        
        List<RBNode> out = new ArrayList<RBNode>();
        int count = 0;
        
        if (node == null) {
            return out;
        }
        
        Stack<RBNode> stack = new Stack<>();
        Stack<RBNode> stack2 = new Stack<>();
        stack.push(node);
        
        while (!stack.isEmpty()) {
            
            node = stack.pop();
            
            stack2.push(node);
            
            if (node.left != null) {
                stack.push(node.left);
            }

            if (node.right != null) {
                stack.push(node.right);
            }            
        }
        
        while (!stack2.isEmpty()) {
            
            node = stack2.pop();
            
            //process(node);
            out.add(node);
            count++;
            //System.out.println(node.key);
        }
         
        return out;
    }
    
}
