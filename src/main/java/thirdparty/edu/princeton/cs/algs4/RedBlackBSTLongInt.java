package thirdparty.edu.princeton.cs.algs4;

/******************************************************************************
   adapted from RedBlackBST.java 
   from the 
   book "Algorithms" by Sedgewick and Wayne
   http://algs4.cs.princeton.edu/33balanced/RedBlackBST.java
   copyright for authors Robert Sedgewick and Kevin Wayne
   is GPLV3, http://algs4.cs.princeton.edu/faq/

*     x.left.key .lte. x.key
*     x.right.key .gte. x.key
* 
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
 *  % more tinyST.txt
 *  S E A R C H E X A M P L E
 *  
 *  % java RedBlackBST left-pipe tinyST.txt
 *  A 8
 *  C 4
 *  E 12
 *  H 5
 *  L 11
 *  M 9
 *  P 10
 *  R 3
 *  S 0
 *  X 7
 *
 ******************************************************************************/

import algorithms.util.ObjectSpaceEstimator;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TLongArrayList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Stack;

/*
note that a reduction in space, especically for large number of nodes,
can be made by replacing nodes with associative arrays, but one would
need to include a parent array also and additional logic to
update the relationships where changes occur for a 
key, left or right field.
*/

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
 */

public class RedBlackBSTLongInt {
    
    private static final boolean RED   = true;
    private static final boolean BLACK = false;

    private Node root;     // root of the BST
    
    // BST helper node data type
    private static class Node {
        // changing Key to long primitive saves factor of 3 on 32 bit platforms,
        //    else factor of 2 on 64 bit platforms.
        //    (long size on stack is 64 bits, or 128 bits, respectively)
        //    Object Long is 16Bytes overhead + contents.
        // changing value from Integer to int saves factor of 5 on 32 bit
        //    platforms else 3 on 64 bit platforms
        private long key;           // key
        private int val;         // associated data
        private Node left, right;  // links to left and right subtrees
        private boolean color;     // color of parent link
        private int size;          // subtree count

        public Node(long key, int val, boolean color, int size) {
            this.key = key;
            this.val = val;
            this.color = color;
            this.size = size;
        }
        
        public static long estimateSizeOnHeap() {
            ObjectSpaceEstimator est = new ObjectSpaceEstimator();
            est.setNBooleanFields(1);
            est.setNLongFields(1);
            est.setNIntFields(2);
            est.setNObjRefsFields(2);
            return est.estimateSizeOnHeap();
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
        
        public String toString(Node p) {
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
    public RedBlackBSTLongInt() {
    }

   /***************************************************************************
    *  Node helper methods.
     @param x
     @return 
    ***************************************************************************/
    // is node x red; false if x is null ?
    private boolean isRed(Node x) {
        if (x == null) return false;
        return x.color == RED;
    }

    // number of node in subtree rooted at x; 0 if x is null
    private int size(Node x) {
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
     @param output if output[0] == -1, then key was not present, else the
     *    returned value is found in output[1]
     */
    public void get(long key, int[] output) {
        if (output == null || output.length != 2) {
            throw new IllegalArgumentException("output must be length 2");
        }
        output[0] = 0;
        get(root, key, output);
    }

    // value associated with the given key in subtree rooted at x; null if no such key
    private void get(Node x, long key, int[] output) {
        while (x != null) {
            int cmp = key < x.key ? -1 : (key > x.key) ? 1 : 0; 
            if      (cmp < 0) x = x.left;
            else if (cmp > 0) x = x.right;
            else {
                output[1] = x.val;
                return;
            }
        }
        output[0] = -1;
    }

    /**
     * Does this symbol table contain the given key?
     @param key the key
     @return {@code true} if this symbol table contains {@code key} and
     *     {@code false} otherwise
     * @throws IllegalArgumentException if {@code key} is {@code null}
     */
    public boolean contains(long key) {
        get(key, cache0);
        return (cache0[0] != -1);
    }
    private int[] cache0 = new int[2];

   /***************************************************************************
    *  Red-black tree insertion.
    ***************************************************************************/

    /**
     * Inserts the specified key-value pair into the symbol table, overwriting the old 
     * value with the new value if the symbol table already contains the specified key.
     *
     @param key the key
     @param val the value
     */
    public void put(long key, int val) {
        
        //System.out.println("put " + key);
        
        root = put(root, key, val);
        root.color = BLACK;
        
        //printPreOrderTraversal();
        //System.out.println("after put " + key + " root=" 
        //    + root.toString(findParent(root)));
                
        // assert check();
    }

    // insert the key-value pair in the subtree rooted at h
    private Node put(Node h, long key, int val) { 
        
        //System.out.println("put h=" + h + " key=" + key);
        
        if (h == null) return new Node(key, val, RED, 1);

        int cmp = (key < h.key) ? -1 : (key > h.key) ? 1 : 0;
        
        //System.out.println("within put key=" + key + " cmp=" + cmp);
        
        if (cmp < 0) {
            h.left  = put(h.left,  key, val);
        } else if (cmp > 0) {
            h.right = put(h.right, key, val);
        } else {
            h.val   = val;
        }
        //NOTE: sizes are not updated yet
        
        // fix-up any right-leaning links
        if (isRed(h.right) && !isRed(h.left)) {
            h = rotateLeft(h);
        }
        if (isRed(h.left)  &&  isRed(h.left.left)) {
            h = rotateRight(h);
        }
        if (isRed(h.left)  &&  isRed(h.right)) {
            flipColors(h);
        }
        h.size = size(h.left) + size(h.right) + 1;

        //printPreOrderTraversal();
        //System.out.println("after put h=" + h + " key=" + key + " root=" 
        //    + root.toString(findParent(root)));
        
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
    private Node deleteMin(Node h) { 
        
        //System.out.println("deleteMin(" + h.toString() + ")");
        
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

        //System.out.println("deleteMax.  root=" + root);
        //printPreOrderTraversal();
        
        // if both children of root are black, set root to red
        if (!isRed(root.left) && !isRed(root.right))
            root.color = RED;

        root = deleteMax(root);
        if (!isEmpty()) root.color = BLACK;
        // assert check();
    }

    // delete the key-value pair with the maximum key rooted at h
    private Node deleteMax(Node h) { 
    
        //TODO: make this iterative
        if (isRed(h.left)) {
            h = rotateRight(h);
        }

        if (h.right == null) {
            // h is max key
            return null;
        }

        if (!isRed(h.right) && !isRed(h.right.left)) {
            h = moveRedRight(h);
        }

        h.right = deleteMax(h.right);

        return balance(h);
    }

    /**
     * Removes the specified key and its associated value from this symbol table     
     * (if the key is in this symbol table).    
     *
     @param  key the key
     */
    public void delete(long key) { 
        
        //System.out.println("delete " + key);
        
        if (!contains(key)) return;

        // if both children of root are black, set root to red
        if (!isRed(root.left) && !isRed(root.right))
            root.color = RED;

        root = delete(root, key);
        if (!isEmpty()) root.color = BLACK;
        // assert check();
        
        //System.out.println("after delete " + key + " root=" + root);
    }

    // delete the key-value pair with the given key rooted at h
    private Node delete(Node h, long key) { 
        // assert get(h, key) != null;

        if (key < h.key)  {
            if (!isRed(h.left) && !isRed(h.left.left)) {
                h = moveRedLeft(h);
            }
            h.left = delete(h.left, key);
        } else {
            if (isRed(h.left)) {
                h = rotateRight(h);
            }
            if (key == h.key && (h.right == null)) {
                return null;
            }
            if (!isRed(h.right) && !isRed(h.right.left)) {
                h = moveRedRight(h);
            }
            //System.out.println("key=" + key + " h.key=" + h.key);
            if (key == h.key) {
                
                Node x = min(h.right);
                
                //System.out.println(
                //    "   x to get h fields except val. " + 
                //    "\n   x=" + x.toString(findParent(x))
                //    + "\n   h=" + h.toString(findParent(h))
                //);
                
                h.key = x.key;
                h.val = x.val;
                // h.val = get(h.right, min(h.right).key);
                // h.key = min(h.right).key;
                                
                h.right = deleteMin(h.right);
                
                //System.out.println("   after deleteMin x =" + 
                //    h.toString(findParent(h)));
              
            } else {
                h.right = delete(h.right, key);
            }
        }
        
        return balance(h);
    }

   /***************************************************************************
    *  Red-black tree helper functions.
     @param h
     @return 
    ***************************************************************************/

    // make a left-leaning link lean to the right
    private Node rotateRight(Node h) {
        // assert (h != null) && isRed(h.left);
        Node x = h.left;
      
        h.left = x.right;
        x.right = h;
        x.color = x.right.color;
        x.right.color = RED;
        x.size = h.size;
        h.size = size(h.left) + size(h.right) + 1;
            
        return x;
    }

    // make a right-leaning link lean to the left
    private Node rotateLeft(Node h) {
                
        // assert (h != null) && isRed(h.right);
        Node x = h.right;
        h.right = x.left;
        x.left = h;
        x.color = x.left.color;
        x.left.color = RED;
        x.size = h.size;
        h.size = size(h.left) + size(h.right) + 1;
                
        return x;
    }

    // flip the colors of a node and its two children
    private void flipColors(Node h) {
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
    private Node moveRedLeft(Node h) {
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
    private Node moveRedRight(Node h) {
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
    private Node balance(Node h) {
        // assert (h != null);

        if (isRed(h.right))   {
            h = rotateLeft(h);
        }
        if (isRed(h.left) && isRed(h.left.left)) {
            h = rotateRight(h);
        }
        if (isRed(h.left) && isRed(h.right))  {
            flipColors(h);
        }

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
    private int height(Node x) {
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
    public long min() {
        if (isEmpty()) throw new NoSuchElementException("called min() with empty symbol table");
        return min(root).key;
    } 

    // the smallest key in subtree rooted at x; null if no such key
    private Node min(Node x) { 
        assert(x != null);
        while (x.left != null) {
            x = x.left;
        }
        return x;
    } 

    /**
     * Returns the largest key in the symbol table.
     @return the largest key in the symbol table
     * @throws NoSuchElementException if the symbol table is empty
     */
    public long max() {
        if (isEmpty()) throw new NoSuchElementException("called max() with empty symbol table");
        return max(root).key;
    } 

    // the largest key in the subtree rooted at x; null if no such key
    private Node max(Node x) { 
        assert(x != null);
        while (x.right != null) {
            x = x.right;
        }
        return x; 
    }

    /**
     * Returns the largest key in the symbol table less than or equal to {@code key}.
     @param key the key
     @param output if output[0] == -1, the key was not present, 
     * else output[1] holds the largest key in the symbol table less than or equal to {@code key}
     @throws NoSuchElementException if the tree is empty
     */
    public void floor(long key, long[] output) {
        if (output == null || output.length != 2) {
            throw new IllegalArgumentException("output must be length 2");
        }
        if (isEmpty()) {
            output[0] = -1;
            throw new NoSuchElementException("called floor() with empty symbol table");
        }
        Node x = floor(root, key);
        if (x == null) {
            output[0] = -1;
        } else {
            output[0] = 0;
            output[1] = x.key;
        }
    }    

    /**
     * Returns the largest key in the symbol table less than {@code key}.
     @param key the key
     @param output if output[0] == -1, the key was not present, 
     * else output[1] holds the largest key in the symbol table less than or equal to {@code key}
     @throws NoSuchElementException if the tree is empty
     */
    public void lower(long key, long[] output) {
        if (output == null || output.length != 2) {
            throw new IllegalArgumentException("output must be length 2");
        }
        if (isEmpty()) {
            output[0] = -1;
            throw new NoSuchElementException("called floor() with empty symbol table");
        }
        List<Node> stack = new ArrayList<Node>();
        Node x = lower(root, key, stack);
        if (x == null) {
            output[0] = -1;
        } else {
            output[0] = 0;
            output[1] = x.key;
        }
    }    
    
    // the largest key in the subtree rooted at x less than or equal to the given key
    private Node floor(Node x, long key) {
        if (x == null) return null;
        int cmp = (key < x.key) ? -1 : (key > x.key) ? 1 : 0;
        if (cmp == 0) {
            return x;
        }
        if (cmp < 0)  {
            return floor(x.left, key);
        }
        Node t = floor(x.right, key);
        if (t != null) {
            return t;
        } else {
            return x;
        }
    }
    
    /** the largest key in the subtree rooted at x less than the given key
     * 
     * The method is symmetric t the successor method called 
     * higher, as suggested by Cormen, Leiserson, Rivest, and Stein in
     * the book "Introduction to Algorithms".
     * 
     @param x
     @param key
     @param stack
     @return 
     */
    private Node lower(Node x, long key, List<Node> stack) { 
        
        //binary search until overshoot
        while (x != null && key != x.key) {
            stack.add(x);
            if (key < x.key) {
                x = x.left;
            } else {
                x = x.right;
            }
        }
        
        Node y = null;
        int yIdx = -1;
        
        if (x == null) {
            if (stack.size() > 1) {
                x = stack.get(stack.size() - 1);
                yIdx = stack.size() - 2;
                y = stack.get(yIdx);
            } else if (stack.size() == 1) {
                x = stack.get(stack.size() - 1);
            } else {
                //x and y are null and stack was empty so method argument x was also null
                return null;
            }
        } else if (!stack.isEmpty()) {
            yIdx = stack.size() - 1;
            y = stack.get(yIdx);
        }
       
        if (x.left != null) {
            return max(x.left);
        }
      
        while (y != null && x == y.left) {
            //System.out.println("lower: y=" + y.key + " q=" + key);
            x = y;
            yIdx--;
            if (yIdx < 0) break;
            y = stack.get(yIdx);
        }
        
        if (x.key < key && x.key > y.key) {
            return x;
        }
        
        //System.out.println("    y=" + y.key + " q=" + key);
        if (y != null && y.key >= key) {
            return null;
        }
        
        return y;
    }

    /**
     * Returns the smallest key in the symbol table greater than or equal to {@code key}.
     @param key the key
     @param output if output[0] == -1, the key was not present, 
     * else output[1] holds
     * the smallest key in the symbol table greater than or equal to {@code key}
     * @throws NoSuchElementException if the tree is empty
     */
    public void ceiling(long key, long[] output) {
        if (output == null || output.length != 2) {
            throw new IllegalArgumentException("output must be length 2");
        }
        if (isEmpty()) {
            output[0] = -1;
            throw new NoSuchElementException("called floor() with empty symbol table");
        }
        Node x = ceiling(root, key);
        if (x == null) {
            output[0] = -1;
        } else   {
            output[0] = 0;
            output[1] = x.key;
        }
    }

    /**
     * Returns the smallest key in the symbol table greater than {@code key}.
     @param key the key
     @param output if output[0] == -1, the key was not present, 
     * else output[1] holds
     * the smallest key in the symbol table greater than or equal to {@code key}
     * @throws NoSuchElementException if the tree is empty
     */
    public void higher(long key, long[] output) {
        if (output == null || output.length != 2) {
            throw new IllegalArgumentException("output must be length 2");
        }
        if (isEmpty()) {
            output[0] = -1;
            throw new NoSuchElementException(
                "called higher() with empty symbol table");
        }
        List<Node> stack = new ArrayList<Node>();
        Node x = higher(root, key, stack);
        if (x == null) {
            output[0] = -1;
        } else   {
            output[0] = 0;
            output[1] = x.key;
        }
    }
    
    // the smallest key in the subtree rooted at x greater than or equal to the given key
    private Node ceiling(Node x, long key) {
        //TODO: make this iterative
        if (x == null) return null;
        int cmp = key < x.key ? -1 : (key > x.key) ? 1 : 0; 
        if (cmp == 0) return x;
        if (cmp > 0)  return ceiling(x.right, key);
        Node t = ceiling(x.left, key);
        if (t != null) return t; 
        else           return x;
    }
    
    /** the smallest key in the subtree rooted at x greater than the given key.
     * 
     * NOTE: the method uses in part, a pattern adapted from the Cormen, Leiserson, Rivest, and Stein
     * book "Introduction to Algorithms" for their Red Black Tree.
     * 
     @param x
     @param key
     @param stack
     @return 
     */
    private Node higher(Node x, long key, List<Node> stack) {  
                        
        /*
                    X
        left .lte.     right .gte.
        */
        
        //binary search until overshoot
        while (x != null && key != x.key) {
            stack.add(x);
            //System.out.println("higher: x=" + x.key + " q=" + key);
            if (key < x.key) {
                x = x.left;
            } else {
                //right is towards numbers larger than x.key
                x = x.right;
            }
        }
        
        Node y = null;
        int yIdx = -1;
        
        if (x == null) {
            if (stack.size() > 1) {
                x = stack.get(stack.size() - 1);
                yIdx = stack.size() - 2;
                y = stack.get(yIdx);
            } else if (stack.size() == 1) {
                x = stack.get(stack.size() - 1);
            } else {
                //x and y are null and stack was empty so method argument x was also null
                return null;
            }
        } else if (!stack.isEmpty()) {
            yIdx = stack.size() - 1;
            y = stack.get(yIdx);
        }
        
        if (x.right != null) {
            return min(x.right);
        }
        
        while (y != null && x == y.right) {
            //System.out.println("higher: y=" + y.key + " q=" + key);
            x = y;
            yIdx--;
            if (yIdx < 0) break;
            y = stack.get(yIdx);
        }
        
        if (x.key > key && x.key < y.key) {
            return x;
        }
        
        if (y != null && y.key <= key) {
            return null;
        }
        
        return y;
    }

    /**
     * Return the kth smallest key in the symbol table.
     @param k the order statistic
     @return the {@code k}th smallest key in the symbol table
     * @throws IllegalArgumentException unless {@code k} is between 0 and
     *     <em>n</em>–1
     */
    public long select(int k) {
        if (k < 0 || k >= size()) {
            throw new IllegalArgumentException("called select() with invalid argument: " + k);
        }
        Node x = select(root, k);
        return x.key;
    }

    // the key of rank k in the subtree rooted at x
    private Node select(Node x, int k) {
        //TODO: make this iterative
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
    public int rank(long key) {
        return rank(key, root);
    } 

    // number of keys less than key in the subtree rooted at x
    private int rank(long key, Node x) {
        //TODO: make this iterative
        if (x == null) return 0; 
        int cmp = key < x.key ? -1 : (key > x.key) ? 1 : 0;  
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
    public TLongList keys() {
        if (isEmpty()) return new TLongArrayList();
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
    public TLongList keys(long lo, long hi) {

        TLongList queue = new TLongArrayList();
        // if (isEmpty() || lo.compareTo(hi) > 0) return queue;
        keys(root, queue, lo, hi);
        return queue;
    } 

    // add the keys between lo and hi in the subtree rooted at x
    // to the queue
    private void keys(Node x, TLongList queue, long lo, long hi) { 
        //TODO: make this iterative
        if (x == null) return; 
        int cmplo = lo < x.key ? -1 : (lo > x.key) ? 1 : 0;
        int cmphi = hi < x.key ? -1 : (hi > x.key) ? 1 : 0; 
        if (cmplo < 0) keys(x.left, queue, lo, hi); 
        if (cmplo <= 0 && cmphi >= 0) queue.add(x.key); 
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
    public int size(long lo, long hi) {

        if (lo > hi) return 0;
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
    private boolean isBST(Node x, Long min, Long max) {
        if (x == null) return true;
        if (min != null && x.key <= min) return false;
        if (max != null && x.key >= max) return false;
        return isBST(x.left, min, Long.valueOf(x.key)) 
            && isBST(x.right, Long.valueOf(x.key), max);
    } 

    // are the size fields correct?
    private boolean isSizeConsistent() { return isSizeConsistent(root); }
    private boolean isSizeConsistent(Node x) {
        if (x == null) return true;
        if (x.size != size(x.left) + size(x.right) + 1) return false;
        return isSizeConsistent(x.left) && isSizeConsistent(x.right);
    } 

    // check that ranks are consistent
    private boolean isRankConsistent() {
        for (int i = 0; i < size(); i++) {
            if (i != rank(select(i))) return false;
        }
        TLongList keys = keys();
        for (int i = 0; i < keys.size(); ++i) {
            long key = keys.get(i);
            if (key != select(rank(key))) {
                return false;
            }
        }
        return true;
    }

    // Does the tree have no red right links, and at most one (left)
    // red links in a row on any path?
    private boolean is23() { return is23(root); }
    private boolean is23(Node x) {
        if (x == null) return true;
        if (isRed(x.right)) return false;
        if (x != root && isRed(x) && isRed(x.left))
            return false;
        return is23(x.left) && is23(x.right);
    } 

    // do all paths from root to leaf have same number of black edges?
    private boolean isBalanced() { 
        int black = 0;     // number of black links on path from root to min
        Node x = root;
        while (x != null) {
            if (!isRed(x)) black++;
            x = x.left;
        }
        return isBalanced(root, black);
    }

    // does every path from the root to a leaf have the given number of black links?
    private boolean isBalanced(Node x, int black) {
        if (x == null) return black == 0;
        if (!isRed(x)) black--;
        return isBalanced(x.left, black) && isBalanced(x.right, black);
    } 
    
    /**
     * left subtree, root, right subtree
     */
    public void printInOrderTraversal() {
        Node[] nodes = getInOrderTraversalIterative(root);
        for (Node node : nodes) {
            System.out.println("node=" + node.toString());
        }
    }
    
    /**
     * root, left subtree, right subtree
     */
    public void printPreOrderTraversal() {
        Node[] nodes = getPreOrderTraversalIterative(root, 0);
        for (Node node : nodes) {
            Node p = findParent(node, nodes);
            System.out.println("node=" + node.toString(p));
        }
    }

    /**
     *
     @param topNode
     */
    public void printPreOrderTraversal2(Node topNode) {
        Node[] nodes = getPreOrderTraversalIterative(topNode, 0);
        for (Node node : nodes) {
            Node p = findParent(node, nodes);
            System.out.println("node=" + node.toString(p));
        }
    }
    private void printPreOrderTraversal(int addExtraToSize) {
        Node[] nodes = getPreOrderTraversalIterative(root, addExtraToSize);
        for (Node node : nodes) {
            if (node == null) {
                continue;
            }
            Node p = findParent(node, nodes);
            System.out.println("node=" + node.toString(p));
        }
    }
    
    private Node findParent(Node h) {
        Node[] nodes = getPreOrderTraversalIterative(root, 1);
        return findParent(h, nodes);
    }
    private Node findParent(Node h, Node[] nodes) {
        
        for (Node node : nodes) {
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
        Node[] nodes = getPostOrderTraversalIterative(root);
        for (Node node : nodes) {
            System.out.println("node=" + node.toString());
        }
    }

    /**
     * visit each node using pattern left subtree, root, right subtree
     * in an iterative manner rather than invoking the method recursively.
     @param node
     @return 
     */
    protected Node[] getInOrderTraversalIterative(Node node) {
        
        if (isEmpty()) {
            return new Node[0];
        }
        
        int sz = size();
       
        Node[] array = new Node[sz];
        int count = 0;
        
        Stack<Node> stack = new Stack<>();
               
        while (!stack.isEmpty() || (node != null)) {
            if (node != null) {
                 
                stack.push(node);
                
                node = node.left;
            
            } else if (count < sz) {
                
                node = stack.pop();
                
                array[count] = node;
                count++;
                
                //System.out.println(node.key);
                
                node = node.right;
            }
        }
        return array;
    }
    
    /**
     * visit each node using pattern: root, left subtree, right subtree
     * in an iterative manner rather than invoking the method recursively.
     @param node
     @param addExtraToSize
     @return 
     */
    protected Node[] getPreOrderTraversalIterative(Node node, int addExtraToSize) {
       
        //NOTE: added additional integer and conditions 
        //   for size because may be printing tree
        //   in the middle of a put where the node size is not yet updated.
        // The count conditionals below are otherwise, not needed.
        
        if (isEmpty()) {
            return new Node[0];
        }
                
        int sz = size() + addExtraToSize;
        
        Node[] array = new Node[sz];
        int count = 0;
        
        Stack<Node> stack = new Stack<>();
        
        while (count < sz && (!stack.isEmpty() || node != null)) {
            if (node != null && count < sz) {
              
                array[count] = node;
                
                count++;
                
                if (count < sz) {
                    stack.push(node);
                }
                
                node = node.left;
            
            } else if (count < sz) {
                
                node = stack.pop();
                
                node = node.right;
            }
        }
        
        if (count < array.length) {
            array = Arrays.copyOf(array, count);
        }
        
        return array;
    }

    /**
     * visit each node using pattern: left subtree, right subtree, root subtree
     * in an iterative manner rather than invoking the method recursively.
     @param node
     @return 
     */
    protected Node[] getPostOrderTraversalIterative(Node node) {
        
        if (isEmpty()) {
            return new Node[0];
        }
        
        Node[] array = new Node[size()];
        int count = 0;
        
        if (node == null) {
            return array;
        }
        
        Stack<Node> stack = new Stack<Node>();
        Stack<Node> stack2 = new Stack<Node>();
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
            array[count] = node;
            count++;
            //System.out.println(node.key);
        }
         
        return array;
    }
    
    /**
     *
     @return
     */
    public static long estimateNodeSizeOnHeap() {
        return Node.estimateSizeOnHeap();
    }
    
    /**
     * estimate the size that an instance of RedBlackBSTLongInt with
     * n entries would occupy in heap space in Bytes.
     * 
     @param numberOfEntries amount of space for this object's instance
     * with n entries in Bytes on the heap.
     * 
     @return 
     */
    public static long estimateSizeOnHeap(int numberOfEntries) {
        
        long total = 0;
        
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNBooleanFields(2);
        est.setNObjRefsFields(1);
       
        total += est.estimateSizeOnHeap();
        
        total += numberOfEntries * Node.estimateSizeOnHeap();
    
        return total;
    }
}
