package thirdparty.edu.princeton.cs.algs4;

/******************************************************************************
   adapted from RedBlackBST.java 
   from the 
   book "Algorithms" by Sedgewick and Wayne
   http://algs4.cs.princeton.edu/33balanced/RedBlackBST.java
   copyright for authors Robert Sedgewick and Kevin Wayne
   is GPLV3, http://algs4.cs.princeton.edu/faq/

This version uses smaller amount of memory by replacing linked nodes with
associative arrays.

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
import gnu.trove.map.TLongIntMap;
import gnu.trove.map.TLongLongMap;
import gnu.trove.map.hash.TLongIntHashMap;
import gnu.trove.map.hash.TLongLongHashMap;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.NoSuchElementException;
import java.util.Stack;

/**
 * NOTE: this implementation of the long-int Red Black tree uses the least
 * amount of memory of all versions by replacing the tree structure of linked
 * lists of nodes as objects with associative arrays that use primitive 
 * arrays internally.
 * 
 * 
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
 *  For additional documentation, see 
 * <a href="http://algs4.cs.princeton.edu/33balanced">Section 3.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *  For other implementations of the same API, see {@link ST}, {@link BinarySearchST},
 *  {@link SequentialSearchST}, {@link BST},
 *  {@link SeparateChainingHashST}, {@link LinearProbingHashST}, and {@link AVLTreeST}.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 
 * edits made to their original code include replacing the linked object Nodes
 * with an associative array holding multiple values.
 */
public class RedBlackBSTLongInt2 {
    
    private static final int RED   = 1;
    private static final int BLACK = 0;

    protected long root = -1;
    protected boolean rootIsSet = false;

    //TODO: as soon as this is debugged,
    //   make a class that extends TLongLongMap for
    //   key long, and values long, int, long, long, int, int
    //   reducing the number of long keys from 6 to 1
    protected final TLongLongMap keyParentMap;
    protected final TLongIntMap keyValMap;
    protected final TLongLongMap keyLeftMap;
    protected final TLongLongMap keyRightMap;
    protected final TLongIntMap keyColorMap;
    protected final TLongIntMap keySizeMap;
   
    /**
     * Initializes an empty symbol table.
     */
    public RedBlackBSTLongInt2() {
        keyParentMap = new TLongLongHashMap();
        keyValMap = new TLongIntHashMap();
        keyLeftMap = new TLongLongHashMap();
        keyRightMap = new TLongLongHashMap();
        keyColorMap = new TLongIntHashMap();
        keySizeMap = new TLongIntHashMap();
    }
    
    /**
     * Initializes an empty symbol table.
     */
    public RedBlackBSTLongInt2(int capacity) {
        keyParentMap = new TLongLongHashMap(capacity);
        keyValMap = new TLongIntHashMap(capacity);
        keyLeftMap = new TLongLongHashMap(capacity);
        keyRightMap = new TLongLongHashMap(capacity);
        keyColorMap = new TLongIntHashMap(capacity);
        keySizeMap = new TLongIntHashMap(capacity);
    }

    private long addNewNode(long key, int val, int color, int size) {
        //TODO: consider adding assert that key doesn't exist
        keyValMap.put(key, val);
        keyColorMap.put(key, color);
        keySizeMap.put(key, size);
        return key;
    }
    
   /***************************************************************************
    *  Node helper methods.
    ***************************************************************************/
    // is node x red; false if x is null ?
    private boolean isRed(long x) {
        if (!keyColorMap.containsKey(x)) return false;
        return keyColorMap.get(x) == RED;
    }
    private boolean isLeftRed(long x) {
        if (!keyLeftMap.containsKey(x)) {
            return false;
        }
        return isRed(keyLeftMap.get(x));
    }
    private boolean isRightRed(long x) {
        if (!keyRightMap.containsKey(x)) {
            return false;
        }
        return isRed(keyRightMap.get(x));
    }
    private boolean isLeftLeftRed(long x) {
        if (!keyLeftMap.containsKey(x)) {
            return false;
        }
        long left = keyLeftMap.get(x);
        return isLeftRed(left);
    }
    private boolean isRightLeftRed(long x) {
        if (!keyRightMap.containsKey(x)) {
            return false;
        }
        long right = keyRightMap.get(x);
        return isLeftRed(right);
    }

    // number of node in subtree rooted at x; 0 if x is null
    private int size(long x) {
        if (!keySizeMap.containsKey(x)) return 0;
        return keySizeMap.get(x);
    }
    private int sizeLeft(long x) {
        if (!keySizeMap.containsKey(x)) return 0;
        if (!keyLeftMap.containsKey(x)) return 0;
        return keySizeMap.get(keyLeftMap.get(x));
    }
    private int sizeRight(long x) {
        if (!keySizeMap.containsKey(x)) return 0;
        if (!keyRightMap.containsKey(x)) return 0;
        return keySizeMap.get(keyRightMap.get(x));
    }

    /**
     * Returns the number of key-value pairs in this symbol table.
     * @return the number of key-value pairs in this symbol table
     */
    public int size() {
        return size(root);
    }

   /**
     * Is this symbol table empty?
     * @return {@code true} if this symbol table is empty and {@code false} otherwise
     */
    public boolean isEmpty() {
        return !rootIsSet;
    }

   /***************************************************************************
    *  Standard BST search.
    ***************************************************************************/

    /**
     * Returns the value associated with the given key.
     * @param key the key
     * @param output if output[0] == -1, then key was not present, else the
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
    private void get(long x, long key, int[] output) {
        //while (x != null) {
        while (keyValMap.containsKey(x)) {
            int cmp = key < x ? -1 : (key > x) ? 1 : 0; 
            if (cmp < 0) {
                if (keyLeftMap.containsKey(x)) {
                    x = keyLeftMap.get(x);
                } else {
                    break;
                }
            } else if (cmp > 0) {
                if (keyRightMap.containsKey(x)) {
                    x = keyRightMap.get(x);
                } else {
                    break;
                }
            } else {
                output[1] = keyValMap.get(x);
                return;
            }
        }
        output[0] = -1;
    }

    /**
     * Does this symbol table contain the given key?
     * @param key the key
     * @return {@code true} if this symbol table contains {@code key} and
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
     * @param key the key
     * @param val the value
     */
    public void put(long key, int val) {
        
        //System.out.println("put " + key + ":");
        
        //System.out.println("before put " + key);
        //printPreOrderTraversal();
        
        root = put(root, key, val);
        rootIsSet = true;
        keyColorMap.put(root, BLACK);
        
        if (keyParentMap.containsKey(root)) {
            keyParentMap.remove(root);
        }
        
        //printPreOrderTraversal();
        //System.out.println("after put " + key + " root=" + root);
        
        assert(check());
    }

    // insert the key-value pair in the subtree rooted at h
    //private Node put(Node h, long key, int val) {
    private long put(long h, long key, int val) {
        
        //System.out.println("put h=" + h + " key=" + key);
        
        if (!rootIsSet || !keyValMap.containsKey(h)) {
            return addNewNode(key, val, RED, 1);
        }
       
        int cmp = (key < h) ? -1 : (key > h) ? 1 : 0;
        if (cmp < 0) {
            //h.left  = put(h.left,  key, val);
            long putKey;
            if (keyLeftMap.containsKey(h)) {
                //h.left = putKey
                //putKey.parent = h
                long hLeft = keyLeftMap.get(h);
                putKey = put(hLeft,  key, val);
            } else {
                putKey = addNewNode(key, val, RED, 1);
            }
            keyLeftMap.put(h, putKey);
            keyParentMap.put(putKey, h);
        } else if (cmp > 0) {
            //h.right = put(h.right, key, val);
            long putKey;
            if (keyRightMap.containsKey(h)) {
                long hRight = keyRightMap.get(h);
                putKey = put(keyRightMap.get(h),  key, val);
            } else {
                putKey = addNewNode(key, val, RED, 1);
            }
            keyRightMap.put(h, putKey);
            keyParentMap.put(putKey, h);
        } else {
            //h.val   = val;
            keyValMap.put(h, val);
        }
        
        // fix-up any right-leaning links
        if (isRightRed(h) && !isLeftRed(h)) {
            h = rotateLeft(h);
        }
        if (isLeftRed(h)  &&  isLeftLeftRed(h)) {
            h = rotateRight(h);
        }
        if (isLeftRed(h)  &&  isRightRed(h))  {
            flipColors(h);
        }
        int size = sizeLeft(h) + sizeRight(h) + 1;
        keySizeMap.put(h, size);
       
        return h;
    }

   /***************************************************************************
    *  Red-black tree deletion.
    ***************************************************************************/

    private void setRootToRedIfChildrenAreBlack() {
        if (!isLeftRed(root) && !isRightRed(root)) {
            keyColorMap.put(root, RED);
        }
    }
    
    /**
     * Removes the smallest key and associated value from the symbol table.
     * @throws NoSuchElementException if the symbol table is empty
     */
    public void deleteMin() {
        if (isEmpty()) throw new NoSuchElementException("BST underflow");

        //System.out.println("deleteMin.  root=" +root);
        //printPreOrderTraversal();
        
        int sz0 = size();
        
        assert(keyValMap.containsKey(root));
        
        // if both children of root are black, set root to red
        setRootToRedIfChildrenAreBlack();
        
        //root = deleteMin(root);
        long[] output = new long[2];
        deleteMin(root, output);
        
        //System.out.println("deleteMin=" + Arrays.toString(output));
        
        if (output[0] == -1) {
            rootIsSet = false;
            deleteFromMaps(root);
            root = -1;
        } else {
            System.out.println("in deleteMin: assigning root=" + output[1]);
            root = output[1];
            keyParentMap.remove(root);
        }
        if (!isEmpty()) {
            //root.color = BLACK;
            keyColorMap.put(root, BLACK);
        }
        
        //System.out.format("AFTER deleteMin()\n");
        //printPreOrderTraversal(1);
        
        assert(check());
        
        assert(sz0 == (size() + 1));
    }

    // delete the key-value pair with the minimum key rooted at h.
    // output is length==2, output[0] holds the error code
    // and is -1 when there is an error.  output[1] holds the top node
    // from which the min was deleted (note that the top node may have
    // been rotated, so might not equal h).
    // if output.length == 3, the third item is set to be the
    // minimum node which was deleted.
    private void deleteMin(long h, long[] output) { 
        
        //System.out.println("deleteMin " + nodeToString(h));
        
        if (!keyValMap.containsKey(h)) {
            output[0] = -1;
            return;
        } else if (!keyLeftMap.containsKey(h)) {
            output[0] = -1;
            if (output.length == 3) {
                output[2] = h;
            }
            return;
        }

        //System.out.format("BEFORE deleteMin(%d)\n", h);
        //printPreOrderTraversal(1);
        
        
        if (!isLeftRed(h) && !isLeftLeftRed(h)) {
            h = moveRedLeft(h);
        }
        
        //h.left = deleteMin(h.left);
        
        if (!keyLeftMap.containsKey(h)) {
            output[0] = -1;
            return;
        }
        long left = keyLeftMap.get(h);
        //assert(parent == h);
        
        if (output.length == 3) {
            output[2] = left;
        }
        
        deleteMin(left, output);
        if (output[0] == -1) {
            //h.left = null
            keyLeftMap.remove(h);
            output[0] = 0;
        } else {
            System.out.println(" deleteMin return is " + output[1] 
               + " h.left gets assigned it");
            keyParentMap.put(output[1], h);
            keyLeftMap.put(h, output[1]);
        }
        output[1] = balance(h);
        
        //printPreOrderTraversal(1);
        //System.out.format("AFTER deleteMin(h)\n");
    }

    /**
     * Removes the largest key and associated value from the symbol table.
     * @throws NoSuchElementException if the symbol table is empty
     */
    public void deleteMax() {
        if (isEmpty()) throw new NoSuchElementException("BST underflow");

        //System.out.println("deleteMax.  root=" +root);
        //printPreOrderTraversal();
        
        // if both children of root are black, set root to red
        setRootToRedIfChildrenAreBlack();

        //root = deleteMax(root);
        
        long[] output = new long[2];
        deleteMax(root, output);
        if (output[0] == -1) {
            deleteFromMaps(root);
            rootIsSet = false;
            root = -1;
            return;
        }
        root = output[1];
        keyParentMap.remove(root);
        
        if (!isEmpty()) {
            //root.color = BLACK;
            keyColorMap.put(root, BLACK);
        }
        
        //System.out.format("AFTER deleteMax()\n");
        //printPreOrderTraversal(1);
        
        assert(check());
    }
    
    private void deleteFromMaps(long key) {
        /*
            p0--> key->left
                     ->right
        */
        
        keyValMap.remove(key);
        if (keyLeftMap.containsKey(key)) {
            //keyParentMap.remove(keyLeftMap.get(key));
            keyLeftMap.remove(key);
        }
        if (keyRightMap.containsKey(key)) {
            //keyParentMap.remove(keyRightMap.get(key));
            keyRightMap.remove(key);
        }
        keyColorMap.remove(key);
        keySizeMap.remove(key);
        
    }

    // delete the key-value pair with the maximum key rooted at h
    private void deleteMax(long h, long[] output) { 
    
        //System.out.format("BEFORE deleteMax(%d)\n", h);
        //printPreOrderTraversal(1);
        
        if (isLeftRed(h)) {
            h = rotateRight(h);
        }
      
        if (!keyRightMap.containsKey(h)) {
            // h is max key
            output[0] = -1;
            return;
        }

        if (!isRightRed(h) && !isRightLeftRed(h)) {
            //move red node down the right spine of the tree
            h = moveRedRight(h);
        }
        
        //h.right = deleteMax(h.right);
        
        assert(keyRightMap.containsKey(h));
        
        deleteMax(keyRightMap.get(h), output);
        if (output[0] == -1) {
            keyRightMap.remove(h);
        } else {
            keyParentMap.put(output[1], h);
            keyRightMap.put(h, output[1]);
        }

        output[0] = 0;
        output[1] = balance(h);
        
        
        //System.out.format("AFTER deleteMax(h)\n");
        //printPreOrderTraversal(1);
        
    }

    /**
     * Removes the specified key and its associated value from this symbol table     
     * (if the key is in this symbol table).    
     *
     * @param  key the key
     */
    public void delete(long key) { 
        
        int sz0 = size();
        
        System.out.println("\nbefore delete " + key + " root=" + root);
        printPreOrderTraversal();
        
        if (!contains(key)) return;

        // if both children of root are black, set root to red
        setRootToRedIfChildrenAreBlack();
        
        ////root = delete(root, key);
        
        long[] output = new long[2];
        delete(root, key, output);
        
        if (output[0] == -1) {
            rootIsSet = false;
            deleteFromMaps(root);
            root = -1;
            return;
        }
        root = output[1];
        keyParentMap.remove(root);
        
        if (!isEmpty()) {
            keyColorMap.put(root, BLACK);
        }
        
        printPreOrderTraversal();
        System.out.println("after delete " + key + " root=" + root);
        
        assert(check());
        
        assert(sz0 == (size() + 1));
    }

    // delete the key-value pair with the given key rooted at h
    private void delete(long h, long key, long[] output) { 
        
        System.out.format("delete(%d, %d)\n", h, key);
        
        assert(keyValMap.containsKey(key));
        assert(keyValMap.containsKey(h));
        
        // assert get(h, key) != null;
        {//DEBUG
            int[] vOutput = new int[2];
            get(h, key, vOutput);
            assert(vOutput[0] != -1);
        }
        
        output[0] = 0;
        
        //NOTE: this method and the methods it uses do not always
        //check for nulls and handle them
                    
        if (key < h)  {
            if (!isLeftRed(h) && !isLeftLeftRed(h)) {
                h = moveRedLeft(h);
            }
            assert(keyLeftMap.containsKey(h));
            //h.left = delete(h.left, key);
            deleteFromLeftAssignLeft(h, h, key, output);
            output[0] = 0;
        } else {
            if (isLeftRed(h)) {
                h = rotateRight(h);
            }
            if (key == h && !keyRightMap.containsKey(h)) {
                output[0] = -1;
                return;
            }
            if (!isRightRed(h) && !isRightLeftRed(h)) {
                h = moveRedRight(h);
            }
            if (key == h) {
               
                /*
                useful in visualizing this case is
                https://www.cs.princeton.edu/~rs/talks/LLRB/RedBlack.pdf
                pg 66, "Deleting an arbitrary node"
                
                To delete 'D':
                
                               H
                        D             L
                    B      F       J     N
                   A C    E G     I K   M
                
                  D.key = min(D.right)
                  D.value = get(D.right, D.key)
                  D.right = deleteMin(D.right)
                  then delete the min node
                  then fix right-leaning red link of F, the parent
                    of the just deleted E node
                              H
                        E             L
                    B      F       J     N
                   A C      G     I K   M   
                
                
                              H
                        E             L
                    B      F       J     N
                   A C    G       I K   M   
                */
                
                //Node x = min(h.right);
                // Note: need to use deleteMin all in one because there
                //   are rotations in it that may change the branch
                //   traversed and hence the minimum of that branch
                //   might not be the same as x.
                long[] output3 = new long[3];
                deleteMin(keyRightMap.get(h), output3);
                        
                System.out.println("deleteMin=" + Arrays.toString(output3));
        
                int hClr = keyColorMap.get(h);
                boolean hLeftExists = keyLeftMap.containsKey(h);
                boolean hRightExists = keyRightMap.containsKey(h);
                boolean hParentExists = keyParentMap.containsKey(h);
                long hLeft = keyLeftMap.get(h);
                long hRight = keyRightMap.get(h);
                long hParent = keyParentMap.get(h);
                int hSize = keySizeMap.get(h);
               
                String hDBG = nodeToString(h);
                
                long x = output3[2];
               
                System.out.println(
                    "   x to get h fields except val. " + 
                    "\n   x=" + nodeToString(x)
                    + "\n   h=" + hDBG
                );
               
                assert(x != h);
                
                //int xClr = keyColorMap.get(x);
                //boolean xLeftExists = keyLeftMap.containsKey(x);
                //boolean xRightExists = keyRightMap.containsKey(x);
                boolean xParentExists = keyParentMap.containsKey(x);
                //long xLeft = keyLeftMap.get(x);
                //long xRight = keyRightMap.get(x);
                long xParent = keyParentMap.get(x);
                int xVal = keyValMap.get(x);
                
                // finish deleting x if not already
                if (keyParentMap.containsKey(x)) {
                    long minParent = keyParentMap.get(x);
                    if (keyLeftMap.containsKey(minParent)
                        && keyLeftMap.get(minParent) == x) {
                        keyLeftMap.remove(minParent);
                    } else if (keyRightMap.containsKey(minParent)
                            && keyRightMap.get(minParent) == x) {
                        keyRightMap.remove(minParent);
                    }
                    keyParentMap.remove(x);
                }
                if (keyValMap.containsKey(x)) {
                    keyValMap.remove(x);
                    keyColorMap.remove(x);
                    keySizeMap.remove(x);
                }
                
                if (hRightExists && hRight == x) {
                    hRightExists = false;
                }
                
                if (xParentExists) {
                    assert(!keyParentMap.containsKey(x));
                }
                
                setHFromX(x, h, xVal, hClr, 
                    hLeftExists, hLeft, 
                    hRightExists, hRight, hParentExists, hParent, hSize - 1);
                    
                h = x;
                
                System.out.println("  merged=" + nodeToString(h));
                
                output[0] = 0;
                 
            } else {
                //h.right = delete(h.right, key);
                deleteFromRightAssignRight(h, h, key, output);
                output[0] = 0;
            }
        }
                
        output[0] = 0;
        output[1] = balance(h);  
    }
    
    /*
    private long deleteInternal(long z) {
        RedBlackTreeNode y = null;
        RedBlackTreeNode x = null;
        if (z.left == null || z.right == null) {
            y = z;
        } else {
            y = treeSuccessor(z);
        }
        if (y.left != null) {
            x = y.left;
        } else {
            x = y.right;
        }

        //NLK: x might be null here.  separate into handle children, handle else
        if (x != null) {
            x.parent = y.parent;
        }

        if (y.parent == null) {
            root = x;
        } else if (y == y.parent.left) {
            y.parent.left = x;
        } else {
            y.parent.right = x;
        }

        if (y != z) {
            // copy y data into z
            z.key = y.key;
        }
        //if (y.color == BLACK) {
        //    // NLK:  x cannot be null
        //    if (x != null) {
        //        deleteFixup(x);
        //    }
        //}
        return y;
    }
    */
    
   /***************************************************************************
    *  Red-black tree helper functions.
    ***************************************************************************/

    /**
     * make a left-leaning link lean to the right.
     Note that the parent link logic is from Cormen et al. "Introduction to
     Algorithms".
     */
    protected long rotateRight(long h) {
        
        System.out.println("  *RR " + h);
        
        // assert (h != null) && isRed(h.left);
        assert(keyValMap.containsKey(h));
        assert(isLeftRed(h));
        
        /*
        Node x = h.left;
        h.left = x.right;
        x.right = h;
        x.color = x.right.color;
        x.right.color = RED;
        x.size = h.size;
        h.size = size(h.left) + size(h.right) + 1;
        */
        
        long x = keyLeftMap.get(h);

        //System.out.println("  before RR h=" + nodeToString(h));
        //System.out.println("  before RR x=" + nodeToString(x));
        
        if (keyRightMap.containsKey(x)) {
            //System.out.println("--0");
            //h.left = x.right;
            long xRight = keyRightMap.get(x);
            
            keyLeftMap.put(h, xRight);
            keyParentMap.put(xRight, h);
        } else {
            //System.out.println("--1");
            keyLeftMap.remove(h);
        }
        if (keyParentMap.containsKey(h)) {
            //System.out.println("--2");
            long hParent = keyParentMap.get(h);
            // assign x as child of its new parent
            
            keyParentMap.put(x, hParent);
            if (keyLeftMap.containsKey(hParent) && keyLeftMap.get(hParent) ==
                h) {
                //System.out.println("--3");
                keyLeftMap.put(hParent, x);
            } else {
                //System.out.println("--4");
                keyRightMap.put(hParent, x);
            }
        } else {
            //System.out.println("--5");
            keyParentMap.remove(x);
            //root = x;
        }
        
        //System.out.println("  in RR after h.left h=" + nodeToString(h));
        
        //x.right = h;
        keyRightMap.put(x, h);
        keyParentMap.put(h, x);
        
        //System.out.println("  in RR after x.right h=" + nodeToString(h));
        //System.out.println("  in RR after x.right x=" + nodeToString(x));
        
        //x.color = x.right.color;
        //x.right.color = RED;
        //x.size = h.size;
        keyColorMap.put(x, keyColorMap.get(keyRightMap.get(x)));
        keyColorMap.put(keyRightMap.get(x), RED);
        
        keySizeMap.put(x, keySizeMap.get(h));
        
        //h.size = size(h.left) + size(h.right) + 1;
        int size = sizeLeft(h) + sizeRight(h) + 1;
        keySizeMap.put(h, size);
        
        //System.out.println("  after RR: h=" + nodeToString(h));
        
        //System.out.println("after rotateRight:");
        //printPreOrderTraversal();
        
        return x;
    }

    /**
     make a right-leaning link lean to the left.
     Note that the parent link logic is from Cormen et al. "Introduction to
     Algorithms".
     */
    protected long rotateLeft(long h) {
        
        System.out.println("  *RL h=" + h);
        
        // assert (h != null) && isRed(h.right);
        assert(keyValMap.containsKey(h));
        assert(isRightRed(h));
        
        //Node x = h.right;
        long x = keyRightMap.get(h);
        
        //System.out.println("  before RL h=" + nodeToString(h));
        //System.out.println("  before RL x=" + nodeToString(x));

        //h.right = x.left;
        if (keyLeftMap.containsKey(x)) {
            //System.out.println("---0");
            long left = keyLeftMap.get(x);
            
            keyRightMap.put(h, left);
            keyParentMap.put(left, h);
            
        } else {
            //System.out.println("---1");
            keyRightMap.remove(h);
        }
       
        //System.out.println("  in RL after h.right h=" + nodeToString(h));
        
        if (keyParentMap.containsKey(h)) {
            //System.out.println("---2");
            long hParent = keyParentMap.get(h);
            keyParentMap.put(x, hParent);
            if (keyLeftMap.containsKey(hParent) && keyLeftMap.get(hParent) ==
                h) {
                //System.out.println("---3");
                keyLeftMap.put(hParent, x);
            } else {
                //System.out.println("---4");
                keyRightMap.put(hParent, x);
            }
        } else {
            //System.out.println("---5");
            keyParentMap.remove(x);
        }
                
        //x.left = h;
        keyLeftMap.put(x, h);
        keyParentMap.put(h, x);
        
        //System.out.println("  in RL after x.left h=" + nodeToString(h));
        //System.out.println("  in RL after x.left x=" + nodeToString(x));

        //x.color = x.left.color;
        keyColorMap.put(x, keyColorMap.get(h));
        //x.left.color = RED;
        keyColorMap.put(h, RED);
        
        keySizeMap.put(x, keySizeMap.get(h));
        
        int size = sizeLeft(h) + sizeRight(h) + 1;
        keySizeMap.put(h, size);
               
        //System.out.println("  after RL: h=" + nodeToString(h));
        
        return x;
    }

    // flip the colors of a node and its two children
    protected void flipColors(long h) {
        // h must have opposite color of its two children
        assert(keyValMap.containsKey(h) && keyLeftMap.containsKey(h) 
            && keyRightMap.containsKey(h));
        assert(!isRed(h) &&  isLeftRed(h) &&  isRightRed(h))
        || (isRed(h)  && !isLeftRed(h) && !isRightRed(h));
        
        /*
        h.color = !h.color;
        h.left.color = !h.left.color;
        h.right.color = !h.right.color;
        */
        int clr = keyColorMap.get(h);
        clr ^= 1;
        keyColorMap.put(h, clr);
        
        clr = keyColorMap.get(keyLeftMap.get(h));
        clr ^= 1;
        keyColorMap.put(keyLeftMap.get(h), clr);
        
        clr = keyColorMap.get(keyRightMap.get(h));
        clr ^= 1;
        keyColorMap.put(keyRightMap.get(h), clr);
    }

    // Assuming that h is red and both h.left and h.left.left
    // are black, make h.left or one of its children red.
    protected long moveRedLeft(long h) {
        assert(keyValMap.containsKey(h));
        assert(isRed(h) && !isLeftRed(h) && !isLeftLeftRed(h));

        //System.out.println("moveRedLeft " + nodeToString(h));
        
        flipColors(h);
        if (isRightLeftRed(h)) {
            //h.right = rotateRight(h.right);
            long rKey = rotateRight(keyRightMap.get(h));
            keyRightMap.put(h, rKey);
            keyParentMap.put(rKey, h);
            h = rotateLeft(h);
            flipColors(h);
        }
        return h;
    }

    //move red node down the right spine of the tree
    // Assuming that h is red and both h.right and h.right.left
    // are black, make h.right or one of its children red.
    private long moveRedRight(long h) {
        assert (keyValMap.containsKey(h));
        assert(isRed(h) && !isRightRed(h) && !isRightLeftRed(h));
        
        flipColors(h);
        if (isLeftLeftRed(h)) { 
            h = rotateRight(h);
            flipColors(h);
        }
        
        //System.out.println("after moveRedRight h=" + nodeToString(h));
        //printPreOrderTraversal();
        
        return h;
    }

    // restore red-black tree invariant
    private long balance(long h) {
        
        //System.out.println("balance " + nodeToString(h));
        
        assert (keyValMap.containsKey(h));

        if (isRightRed(h))   {
            h = rotateLeft(h);
        }
        if (isLeftRed(h) && isLeftLeftRed(h)) {
            h = rotateRight(h);
        }
        if (isLeftRed(h) && isRightRed(h))  {
            flipColors(h);
        }

        int size = sizeLeft(h) + sizeRight(h) + 1;
        keySizeMap.put(h, size);
        
        return h;
    }


   /***************************************************************************
    *  Utility functions.
    ***************************************************************************/

    /**
     * Returns the height of the BST (for debugging).
     * @return the height of the BST (a 1-node tree has height 0)
     */
    public int height() {
        return height(root);
    }
    private int height(long x) {
        if (!keyValMap.containsKey(x)) return -1;
        return 1 + Math.max(heightLeft(x), heightRight(x));
    }
    private int heightLeft(long x) {
        if (!keyLeftMap.containsKey(x)) return -1;
        return height(keyLeftMap.get(x));
    }
    private int heightRight(long x) {
        if (!keyRightMap.containsKey(x)) return -1;
        return height(keyRightMap.get(x));
    }

   /***************************************************************************
    *  Ordered symbol table methods.
    ***************************************************************************/

    /**
     * Returns the smallest key in the symbol table.
     * @param output if output[0] == -1 no minimum was present,
     * else output[1] holds the smallest key in the symbol table
     * @throws NoSuchElementException if the symbol table is empty
     */
    public void min(long[] output) {
        if (isEmpty()) throw new NoSuchElementException("called min() with empty symbol table");
        output[0] = 0;
        min(root, output);
    }

    // the smallest key in subtree rooted at x; null if no such key
    private void min(long x, long[] output) { 
        //assert(keyValMap.containsKey(x));
        if (!keyValMap.containsKey(x)) {
            output[0] = -1;
            return;
        }
        while (keyLeftMap.containsKey(x)) {
            x = keyLeftMap.get(x);
        }
        output[1] = x;
    }

    /**
     * Returns the largest key in the symbol table.
     * @param output if output[0] == -1 no minimum was present,
     * else output[1] holds the largest key in the symbol table
     * @throws NoSuchElementException if the symbol table is empty
     */
    public void max(long[] output) {
        if (isEmpty()) throw new NoSuchElementException("called max() with empty symbol table");
        output[0] = 0;
        max(root, output);
    } 

    // the largest key in the subtree rooted at x; null if no such key
    private void max(long x, long[] output) { 
        //assert(keyValMap.containsKey(x));
        if (!keyValMap.containsKey(x)) {
            output[0] = -1;
            return;
        }
        while (keyRightMap.containsKey(x)) {
            x = keyRightMap.get(x);
        }
        output[1] = x;
    }

    /**
     * Returns the largest key in the symbol table less than or equal to {@code key}.
     * @param key the key
     * @param output if output[0] == -1, the key was not present, 
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
        floor(root, key, output);
    }    

    /**
     * Returns the largest key in the symbol table less than {@code key}.
     * @param key the key
     * @param output if output[0] == -1, the key was not present, 
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
        TLongList stack = new TLongArrayList();
        output[0] = 0;
        lower(root, key, stack, output);
    }    
    
    // the largest key in the subtree rooted at x less than or equal to the given key
    private void floor(long x, long key, long[] output) {
        if (!keyValMap.containsKey(x)) {
            output[0] = -1;
            return;
        }
        output[0] = 0;
        int cmp = (key < x) ? -1 : (key > x) ? 1 : 0;
        if (cmp == 0) {
            output[1] = x;
            return;
        }
        if (cmp < 0)  {
            if (keyLeftMap.containsKey(x)) {
                floor(keyLeftMap.get(x), key, output);
                return;
            }
            output[0] = -1;
            return;
        }
        if (keyRightMap.containsKey(x)) {
            floor(keyRightMap.get(x), key, output);
            if (output[0] == -1) {
                output[0] = 0;
                output[1] = x;
            }
            return;
        }
        output[0] = 0;
        output[1] = x;
    }
    
    /** the largest key in the subtree rooted at x less than the given key
     * 
     * The method is symmetric t the successor method called 
     * higher, as suggested by Cormen et al. in
     * the book "Introduction to Algorithms".
     * 
     * @param x
     * @param key
     * @param stack
     * @return 
     */
    private void lower(long x, long key, TLongList stack, long[] output) { 
        
        //binary search until overshoot
        while (keyValMap.containsKey(x) && key != x) {
            stack.add(x);
            //System.out.println("lower: x=" + x + " q=" + key);
            if (key < x) {
                if (!keyLeftMap.containsKey(x)) {
                    break;
                }
                x = keyLeftMap.get(x);
            } else {
                //System.out.println("   x=" + x.key);
                if (!keyRightMap.containsKey(x)) {
                    break;
                }
                x = keyRightMap.get(x);
            }
        }
        
        long y = -1;
        int yIdx = -1;
        boolean yIsSet = false;
        
        if (!keyValMap.containsKey(x)) {
            if (stack.size() > 1) {
                x = stack.get(stack.size() - 1);
                yIdx = stack.size() - 2;
                y = stack.get(yIdx);
                yIsSet = true;
            } else if (stack.size() == 1) {
                x = stack.get(stack.size() - 1);
            }
        } else if (!stack.isEmpty()) {
            yIdx = stack.size() - 1;
            y = stack.get(yIdx);
            yIsSet = true;
        }
        
        if (keyLeftMap.containsKey(x)) {
            max(keyLeftMap.get(x), output);
            return;
        }
      
        //while (y != null && x == y.left) {
        while (yIsSet && 
            (
            (keyLeftMap.containsKey(y) && keyValMap.containsKey(x) &&
            keyLeftMap.get(y) == keyValMap.get(x)) 
            || 
            (!keyLeftMap.containsKey(y) && !keyValMap.containsKey(x))
            )
            ) {
            
            //System.out.println("lower: y=" + y.key + " q=" + key);
            x = y;
            yIdx--;
            if (yIdx < 0) break;
            y = stack.get(yIdx);
            yIsSet = true;
        }
        
        //System.out.println("    y=" + y.key + " q=" + key);
        if (yIsSet && y >= key) {
            output[0] = -1;
            return;
        }
        if (yIsSet) {
            output[0] = 0;
            output[1] = y;
        } else {
            output[0] = -1;
        }
    }

    /**
     * Returns the smallest key in the symbol table greater than or equal to {@code key}.
     * @param key the key
     * @param output if output[0] == -1, the key was not present, 
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
        output[0] = 0;
        ceiling(root, key, output);
    }

    /**
     * Returns the smallest key in the symbol table greater than {@code key}.
     * @param key the key
     * @param output if output[0] == -1, the key was not present, 
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
            throw new NoSuchElementException("called floor() with empty symbol table");
        }
        higher(root, key, output);
    }
    
    // the smallest key in the subtree rooted at x greater than or equal to the given key
    private void ceiling(long x, long key, long[] output) {
        
        //TODO: make this iterative
        
        if (!keyValMap.containsKey(x)) {
            output[0] = -1;
            return;
        }
        output[0] = 0;
        int cmp = key < x ? -1 : (key > x) ? 1 : 0; 
        if (cmp == 0) {
            output[1] = x;
            return;
        }
        if (cmp > 0)  {
            if (keyRightMap.containsKey(x)) {
                ceiling(keyRightMap.get(x), key, output);
            } else {
                output[0] = -1;
            }
            return;
        }
        /*
        Node t = ceiling(x.left, key);
        if (t != null) {
            return t;
        } else {
            return x;
        }*/
        if (keyLeftMap.containsKey(x)) {
            ceiling(keyLeftMap.get(x), key, output);
            if (output[0] == -1) {
                output[0] = 0;
                output[1] = x;
            }
        } else {
            output[0] = 0;
            output[1] = x;
        }
    }
    
    /** the smallest key in the subtree rooted at x greater than the given key.
     * 
     * NOTE: the method uses in part, a pattern adapted from the Cormen et al.
     * book "Introduction to Algorithms" for their Red Black Tree.
     * 
     * @param x
     * @param key
     * @return 
     */
    private void higher(long x, long key, long[] output) {  
        
        //System.out.println("x=" + x + " key=" + key);
                
        /*
                    X
        left .lte.     right .gte.
        */
        
        while (keyValMap.containsKey(x)) {
            int cmp = key < x ? -1 : (key > x) ? 1 : 0; 
            if (cmp < 0) {
                if (keyLeftMap.containsKey(x)) {
                    x = keyLeftMap.get(x);
                } else {
                    break;
                }
            } else if (cmp > 0) {
                if (keyRightMap.containsKey(x)) {
                    x = keyRightMap.get(x);
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        
        if (keyRightMap.containsKey(x)) {
            output[0] = 0;
            min(keyRightMap.get(x), output);
            //System.out.println("min=" + Arrays.toString(output));
            return;
        }
        Long y = null;
        if (keyParentMap.containsKey(x)) {
            y = keyParentMap.get(x);
            //System.out.println("y=" + y);
        }
        while (y != null && keyRightMap.containsKey(y) &&
            x == keyRightMap.get(y)) {
            //System.out.println("y=" + y);
            x = y;
            if (keyParentMap.containsKey(x)) {
                y = keyParentMap.get(x);
            } else {
                y = null;
            }
        }
        if (y == null) {
            output[0] = -1;
        } else {
            //System.out.println("y=" + y);
            output[0] = 0;
            output[1] = y;
        }
    }

    /**
     * Return the kth smallest key in the symbol table.
     * @param k the order statistic
     * @return the {@code k}th smallest key in the symbol table
     * @throws IllegalArgumentException unless {@code k} is between 0 and
     *     <em>n</em>–1
     */
    public long select(int k) {
        if (k < 0 || k >= size()) {
            throw new IllegalArgumentException("called select() with invalid argument: " + k);
        }
        long[] output = new long[2];
        select(root, k, output);
        if (output[0] == -1) {
            return -1;
        }
        return output[1];
    }

    // the key of rank k in the subtree rooted at x
    private void select(long x, int k, long[] output) {
        
        //TODO: make this iterative
        
        assert(keyValMap.containsKey(x));
        assert(k >= 0 && k < size(x));
        
        output[0] = 0;
        
        int t = sizeLeft(x); 
        if (t > k) {
            select(keyLeftMap.get(x),  k, output);
        } else if (t < k) {
            if (keyRightMap.containsKey(x)) {
                select(keyRightMap.get(x), k-t-1, output);
            } else {
                output[0] = -1;
            }
        } else {
            output[1] = x;
        }
    } 

    /**
     * Return the number of keys in the symbol table strictly less than {@code key}.
     * @param key the key
     * @return the number of keys in the symbol table strictly less than {@code key}
     * @throws IllegalArgumentException if {@code key} is {@code null}
     */
    public int rank(long key) {
        int[] output = new int[2];
        rank(key, root, output);
        if (output[0] == -1) {
            return 0;
        }
        return output[1];
    } 

    // number of keys less than key in the subtree rooted at x
    private void rank(long key, long x, int[] output) {
        
        //TODO: make this iterative
        output[0] = 0;
        
        if (!keyValMap.containsKey(x)) {
            output[1] = 0;
            return;
        } 
        int cmp = key < x ? -1 : (key > x) ? 1 : 0;  
        if (cmp < 0) {
            if (keyLeftMap.containsKey(x)) {
                rank(key, keyLeftMap.get(x), output);
            } else {
                output[1] = 0;
            }
        } else if (cmp > 0) {
            //1 + size(x.left) + rank(key, x.right);
            int sz;
            if (keyRightMap.containsKey(x)) {
                output[0] = 0;
                rank(key, keyRightMap.get(x), output);
                if (output[0] == -1) {
                    sz = 1 + sizeLeft(x);
                } else {
                    sz = 1 + sizeLeft(x) + output[1];
                }
            } else {
                sz = 1 + sizeLeft(x);
            }
            output[1] = sz;
        } else {
            int sz = sizeLeft(x);
            output[1] = sz;
        } 
    } 

   /***************************************************************************
    *  Range count and range search.
    ***************************************************************************/

    /**
     * Returns all keys in the symbol table as an {@code Iterable}.
     * To iterate over all of the keys in the symbol table named {@code st},
     * use the foreach notation: {@code for (Key key : st.keys())}.
     * @return all keys in the symbol table as an {@code Iterable}
     */
    public TLongList keys() {
        
        if (isEmpty()) return new TLongArrayList();
        
        long[] output = new long[2];
        min(output);
        assert(output[0] != -1);
        long lo = output[1];
        
        max(output);
        assert(output[0] != -1);
        long hi = output[1];
        
        return keys(lo, hi);
    }

    /**
     * Returns all keys in the symbol table in the given range,
     * as an {@code Iterable}.
     *
     * @param  lo minimum endpoint
     * @param  hi maximum endpoint
     * @return all keys in the sybol table between {@code lo} 
     *    (inclusive) and {@code hi} (inclusive) as an {@code Iterable}
     * @throws IllegalArgumentException if either {@code lo} or {@code hi}
     *    is {@code null}
     */
    public TLongList keys(long lo, long hi) {

        TLongList queue = new TLongArrayList();
        if (isEmpty() || lo > hi) return queue;
        
        keys(root, queue, lo, hi);
        
        return queue;
    } 

    // add the keys between lo and hi in the subtree rooted at x
    // to the queue
    private void keys(long x, TLongList queue, long lo, long hi) {
        
        //TODO: make this iterative
        
        if (!keyValMap.containsKey(x)) {
            return;
        } 
        int cmplo = lo < x ? -1 : (lo > x) ? 1 : 0;
        int cmphi = hi < x ? -1 : (hi > x) ? 1 : 0; 
        if (cmplo < 0) {
            if (keyLeftMap.containsKey(x)) {
                keys(keyLeftMap.get(x), queue, lo, hi);
            }
        } 
        if (cmplo <= 0 && cmphi >= 0) {
            queue.add(x);
        } 
        if (cmphi > 0) {
            if (keyRightMap.containsKey(x)) {
                keys(keyRightMap.get(x), queue, lo, hi);
            }
        } 
    } 

    /**
     * Returns the number of keys in the symbol table in the given range.
     *
     * @param  lo minimum endpoint
     * @param  hi maximum endpoint
     * @return the number of keys in the sybol table between {@code lo} 
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
    ***************************************************************************/
    private boolean check() {
        boolean t1 = isParentChildConsistent();
        boolean t2 = isBST();
        boolean t3 = isSizeConsistent();
        boolean t4 = isRankConsistent();
        boolean t5 = is23();
        boolean t6 = isBalanced();
        if (!t1) System.out.println("Not consistent parent child relationships");
        if (!t2) System.out.println("Not in symmetric order");
        if (!t3) System.out.println("Subtree counts not consistent");
        if (!t4) System.out.println("Ranks not consistent");
        if (!t5) System.out.println("Not a 2-3 tree");
        if (!t6) System.out.println("Not balanced");
        System.out.flush();
        return t1 && t2 && t3 && t4 && t5 && t6;
    }

    // does this binary tree satisfy symmetric order?
    // Note: this test also ensures that data structure is a binary tree since order is strict
    private boolean isBST() {
        return isBST(root, null, null);
    }

    // is the tree rooted at x a BST with all keys strictly between min and max
    // (if min or max is null, treat as empty constraint)
    // Credit: Bob Dondero's elegant solution
    private boolean isBST(long x, Long min, Long max) {
        if (!keyValMap.containsKey(x)) {
            return true;
        }
        if (min != null && x <= min.longValue()) {
            return false;
        }
        if (max != null && x >= max.longValue()) {
            return false;
        }
        Long key = Long.valueOf(x);
        return isBSTLeft(x, min, key) && isBSTRight(x, key, max);
    }
    private boolean isBSTLeft(long x, Long min, Long max) {
        if (!keyLeftMap.containsKey(x)) {
            return true;
        }
        long key = keyLeftMap.get(x);
        return isBST(key, min, max);
    }
    private boolean isBSTRight(long x, Long min, Long max) {
        if (!keyRightMap.containsKey(x)) {
            return true;
        }
        long key = keyRightMap.get(x);
        return isBST(key, min, max);
    }

    private boolean isSizeConsistent() { return isSizeConsistent(root); }
    private boolean isSizeConsistent(long x) {
        if (!keyValMap.containsKey(x)) return true;
        if (keySizeMap.get(x) != sizeLeft(x) + sizeRight(x) + 1) return false;
        return isSizeConsistentLeft(x) && isSizeConsistentRight(x);
    }
    private boolean isSizeConsistentLeft(long x) {
        if (!keyLeftMap.containsKey(x)) return true;
        long key = keyLeftMap.get(x);
        return isSizeConsistent(key);
    }
    private boolean isSizeConsistentRight(long x) {
        if (!keyRightMap.containsKey(x)) return true;
        long key = keyRightMap.get(x);
        return isSizeConsistent(key);
    }

    // check that ranks are consistent
    private boolean isRankConsistent() {
        for (int i = 0; i < size(); i++) {
            int r = rank(select(i));
            //System.out.println("i=" + i + " r=" + r + " size=" + size());
            if (i != r) return false;
        }
        TLongList keys = keys();
        for (int i = 0; i < keys.size(); ++i) {
            long key = keys.get(i);
            int r = rank(key);
            long s = select(r);
            //System.out.println("i=" + i + " key=" + key + " r=" + r + " s=" + s);
            if (key != s) {
                return false;
            }
        }
        return true;
    }

    // Does the tree have no red right links, and at most one (left)
    // red links in a row on any path?
    private boolean is23() { return is23(root); }
    private boolean is23(long x) {
        if (!keyValMap.containsKey(x)) {
            return true;
        }
        if (isRightRed(x)) {
            return false;
        }
        if (x != root && isRed(x) && isLeftRed(x)) {
            return false;
        }
        return is23Left(x) && is23Right(x);
    }
    private boolean is23Left(long x) {
        if (!keyLeftMap.containsKey(x)) return true;
        return is23(keyLeftMap.get(x));
    }
    private boolean is23Right(long x) {
        if (!keyRightMap.containsKey(x)) return true;
        return is23(keyRightMap.get(x));
    }

    // do all paths from root to leaf have same number of black edges?
    private boolean isBalanced() { 
        if (!rootIsSet) {
            return true;
        }
        int black = 0;     // number of black links on path from root to min
        long x = root;
        while (keyValMap.containsKey(x)) {
            if (!isRed(x)) {
                black++;
            }
            if (!keyLeftMap.containsKey(x)) {
                break;
            }
            x = keyLeftMap.get(x);
        }
        return isBalanced(root, black);
    }

    // does every path from the root to a leaf have the given number of black links?
    private boolean isBalanced(long x, int black) {
        if (!keyValMap.containsKey(x)) {
            return black == 0;
        }
        if (!isRed(x)) {
            black--;
        }
        return isLeftBalanced(x, black) && isRightBalanced(x, black);
    }
    private boolean isLeftBalanced(long x, int black) {
        if (!keyLeftMap.containsKey(x)) {
            return black == 0;
        }
        long key = keyLeftMap.get(x);
        return isBalanced(key, black);
    } 
    private boolean isRightBalanced(long x, int black) {
        if (!keyRightMap.containsKey(x)) {
            return black == 0;
        }
        long key = keyRightMap.get(x);
        return isBalanced(key, black);
    } 
    
    /**
     * left subtree, root, right subtree
     */
    public void printInOrderTraversal() {
        System.out.print("root=");
        if (rootIsSet) {
            System.out.print(root);
        }
        System.out.println("");
        long[] nodes = getInOrderTraversalIterative(root);
        for (long node : nodes) {
            System.out.println("node=" + nodeToString(node));
        }
    }
    
    /**
     * root, left subtree, right subtree
     */
    public void printPreOrderTraversal() {
        System.out.print("root=");
        if (rootIsSet) {
            System.out.print(root);
        }
        System.out.println(" size=" + size());
        long[] nodes = getPreOrderTraversalIterative(root, 1);
        for (long node : nodes) {
            System.out.println("  node=" + nodeToString(node));
        }
    }
    public void printPreOrderTraversal2(long topNode) {
        System.out.print("root=");
        if (rootIsSet) {
            System.out.print(root);
        }
        System.out.println(" size=" + size());
        long[] nodes = getPreOrderTraversalIterative(topNode, 1);
        for (long node : nodes) {
            System.out.println("  node=" + nodeToString(node));
        }
    }
    private void printPreOrderTraversal(int addExtraToSize) {
        long[] nodes = getPreOrderTraversalIterative(root, 
            addExtraToSize);
        System.out.print("root=");
        if (rootIsSet) {
            System.out.print(root);
        }
        System.out.println("");
        for (long node : nodes) {
            System.out.println("  node=" + nodeToString(node));
        }
    }
    
    /**
     * left subtree, right subtree, root subtree
     */
    public void printPostOrderTraversal() {
        long[] nodes = getPostOrderTraversalIterative(root);
        for (long node : nodes) {
            System.out.println("  node=" + nodeToString(node));
        }
    }

    /**
     * visit each node using pattern left subtree, root, right subtree
     * in an iterative manner rather than invoking the method recursively.
     */
    protected long[] getInOrderTraversalIterative(Long node) {
       
        if (isEmpty()) {
            return new long[0];
        }
        
        int sz = size();
        
        long[] array = new long[sz];
        int count = 0;
        
        Stack<Long> stack = new Stack<>();
               
        while (!stack.isEmpty() || (node != null)) {
            if (node != null) {
                 
                stack.push(node);
                
                node = keyLeftMap.containsKey(node.longValue()) ?
                    keyLeftMap.get(node.longValue()) : null;
            
            } else if (count < sz) {
                
                node = stack.pop();
                
                array[count] = node;
                count++;
                
                //System.out.println(node.key);
                
                node = keyRightMap.containsKey(node.longValue()) ?
                    keyRightMap.get(node.longValue()) : null;
            }
        }
        if (count < sz) {
            // can happen during debugging when insert is not complete yet
            array = Arrays.copyOf(array, count);
        }
        
        return array;
    }
    
    /**
     * visit each node using pattern: root, left subtree, right subtree
     * in an iterative manner rather than invoking the method recursively.
     */
    protected long[] getPreOrderTraversalIterative(Long node, int addExtraToSize) {
       
        //NOTE: added additional integer and conditions 
        //   for size because may be printing tree
        //   in the middle of a put where the node size is not yet updated.
        // The count conditionals below are otherwise, not needed.
        
        if (isEmpty()) {
            return new long[0];
        }
                
        int sz = size(node) + addExtraToSize;
        
        long[] array = new long[sz];
        int count = 0;
        
        Stack<Long> stack = new Stack<>();
        
        while (count < sz && (!stack.isEmpty() || node != null)) {
            if (node != null && count < sz) {
                
                array[count] = node;
                count++;
                //System.out.println(node);
                
                if (count < sz) {
                    stack.push(node);
                }
                
                node = keyLeftMap.containsKey(node.longValue()) ?
                    keyLeftMap.get(node.longValue()) : null;
            
            } else if (count < sz) {
                
                node = stack.pop();
                
                node = keyRightMap.containsKey(node.longValue()) ?
                    keyRightMap.get(node.longValue()) : null;
            }
        }
        
        if (count < sz) {
            array = Arrays.copyOf(array, count);
        }
        
        return array;
    }

    /**
     * visit each node using pattern: left subtree, right subtree, root subtree
     * in an iterative manner rather than invoking the method recursively.
     */
    protected long[] getPostOrderTraversalIterative(Long node) {
    
        if (isEmpty()) {
            return new long[0];
        }
        
        int sz = size();
        
        long[] array = new long[sz];
        int count = 0;
        
        if (node == null) {
            return array;
        }
        
        Stack<Long> stack = new Stack<Long>();
        Stack<Long> stack2 = new Stack<Long>();
        stack.push(node);
        
        while (!stack.isEmpty()) {
            
            node = stack.pop();
            
            stack2.push(node);
            
            if (keyLeftMap.containsKey(node.longValue())) {
                stack.push(keyLeftMap.get(node.longValue()));
            }

            if (keyRightMap.containsKey(node.longValue())) {
                stack.push(keyRightMap.get(node.longValue()));
            }
        }
        
        while (!stack2.isEmpty() && count < sz) {
            
            node = stack2.pop();
            
            //process(node);
            array[count] = node;
            count++;
            //System.out.println(node);
        }
        
        if (count < sz) {
            // can happen during debugging when insert is not complete yet
            array = Arrays.copyOf(array, count);
        }
         
        return array;
    }
   
    /**
     * estimate the size that an instance of RedBlackBSTLongInt with
     * n entries would occupy in heap space in Bytes.
     * 
     * @param numberOfEntries amount of space for this object's instance
     * with n entries in Bytes on the heap.
     * 
     * @return 
     */
    public static long estimateSizeOnHeap(int numberOfEntries) {
        
        long total = 0;
       
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNBooleanFields(1);
        est.setNLongFields(1);
        est.setNArrayRefsFields(1);
        est.setNIntFields(4);
       
        total += est.estimateSizeOnHeap();
        
        // add the maps
        total += 5 * ObjectSpaceEstimator.estimateTLongLongHashMap();
        
        // add the entries in the maps
        // 3 maps of longlong, 3 maps of longint
        // 9 longs, 3 ints
        
        total += numberOfEntries * 
            (9 * ObjectSpaceEstimator.estimateLongSize()
            + 3 * ObjectSpaceEstimator.estimateIntSize());
            
        return total;
    }
    
    private String nodeToString(long key) {
        //assert(keyValMap.containsKey(key));
        if (!keyValMap.containsKey(key)) {
            // this can happen in the middle of a method, for example,
            // when root has been removed and new is not yet assigned
            System.out.println("ERROR: key " + key + " not in maps");
            return "";
        }
        
        StringBuilder sb = new StringBuilder();
        //node=key=0 val=0 color=false size=1
        sb.append("key=").append(key).append(" val=").append(keyValMap.get(key));
        sb.append(" color=").append(keyColorMap.get(key));
        sb.append(" size=").append(keySizeMap.get(key));
        sb.append(" p=");
        if (keyParentMap.containsKey(key)) {
           sb.append(keyParentMap.get(key));
        }
        sb.append(" l=");
        if (keyLeftMap.containsKey(key)) {
           sb.append(keyLeftMap.get(key));
        }
        sb.append(" r=");
        if (keyRightMap.containsKey(key)) {
           sb.append(keyRightMap.get(key));
        }
        return sb.toString();
    }

    private boolean isParentChildConsistent() {
        
        boolean passed = true;
        
        if (!rootIsSet) {
            if (!(keyValMap.isEmpty() && keyColorMap.isEmpty() && 
                keyParentMap.isEmpty() && keyLeftMap.isEmpty() &&
                keyRightMap.isEmpty() && keySizeMap.isEmpty())) {
                System.err.println("maps not empty, but root is");
                passed = false;
            }
        }
        
        //System.out.println("root=" + nodeToString(root));
        
        if (keyParentMap.containsKey(root)) {
            System.err.println("root should not have parent key");
            passed = false;
        }
        
        long[] nodes = getPreOrderTraversalIterative(root, 0);
        for (long key : nodes) {
            if (keyLeftMap.containsKey(key)) {
                long child = keyLeftMap.get(key);
                if (!keyParentMap.containsKey(child)) {
                    System.err.format(
                        "error in %d.left=%d has no parent\n",
                        key, child);
                    passed = false;
                } else if (keyParentMap.get(child) != key) {
                    System.err.format(
                        "error in left: %d.left=%d but %d.parent=%d\n",
                        key, child, child, keyParentMap.get(child));
                    passed = false;
                }
            }
            if (keyRightMap.containsKey(key)) {
                long child = keyRightMap.get(key);
                if (!keyParentMap.containsKey(child)) {
                    System.err.format(
                        "error in %d.right=%d has no parent\n",
                        key, child);
                    passed = false;
                } else if (keyParentMap.get(child) != key) {
                    System.err.format(
                        "error in right: %d.right=%d but %d.parent=%d\n",
                        key, child, child, keyParentMap.get(child));
                    passed = false;
                }
            }
        }
        if (nodes.length != size()) {
            System.out.println("ERROR in nodes extraction");
            passed = false;
        }
        
        return passed;
    }
    
    /**
     Print node and its left and right subtrees, but note that one can
     only print if the height and maxValue of the numbers fit into a line 
     limited to 100 characters.
     Each node takes w=log10(maxValue)+1 characters + 1 space.
     The number of characters of the bottom leaves should be.lt. 100 total.
         (w+1) * (1 left-shift (h-1)) .lt. 100.
      
     The restriction is to make it easy to read an ascii tree on a text terminal. 
     
     Example use: for maxValue=99, the number of base-10 digits is 2, so
     the maximum height printable by this method would be a tree 
     with 6 levels (results in leaf level using .lte.  100 characters).
     
     <pre>
                1              
          2           2       
       3     3     3     3     
      4  4  4  4  4  4  4  4   
     </pre>
     @param node
     @param maxValue
     */
    public void printSmallTree(long node, long maxValue) {
       
        int w = (int)Math.ceil(Math.log(maxValue)/Math.log(10));
        
        int n = keySizeMap.get(node);
        int h = (int)Math.ceil(Math.log(n + 1)/Math.log(2));
        
        //System.out.println("w=" + w + " n=" + n + " h=" + h);
        
        int baselineLength = (w + 1) * (1 << (h - 1));
       
        if (baselineLength > 100) {
            throw new IllegalArgumentException("the number of characters needed"
                + " for the leaves is > 100");
        }
        
        ArrayDeque<Long> levelQ = new ArrayDeque();
        ArrayDeque<Long> nextLevelQ = new ArrayDeque<Long>();
        levelQ.add(node);
        int level = 0;
        int hw = (w+1)/2;
        int indent = (baselineLength - (w+1))/2;
        while (!levelQ.isEmpty() || !nextLevelQ.isEmpty()) {
            level++;
            if (level > h) {
                break;
            }
            if (level > 1) {
                indent /= 2;
                if (indent - hw > 1) {
                    indent -= hw;
                }
            }
            int nn = 1 << (level - 1);
            int nodeSpace = (w+1) * nn;
            int spacing = baselineLength - nodeSpace - 2*indent;
            if ((nn - 1) > 0) {
                spacing /= (nn-1);
            }
            //System.out.println("spacing=" + spacing + " indent=" + indent
            //    + " bl=" + baselineLength + " ns=" + nodeSpace);
            int prevPos = -1;
            StringBuilder sb0 = new StringBuilder(baselineLength);
            StringBuilder sb1 = new StringBuilder(baselineLength);
            while (!levelQ.isEmpty()) {
                long z = levelQ.pop();
                int nSpaces;
                if (level == 1) {
                    nSpaces = indent;
                } else if (prevPos == -1) {
                    nSpaces = indent;
                } else {
                    nSpaces = spacing;
                }
                addSpaces(nSpaces, sb0);
                addSpaces(nSpaces, sb1);
                prevPos += nSpaces;
                if (z == Long.MIN_VALUE) {
                    addSpaces(w + 1, sb0);
                    addSpaces(w + 1, sb1);
                    // add 2 empty placeholders
                    nextLevelQ.add(Long.MIN_VALUE);
                    nextLevelQ.add(Long.MIN_VALUE);
                    continue;                    
                }
                if (!keyValMap.containsKey(z)) {
                    continue;
                }
                int clr = keyColorMap.get(z);
                String keyC = Long.toString(z);
                String clrC = Integer.toString(clr);
            
                sb0.append(keyC);
                sb1.append(clrC);
                addSpaces(1, sb0);
                addSpaces(w + 1 - clrC.length(), sb1);
                
                if (keyLeftMap.containsKey(z)) {
                    nextLevelQ.add(keyLeftMap.get(z));
                } else {
                    nextLevelQ.add(Long.MIN_VALUE);
                }
                if (keyRightMap.containsKey(z)) {
                    nextLevelQ.add(keyRightMap.get(z));
                } else {
                    nextLevelQ.add(Long.MIN_VALUE);
                }
            }
            System.out.println(sb0.toString());
            System.out.println(sb1.toString());
            System.out.println("");
            assert(levelQ.isEmpty());
            levelQ.addAll(nextLevelQ);
            nextLevelQ.clear();
        }
    } 
    private void addSpaces(int nSpaces, StringBuilder sb) {
        for (int i = 0; i < nSpaces; ++i) {
            sb.append(" ");
        }
    }
    
    private boolean containsLeft(long key) {
        return keyLeftMap.containsKey(key);
    }
    private boolean containsRight(long key) {
        return keyRightMap.containsKey(key);
    }
    private boolean containsLeftLeft(long key) {
        if (containsLeft(key)) {
            return containsLeft(keyLeftMap.get(key));
        }
        return false;
    }
    private boolean containsLeftRight(long key) {
        if (containsLeft(key)) {
            return containsRight(keyLeftMap.get(key));
        }
        return false;
    }
    private boolean containsLeftLeftLeft(long key) {
        if (containsLeft(key)) {
            return containsLeftLeft(keyLeftMap.get(key));
        }
        return false;
    }
    private boolean containsLeftLeftRight(long key) {
        if (containsLeft(key)) {
            return containsLeftRight(keyLeftMap.get(key));
        }
        return false;
    }
    private boolean containsParent(long key) {
        return keyParentMap.containsKey(key);
    }
    private boolean containsParentLeft(long key) {
        if (containsParent(key)) {
            return containsLeft(keyParentMap.get(key));
        }
        return false;
    }
    private boolean containsParentRight(long key) {
        if (containsParent(key)) {
            return containsRight(keyParentMap.get(key));
        }
        return false;
    }
    private boolean containsParentParent(long key) {
        if (containsParent(key)) {
            return containsParent(keyParentMap.get(key));
        }
        return false;
    }
    private boolean containsParentParentLeft(long key) {
        if (containsParent(key)) {
            return containsParentLeft(keyParentMap.get(key));
        }
        return false;
    }
    private boolean containsParentParentRight(long key) {
        if (containsParent(key)) {
            return containsParentRight(keyParentMap.get(key));
        }
        return false;
    }

    private void deleteFromLeftAssignLeft(long delH, long asnH, long key, 
        long[] output) {
       
        delete(keyLeftMap.get(delH), key, output);
        if (output[0] == -1) {
            keyLeftMap.remove(asnH);
        } else {
            keyLeftMap.put(delH, output[1]);
            keyParentMap.put(output[1], delH);
        }
    }

    /**
     * Note: output[0] == -1 for these cases:
     *     (1) h.right == min and min has no children.
     *         in this case min and its links are removed in method.
     *     (2)
     * else, output[1] is the value h.right and min and its links have
     * been removed in method.
     * @param hDel
     * @param min
     * @param key
     * @param output 
     */
    private void deleteMinFromRightAssignRight(
        long hDel, long min, long[] output) {
       
        //TODO: delete this method
        
        assert(!keyLeftMap.containsKey(min));
        assert (!keyRightMap.containsKey(min));
        assert(keyRightMap.containsKey(hDel));
        
        long hRight = keyRightMap.get(hDel);
        
        //h.right = deleteMin(h.right);
        
        if (hRight == min) {
            if (keyParentMap.containsKey(min)) {
                long minParent = keyParentMap.get(min);
                if (keyLeftMap.containsKey(minParent) 
                    && keyLeftMap.get(minParent) == min) {
                    keyLeftMap.remove(minParent);
                } else {
                    assert(keyRightMap.containsKey(minParent) 
                    && keyRightMap.get(minParent) == min);
                    keyRightMap.remove(minParent);
                }
                keyParentMap.remove(min);
            }
            keyValMap.remove(min);
            keyColorMap.remove(min);
            keySizeMap.remove(min);
            output[0] = -1;
            return;            
        }
        
        deleteMin(hRight, output);
        System.out.println("deleteMin=" + Arrays.toString(output));
        assert(output[0] != -1);
        assert(output[1] == hRight);
        //assert(!keyValMap.containsKey(min));
               
        if (keyParentMap.containsKey(min)) {
            long hParent = keyParentMap.get(min);
            keyRightMap.remove(hParent);
            keyParentMap.remove(min);
        }
        keyValMap.remove(min);
        keyColorMap.remove(min);
        keySizeMap.remove(min);

        keyRightMap.put(hDel, hRight);
        keyParentMap.put(hRight, hDel);
        
    }

    private void deleteFromRightAssignRight(long hDel, long hAsn, long key, 
        long[] output) {
        
        assert(keyRightMap.containsKey(hDel));
        
        //h.right = delete(h.right, key);

        delete(keyRightMap.get(hDel), key, output);
        
        if (output[0] == -1) {
            keyRightMap.remove(hAsn);
        } else {
            keyRightMap.put(hAsn, output[1]);
            keyParentMap.put(output[1], hAsn);
        }
    }
             
    private void setHFromX(long x, long h, int xVal, int hClr, 
        boolean hLeftExists, long hLeft, 
        boolean hRightExists, long hRight, 
        boolean hParentExists, long hParent, int hSize) {
        
        assert(x != h);
        
        if (hParentExists) {
            keyParentMap.put(x, hParent);
            if (keyLeftMap.containsKey(hParent) 
                && keyLeftMap.get(hParent) == h) {
                keyLeftMap.put(hParent, x);
            } else {
                assert(keyRightMap.containsKey(hParent) 
                    && keyRightMap.get(hParent) == h);
                keyRightMap.put(hParent, x);
            }
            keyParentMap.remove(h);
        } else {
            //h is root
        }
        
        if (hLeftExists) {
            keyParentMap.put(hLeft, x);
            keyLeftMap.put(x, hLeft);
            keyLeftMap.remove(h);
        }
        if (hRightExists) {
            keyParentMap.put(hRight, x);
            keyRightMap.put(x, hRight);
            keyRightMap.remove(h);
        }
        keyValMap.put(x, xVal);
        keyColorMap.put(x, hClr);
        keySizeMap.put(x, hSize);
     
        keyValMap.remove(h);
        keyColorMap.remove(h);
        keySizeMap.remove(h); 
        
        if (h == root) {
            root = x;
        }
    }
    
}
