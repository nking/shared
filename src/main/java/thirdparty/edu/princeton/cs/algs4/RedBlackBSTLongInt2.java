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
import algorithms.util.NodeMap;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TLongArrayList;
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
    protected final NodeMap nodeMap;
    
    /**
     * Initializes an empty symbol table.
     */
    public RedBlackBSTLongInt2() {
        nodeMap = new NodeMap();
    }
    
    /**
     * Initializes an empty symbol table.
     */
    public RedBlackBSTLongInt2(int capacity) {
        nodeMap = new NodeMap(capacity);
    }

    private long addNewNode(long key, int val, int color, int size) {
        nodeMap.put(key, val, color, size);
        return key;
    }
    
   /***************************************************************************
    *  Node helper methods.
    ***************************************************************************/
    // is node x red; false if x is null ?
    private boolean isRed(long x) {
        if (!nodeMap.containsKey(x)) return false;
        return nodeMap.getNodeColor(x) == RED;
    }
    private boolean isLeftRed(long x) {
        if (!nodeMap.leftIsSet(x)) {
            return false;
        }
        return isRed(nodeMap.getLeft(x));
    }
    private boolean isRightRed(long x) {
        if (!nodeMap.rightIsSet(x)) {
            return false;
        }
        return isRed(nodeMap.getRight(x));
    }
    private boolean isLeftLeftRed(long x) {
        if (!nodeMap.leftIsSet(x)) {
            return false;
        }
        long left = nodeMap.getLeft(x);
        return isLeftRed(left);
    }
    private boolean isRightLeftRed(long x) {
        if (!nodeMap.rightIsSet(x)) {
            return false;
        }
        long right = nodeMap.getRight(x);
        return isLeftRed(right);
    }

    // number of node in subtree rooted at x; 0 if x is null
    private int size(long x) {
        if (!nodeMap.containsKey(x)) return 0;
        return nodeMap.getNodeSize(x);
    }
    private int sizeLeft(long x) {
        if (!nodeMap.containsKey(x)) return 0;
        if (!nodeMap.leftIsSet(x)) return 0;
        return nodeMap.getNodeSize(nodeMap.getLeft(x));
    }
    private int sizeRight(long x) {
        if (!nodeMap.containsKey(x)) return 0;
        if (!nodeMap.rightIsSet(x)) return 0;
        return nodeMap.getNodeSize(nodeMap.getRight(x));
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
        while (nodeMap.containsKey(x)) {
            int cmp = key < x ? -1 : (key > x) ? 1 : 0; 
            if (cmp < 0) {
                if (nodeMap.leftIsSet(x)) {
                    x = nodeMap.getLeft(x);
                } else {
                    break;
                }
            } else if (cmp > 0) {
                if (nodeMap.rightIsSet(x)) {
                    x = nodeMap.getRight(x);
                } else {
                    break;
                }
            } else {
                output[1] = nodeMap.getNodeValue(x);
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
        nodeMap.updateNodeColor(root, BLACK);
        
        if (nodeMap.parentIsSet(root)) {
            nodeMap.unsetParent(root);
        }
        
        //printPreOrderTraversal();
        //System.out.println("after put " + key + " root=" + root);
        
        assert(check());
    }

    // insert the key-value pair in the subtree rooted at h
    //private Node put(Node h, long key, int val) {
    private long put(long h, long key, int val) {
        
        //System.out.println("put h=" + h + " key=" + key);
        
        if (!rootIsSet || !nodeMap.containsKey(h)) {
            return addNewNode(key, val, RED, 1);
        }
       
        int cmp = (key < h) ? -1 : (key > h) ? 1 : 0;
        if (cmp < 0) {
            //h.left  = put(h.left,  key, val);
            long putKey;
            if (nodeMap.leftIsSet(h)) {
                //h.left = putKey
                //putKey.parent = h
                long hLeft = nodeMap.getLeft(h);
                putKey = put(hLeft,  key, val);
            } else {
                putKey = addNewNode(key, val, RED, 1);
            }
            nodeMap.updateLeft(h, putKey);
            nodeMap.updateParent(putKey, h);
        } else if (cmp > 0) {
            //h.right = put(h.right, key, val);
            long putKey;
            if (nodeMap.rightIsSet(h)) {
                long hRight = nodeMap.getRight(h);
                putKey = put(nodeMap.getRight(h),  key, val);
            } else {
                putKey = addNewNode(key, val, RED, 1);
            }
            nodeMap.updateRight(h, putKey);
            nodeMap.updateParent(putKey, h);
        } else {
            //h.val   = val;
            assert(nodeMap.containsKey(h));
            nodeMap.updateNodeValue(h, val);
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
        nodeMap.updateNodeSize(h, size);
       
        return h;
    }

   /***************************************************************************
    *  Red-black tree deletion.
    ***************************************************************************/

    private void setRootToRedIfChildrenAreBlack() {
        if (!isLeftRed(root) && !isRightRed(root)) {
            nodeMap.updateNodeColor(root, RED);
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
        
        assert(nodeMap.containsKey(root));
        
        // if both children of root are black, set root to red
        setRootToRedIfChildrenAreBlack();
        
        //root = deleteMin(root);
        long[] output = new long[2];
        deleteMin(root, output);
        
        //System.out.println("deleteMin=" + Arrays.toString(output));
        
        if (output[0] == -1) {
            rootIsSet = false;
            nodeMap.remove(root);
            root = -1;
        } else {
            //System.out.println("in deleteMin: assigning root=" + output[1]);
            root = output[1];
            nodeMap.unsetParent(root);
        }
        if (!isEmpty()) {
            //root.color = BLACK;
            nodeMap.updateNodeColor(root, BLACK);
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
        
        if (!nodeMap.containsKey(h)) {
            output[0] = -1;
            return;
        } else if (!nodeMap.leftIsSet(h)) {
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
        
        if (!nodeMap.leftIsSet(h)) {
            output[0] = -1;
            return;
        }
        long left = nodeMap.getLeft(h);
        //assert(parent == h);
        
        if (output.length == 3) {
            output[2] = left;
        }
        
        deleteMin(left, output);
        if (output[0] == -1) {
            //h.left = null
            nodeMap.unsetLeft(h);
            output[0] = 0;
        } else {
            //System.out.println(" deleteMin return is " + output[1] 
            //   + " h.left gets assigned it");
            nodeMap.updateParent(output[1], h);
            nodeMap.updateLeft(h, output[1]);
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
            nodeMap.remove(root);
            rootIsSet = false;
            root = -1;
            return;
        }
        root = output[1];
        nodeMap.unsetParent(root);
        
        if (!isEmpty()) {
            //root.color = BLACK;
            nodeMap.updateNodeColor(root, BLACK);
        }
        
        //System.out.format("AFTER deleteMax()\n");
        //printPreOrderTraversal(1);
        
        assert(check());
    }

    // delete the key-value pair with the maximum key rooted at h
    private void deleteMax(long h, long[] output) { 
    
        //System.out.format("BEFORE deleteMax(%d)\n", h);
        //printPreOrderTraversal(1);
        
        if (isLeftRed(h)) {
            h = rotateRight(h);
        }
      
        if (!nodeMap.rightIsSet(h)) {
            // h is max key
            output[0] = -1;
            return;
        }

        if (!isRightRed(h) && !isRightLeftRed(h)) {
            //move red node down the right spine of the tree
            h = moveRedRight(h);
        }
        
        //h.right = deleteMax(h.right);
        
        assert(nodeMap.rightIsSet(h));
        
        deleteMax(nodeMap.getRight(h), output);
        if (output[0] == -1) {
            nodeMap.unsetRight(h);
        } else {
            nodeMap.updateParent(output[1], h);
            nodeMap.updateRight(h, output[1]);
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
        
        //System.out.println("\nbefore delete " + key + " root=" + root);
        //printPreOrderTraversal();
        
        if (!contains(key)) return;

        // if both children of root are black, set root to red
        setRootToRedIfChildrenAreBlack();
        
        ////root = delete(root, key);
        
        long[] output = new long[2];
        delete(root, key, output);
        
        if (output[0] == -1) {
            rootIsSet = false;
            nodeMap.remove(root);
            root = -1;
            return;
        }
        root = output[1];
        nodeMap.unsetParent(root);
        
        if (!isEmpty()) {
            nodeMap.updateNodeColor(root, BLACK);
        }
        
        //printPreOrderTraversal();
        //System.out.println("after delete " + key + " root=" + root);
        
        assert(check());
        
        assert(sz0 == (size() + 1));
    }

    // delete the key-value pair with the given key rooted at h
    private void delete(long h, long key, long[] output) { 
        
        //System.out.format("delete(%d, %d)\n", h, key);
        
        assert(nodeMap.containsKey(key));
        assert(nodeMap.containsKey(h));
        
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
            assert(nodeMap.leftIsSet(h));
            //h.left = delete(h.left, key);
            deleteFromLeftAssignLeft(h, h, key, output);
            output[0] = 0;
        } else {
            if (isLeftRed(h)) {
                h = rotateRight(h);
            }
            if (key == h && !nodeMap.rightIsSet(h)) {
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
                deleteMin(nodeMap.getRight(h), output3);
                        
                //System.out.println("deleteMin=" + Arrays.toString(output3));
        
                int hClr = nodeMap.getNodeColor(h);
                boolean hLeftExists = nodeMap.leftIsSet(h);
                boolean hRightExists = nodeMap.rightIsSet(h);
                boolean hParentExists = nodeMap.parentIsSet(h);
                long hLeft = nodeMap.getLeft(h);
                long hRight = nodeMap.getRight(h);
                long hParent = nodeMap.getParent(h);
                int hSize = nodeMap.getNodeSize(h);
               
                String hDBG = nodeToString(h);
                
                long x = output3[2];
               
                //System.out.println(
                //    "   x to get h fields except val. " + 
                //    "\n   x=" + nodeToString(x)
                //    + "\n   h=" + hDBG
                //);
               
                assert(x != h);
                
                //int xClr = nodeMap.getNodeColor(x);
                //boolean xLeftExists = nodeMap.leftIsSet(x);
                //boolean xRightExists = nodeMap.rightIsSet(x);
                boolean xParentExists = nodeMap.parentIsSet(x);
                //long xLeft = nodeMap.getLeft(x);
                //long xRight = nodeMap.getRight(x);
                long xParent = nodeMap.getParent(x);
                int xVal = nodeMap.getNodeValue(x);
                
                // finish deleting x if not already
                if (nodeMap.parentIsSet(x)) {
                    long minParent = nodeMap.getParent(x);
                    if (nodeMap.leftIsSet(minParent)
                        && nodeMap.getLeft(minParent) == x) {
                        nodeMap.unsetLeft(minParent);
                    } else if (nodeMap.rightIsSet(minParent)
                            && nodeMap.getRight(minParent) == x) {
                        nodeMap.unsetRight(minParent);
                    }
                    nodeMap.unsetParent(x);
                }
                if (nodeMap.containsKey(x)) {
                    nodeMap.remove(x);
                }
                
                if (hRightExists && hRight == x) {
                    hRightExists = false;
                }
                
                if (xParentExists) {
                    assert(!nodeMap.parentIsSet(x));
                }
                
                setHFromX(x, h, xVal, hClr, 
                    hLeftExists, hLeft, 
                    hRightExists, hRight, hParentExists, hParent, hSize - 1);
                    
                h = x;
                
                //System.out.println("  merged=" + nodeToString(h));
                
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
    
   /***************************************************************************
    *  Red-black tree helper functions.
    ***************************************************************************/

    /**
     * make a left-leaning link lean to the right.
     Note that the parent link logic is from Cormen et al. "Introduction to
     Algorithms".
     */
    protected long rotateRight(long h) {
        
        //System.out.println("  *RR " + h);
        
        // assert (h != null) && isRed(h.left);
        assert(nodeMap.containsKey(h));
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
        
        long x = nodeMap.getLeft(h);

        //System.out.println("  before RR h=" + nodeToString(h));
        //System.out.println("  before RR x=" + nodeToString(x));
        
        if (nodeMap.rightIsSet(x)) {
            //System.out.println("--0");
            //h.left = x.right;
            long xRight = nodeMap.getRight(x);
            
            nodeMap.updateLeft(h, xRight);
            nodeMap.updateParent(xRight, h);
        } else {
            //System.out.println("--1");
            nodeMap.unsetLeft(h);
        }
        if (nodeMap.parentIsSet(h)) {
            //System.out.println("--2");
            long hParent = nodeMap.getParent(h);
            // assign x as child of its new parent
            
            nodeMap.updateParent(x, hParent);
            if (nodeMap.leftIsSet(hParent) && nodeMap.getLeft(hParent) ==
                h) {
                //System.out.println("--3");
                nodeMap.updateLeft(hParent, x);
            } else {
                //System.out.println("--4");
                nodeMap.updateRight(hParent, x);
            }
        } else {
            //System.out.println("--5");
            nodeMap.unsetParent(x);
            //root = x;
        }
        
        //System.out.println("  in RR after h.left h=" + nodeToString(h));
        
        //x.right = h;
        nodeMap.updateRight(x, h);
        nodeMap.updateParent(h, x);
        
        //System.out.println("  in RR after x.right h=" + nodeToString(h));
        //System.out.println("  in RR after x.right x=" + nodeToString(x));
        
        //x.color = x.right.color;
        //x.right.color = RED;
        //x.size = h.size;
        nodeMap.updateNodeColor(x, nodeMap.getNodeColor(nodeMap.getRight(x)));
        nodeMap.updateNodeColor(nodeMap.getRight(x), RED);
        
        nodeMap.updateNodeSize(x, nodeMap.getNodeSize(h));
        
        //h.size = size(h.left) + size(h.right) + 1;
        int size = sizeLeft(h) + sizeRight(h) + 1;
        nodeMap.updateNodeSize(h, size);
        
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
        
        //System.out.println("  *RL h=" + h);
        
        // assert (h != null) && isRed(h.right);
        assert(nodeMap.containsKey(h));
        assert(isRightRed(h));
        
        //Node x = h.right;
        long x = nodeMap.getRight(h);
        
        //System.out.println("  before RL h=" + nodeToString(h));
        //System.out.println("  before RL x=" + nodeToString(x));

        //h.right = x.left;
        if (nodeMap.leftIsSet(x)) {
            //System.out.println("---0");
            long left = nodeMap.getLeft(x);
            
            nodeMap.updateRight(h, left);
            nodeMap.updateParent(left, h);
            
        } else {
            //System.out.println("---1");
            nodeMap.unsetRight(h);
        }
       
        //System.out.println("  in RL after h.right h=" + nodeToString(h));
        
        if (nodeMap.parentIsSet(h)) {
            //System.out.println("---2");
            long hParent = nodeMap.getParent(h);
            nodeMap.updateParent(x, hParent);
            if (nodeMap.leftIsSet(hParent) && nodeMap.getLeft(hParent) ==
                h) {
                //System.out.println("---3");
                nodeMap.updateLeft(hParent, x);
            } else {
                //System.out.println("---4");
                nodeMap.updateRight(hParent, x);
            }
        } else {
            //System.out.println("---5");
            nodeMap.unsetParent(x);
        }
                
        //x.left = h;
        nodeMap.updateLeft(x, h);
        nodeMap.updateParent(h, x);
        
        //System.out.println("  in RL after x.left h=" + nodeToString(h));
        //System.out.println("  in RL after x.left x=" + nodeToString(x));

        //x.color = x.left.color;
        nodeMap.updateNodeColor(x, nodeMap.getNodeColor(h));
        //x.left.color = RED;
        nodeMap.updateNodeColor(h, RED);
        
        nodeMap.updateNodeSize(x, nodeMap.getNodeSize(h));
        
        int size = sizeLeft(h) + sizeRight(h) + 1;
        nodeMap.updateNodeSize(h, size);
               
        //System.out.println("  after RL: h=" + nodeToString(h));
        
        return x;
    }

    // flip the colors of a node and its two children
    protected void flipColors(long h) {
        // h must have opposite color of its two children
        assert(nodeMap.containsKey(h) && nodeMap.leftIsSet(h) 
            && nodeMap.rightIsSet(h));
        assert(!isRed(h) &&  isLeftRed(h) &&  isRightRed(h))
        || (isRed(h)  && !isLeftRed(h) && !isRightRed(h));
        
        /*
        h.color = !h.color;
        h.left.color = !h.left.color;
        h.right.color = !h.right.color;
        */
        int clr = nodeMap.getNodeColor(h);
        clr ^= 1;
        nodeMap.updateNodeColor(h, clr);
        
        clr = nodeMap.getNodeColor(nodeMap.getLeft(h));
        clr ^= 1;
        nodeMap.updateNodeColor(nodeMap.getLeft(h), clr);
        
        clr = nodeMap.getNodeColor(nodeMap.getRight(h));
        clr ^= 1;
        nodeMap.updateNodeColor(nodeMap.getRight(h), clr);
    }

    // Assuming that h is red and both h.left and h.left.left
    // are black, make h.left or one of its children red.
    protected long moveRedLeft(long h) {
        assert(nodeMap.containsKey(h));
        assert(isRed(h) && !isLeftRed(h) && !isLeftLeftRed(h));

        //System.out.println("moveRedLeft " + nodeToString(h));
        
        flipColors(h);
        if (isRightLeftRed(h)) {
            //h.right = rotateRight(h.right);
            long rKey = rotateRight(nodeMap.getRight(h));
            nodeMap.updateRight(h, rKey);
            nodeMap.updateParent(rKey, h);
            h = rotateLeft(h);
            flipColors(h);
        }
        return h;
    }

    //move red node down the right spine of the tree
    // Assuming that h is red and both h.right and h.right.left
    // are black, make h.right or one of its children red.
    private long moveRedRight(long h) {
        assert (nodeMap.containsKey(h));
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
        
        assert (nodeMap.containsKey(h));

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
        nodeMap.updateNodeSize(h, size);
        
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
        if (!nodeMap.containsKey(x)) return -1;
        return 1 + Math.max(heightLeft(x), heightRight(x));
    }
    private int heightLeft(long x) {
        if (!nodeMap.leftIsSet(x)) return -1;
        return height(nodeMap.getLeft(x));
    }
    private int heightRight(long x) {
        if (!nodeMap.rightIsSet(x)) return -1;
        return height(nodeMap.getRight(x));
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
        //assert(nodeMap.containsKey(x));
        if (!nodeMap.containsKey(x)) {
            output[0] = -1;
            return;
        }
        while (nodeMap.leftIsSet(x)) {
            x = nodeMap.getLeft(x);
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
        //assert(nodeMap.containsKey(x));
        if (!nodeMap.containsKey(x)) {
            output[0] = -1;
            return;
        }
        while (nodeMap.rightIsSet(x)) {
            x = nodeMap.getRight(x);
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
        if (!nodeMap.containsKey(x)) {
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
            if (nodeMap.leftIsSet(x)) {
                floor(nodeMap.getLeft(x), key, output);
                return;
            }
            output[0] = -1;
            return;
        }
        if (nodeMap.rightIsSet(x)) {
            floor(nodeMap.getRight(x), key, output);
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
        long maxLower = Long.MAX_VALUE;
        while (nodeMap.containsKey(x) && key != x) {
            stack.add(x);
            //System.out.println("lower: x=" + x + " q=" + key);
            if (x < key) {
                if (maxLower == Long.MAX_VALUE) {
                    maxLower = x;
                } else if (x > maxLower) {
                    maxLower = x;
                }
            }
            if (key < x) {
                if (!nodeMap.leftIsSet(x)) {
                    break;
                }
                x = nodeMap.getLeft(x);
            } else {
                //System.out.println("   x=" + x.key);
                if (!nodeMap.rightIsSet(x)) {
                    break;
                }
                x = nodeMap.getRight(x);
            }
        }
        
        long y = -1;
        int yIdx = -1;
        boolean yIsSet = false;
        
        if (!nodeMap.containsKey(x)) {
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
        
        if (nodeMap.leftIsSet(x)) {
            max(nodeMap.getLeft(x), output);
            return;
        }
      
        //while (y != null && x == y.left) {
        while (yIsSet && 
            (
            (nodeMap.leftIsSet(y) && nodeMap.containsKey(x) &&
            nodeMap.getLeft(y) == nodeMap.getNodeValue(x)) 
            || 
            (!nodeMap.leftIsSet(y) && !nodeMap.containsKey(x))
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
            if (maxLower < Long.MAX_VALUE) {
                output[1] = maxLower;
                return;
            }
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
    public void higher(long key, final long[] output) {
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
        
        if (!nodeMap.containsKey(x)) {
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
            if (nodeMap.rightIsSet(x)) {
                ceiling(nodeMap.getRight(x), key, output);
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
        if (nodeMap.leftIsSet(x)) {
            ceiling(nodeMap.getLeft(x), key, output);
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
    private void higher(long x, long key, final long[] output) {  
        
        //System.out.println("higher: x=" + x + " key=" + key);
        
        /*
                    X
        left .lte.     right .gte.
        */
        
        while (nodeMap.containsKey(x)) {
            int cmp = key < x ? -1 : (key > x) ? 1 : 0; 
            //System.out.println("higher: x=" + x + " cmp=" + cmp);
            if (cmp < 0) {
                if (nodeMap.leftIsSet(x)) {
                    x = nodeMap.getLeft(x);
                } else {
                    break;
                }
            } else if (cmp > 0) {
                if (nodeMap.rightIsSet(x)) {
                    x = nodeMap.getRight(x);
                } else {
                    break;
                }
            } else {
                break;
            }
        }
                
        // right node has larger key
        if (nodeMap.rightIsSet(x)) {
            output[0] = 0;
            min(nodeMap.getRight(x), output);
            //System.out.println("min=" + Arrays.toString(output));
            return;
        }
        boolean yIsSet = false;
        long y = -1;
        if (nodeMap.parentIsSet(x)) {
            y = nodeMap.getParent(x);
            yIsSet = true;
            //System.out.println("y=" + y);
        }
        while (yIsSet && nodeMap.rightIsSet(y) &&
            x == nodeMap.getRight(y)) {
            //System.out.println("y=" + y + " setting x=" + x + " to y");
            x = y;
            if (nodeMap.parentIsSet(x)) {
                y = nodeMap.getParent(x);
            } else {
                yIsSet = false;
                y = -1;
            }
        }
        if (x > key && yIsSet && y > key) {
            output[0] = 0;
            if (x < y) {
                output[1] = x;
            } else {
                output[1] = y;
            }
        } else if (x > key) {
            output[0] = 0;
            output[1] = x;
        } else if (yIsSet && y > key) {
            output[0] = 0;
            output[1] = y;
        } else {
            output[0] = -1;
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
        
        assert(nodeMap.containsKey(x));
        assert(k >= 0 && k < size(x));
        
        output[0] = 0;
        
        int t = sizeLeft(x); 
        if (t > k) {
            select(nodeMap.getLeft(x),  k, output);
        } else if (t < k) {
            if (nodeMap.rightIsSet(x)) {
                select(nodeMap.getRight(x), k-t-1, output);
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
        
        if (!nodeMap.containsKey(x)) {
            output[1] = 0;
            return;
        } 
        int cmp = key < x ? -1 : (key > x) ? 1 : 0;  
        if (cmp < 0) {
            if (nodeMap.leftIsSet(x)) {
                rank(key, nodeMap.getLeft(x), output);
            } else {
                output[1] = 0;
            }
        } else if (cmp > 0) {
            //1 + size(x.left) + rank(key, x.right);
            int sz;
            if (nodeMap.rightIsSet(x)) {
                output[0] = 0;
                rank(key, nodeMap.getRight(x), output);
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
        
        if (!nodeMap.containsKey(x)) {
            return;
        } 
        int cmplo = lo < x ? -1 : (lo > x) ? 1 : 0;
        int cmphi = hi < x ? -1 : (hi > x) ? 1 : 0; 
        if (cmplo < 0) {
            if (nodeMap.leftIsSet(x)) {
                keys(nodeMap.getLeft(x), queue, lo, hi);
            }
        } 
        if (cmplo <= 0 && cmphi >= 0) {
            queue.add(x);
        } 
        if (cmphi > 0) {
            if (nodeMap.rightIsSet(x)) {
                keys(nodeMap.getRight(x), queue, lo, hi);
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
        if (!nodeMap.containsKey(x)) {
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
        if (!nodeMap.leftIsSet(x)) {
            return true;
        }
        long key = nodeMap.getLeft(x);
        return isBST(key, min, max);
    }
    private boolean isBSTRight(long x, Long min, Long max) {
        if (!nodeMap.rightIsSet(x)) {
            return true;
        }
        long key = nodeMap.getRight(x);
        return isBST(key, min, max);
    }

    private boolean isSizeConsistent() { return isSizeConsistent(root); }
    private boolean isSizeConsistent(long x) {
        if (!nodeMap.containsKey(x)) return true;
        if (nodeMap.getNodeSize(x) != sizeLeft(x) + sizeRight(x) + 1) return false;
        return isSizeConsistentLeft(x) && isSizeConsistentRight(x);
    }
    private boolean isSizeConsistentLeft(long x) {
        if (!nodeMap.leftIsSet(x)) return true;
        long key = nodeMap.getLeft(x);
        return isSizeConsistent(key);
    }
    private boolean isSizeConsistentRight(long x) {
        if (!nodeMap.rightIsSet(x)) return true;
        long key = nodeMap.getRight(x);
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
        if (!nodeMap.containsKey(x)) {
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
        if (!nodeMap.leftIsSet(x)) return true;
        return is23(nodeMap.getLeft(x));
    }
    private boolean is23Right(long x) {
        if (!nodeMap.rightIsSet(x)) return true;
        return is23(nodeMap.getRight(x));
    }

    // do all paths from root to leaf have same number of black edges?
    private boolean isBalanced() { 
        if (!rootIsSet) {
            return true;
        }
        int black = 0;     // number of black links on path from root to min
        long x = root;
        while (nodeMap.containsKey(x)) {
            if (!isRed(x)) {
                black++;
            }
            if (!nodeMap.leftIsSet(x)) {
                break;
            }
            x = nodeMap.getLeft(x);
        }
        return isBalanced(root, black);
    }

    // does every path from the root to a leaf have the given number of black links?
    private boolean isBalanced(long x, int black) {
        if (!nodeMap.containsKey(x)) {
            return black == 0;
        }
        if (!isRed(x)) {
            black--;
        }
        return isLeftBalanced(x, black) && isRightBalanced(x, black);
    }
    private boolean isLeftBalanced(long x, int black) {
        if (!nodeMap.leftIsSet(x)) {
            return black == 0;
        }
        long key = nodeMap.getLeft(x);
        return isBalanced(key, black);
    } 
    private boolean isRightBalanced(long x, int black) {
        if (!nodeMap.rightIsSet(x)) {
            return black == 0;
        }
        long key = nodeMap.getRight(x);
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
                
                node = nodeMap.leftIsSet(node.longValue()) ?
                    nodeMap.getLeft(node.longValue()) : null;
            
            } else if (count < sz) {
                
                node = stack.pop();
                
                array[count] = node;
                count++;
                
                //System.out.println(node.key);
                
                node = nodeMap.rightIsSet(node.longValue()) ?
                    nodeMap.getRight(node.longValue()) : null;
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
                
                node = nodeMap.leftIsSet(node.longValue()) ?
                    nodeMap.getLeft(node.longValue()) : null;
            
            } else if (count < sz) {
                
                node = stack.pop();
                
                node = nodeMap.rightIsSet(node.longValue()) ?
                    nodeMap.getRight(node.longValue()) : null;
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
            
            if (nodeMap.leftIsSet(node.longValue())) {
                stack.push(nodeMap.getLeft(node.longValue()));
            }

            if (nodeMap.rightIsSet(node.longValue())) {
                stack.push(nodeMap.getRight(node.longValue()));
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
        
        total += NodeMap.estimateSizeOnHeap(numberOfEntries);
             
        return total;
    }
    
    private String nodeToString(long key) {
        //assert(nodeMap.containsKey(key));
        if (!nodeMap.containsKey(key)) {
            // this can happen in the middle of a method, for example,
            // when root has been removed and new is not yet assigned
            System.out.println("ERROR: key " + key + " not in maps");
            return "";
        }
        
        StringBuilder sb = new StringBuilder();
        //node=key=0 val=0 color=false size=1
        sb.append("key=").append(key).append(" val=").append(nodeMap.getNodeValue(key));
        sb.append(" color=").append(nodeMap.getNodeColor(key));
        sb.append(" size=").append(nodeMap.getNodeSize(key));
        sb.append(" p=");
        if (nodeMap.parentIsSet(key)) {
           sb.append(nodeMap.getParent(key));
        }
        sb.append(" l=");
        if (nodeMap.leftIsSet(key)) {
           sb.append(nodeMap.getLeft(key));
        }
        sb.append(" r=");
        if (nodeMap.rightIsSet(key)) {
           sb.append(nodeMap.getRight(key));
        }
        return sb.toString();
    }

    private boolean isParentChildConsistent() {
        
        boolean passed = true;
        
        if (!rootIsSet) {
            if (!nodeMap.isEmpty()) {
                System.err.println("maps not empty, but root is");
                passed = false;
            }
        }
        
        //System.out.println("root=" + nodeToString(root));
        
        if (nodeMap.parentIsSet(root)) {
            System.err.println("root should not have parent key");
            passed = false;
        }
        
        long[] nodes = getPreOrderTraversalIterative(root, 0);
        for (long key : nodes) {
            if (nodeMap.leftIsSet(key)) {
                long child = nodeMap.getLeft(key);
                if (!nodeMap.parentIsSet(child)) {
                    System.err.format(
                        "error in %d.left=%d has no parent\n",
                        key, child);
                    passed = false;
                } else if (nodeMap.getParent(child) != key) {
                    System.err.format(
                        "error in left: %d.left=%d but %d.parent=%d\n",
                        key, child, child, nodeMap.getParent(child));
                    passed = false;
                }
            }
            if (nodeMap.rightIsSet(key)) {
                long child = nodeMap.getRight(key);
                if (!nodeMap.parentIsSet(child)) {
                    System.err.format(
                        "error in %d.right=%d has no parent\n",
                        key, child);
                    passed = false;
                } else if (nodeMap.getParent(child) != key) {
                    System.err.format(
                        "error in right: %d.right=%d but %d.parent=%d\n",
                        key, child, child, nodeMap.getParent(child));
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
        
        int n = nodeMap.getNodeSize(node);
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
                if (!nodeMap.containsKey(z)) {
                    continue;
                }
                int clr = nodeMap.getNodeColor(z);
                String keyC = Long.toString(z);
                String clrC = Integer.toString(clr);
            
                sb0.append(keyC);
                sb1.append(clrC);
                addSpaces(1, sb0);
                addSpaces(w + 1 - clrC.length(), sb1);
                
                if (nodeMap.leftIsSet(z)) {
                    nextLevelQ.add(nodeMap.getLeft(z));
                } else {
                    nextLevelQ.add(Long.MIN_VALUE);
                }
                if (nodeMap.rightIsSet(z)) {
                    nextLevelQ.add(nodeMap.getRight(z));
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
        return nodeMap.leftIsSet(key);
    }
    private boolean containsRight(long key) {
        return nodeMap.rightIsSet(key);
    }
    private boolean containsParent(long key) {
        return nodeMap.parentIsSet(key);
    }
    private boolean containsParentLeft(long key) {
        if (containsParent(key)) {
            return containsLeft(nodeMap.getParent(key));
        }
        return false;
    }
    private boolean containsParentRight(long key) {
        if (containsParent(key)) {
            return containsRight(nodeMap.getParent(key));
        }
        return false;
    }

    private void deleteFromLeftAssignLeft(long delH, long asnH, long key, 
        long[] output) {
       
        delete(nodeMap.getLeft(delH), key, output);
        if (output[0] == -1) {
            nodeMap.unsetLeft(asnH);
        } else {
            nodeMap.updateLeft(delH, output[1]);
            nodeMap.updateParent(output[1], delH);
        }
    }

    private void deleteFromRightAssignRight(long hDel, long hAsn, long key, 
        long[] output) {
        
        assert(nodeMap.rightIsSet(hDel));
        
        //h.right = delete(h.right, key);

        delete(nodeMap.getRight(hDel), key, output);
        
        if (output[0] == -1) {
            nodeMap.unsetRight(hAsn);
        } else {
            nodeMap.updateRight(hAsn, output[1]);
            nodeMap.updateParent(output[1], hAsn);
        }
    }
             
    private void setHFromX(long x, long h, int xVal, int hClr, 
        boolean hLeftExists, long hLeft, 
        boolean hRightExists, long hRight, 
        boolean hParentExists, long hParent, int hSize) {
        
        assert(x != h);
        nodeMap.put(x, xVal, hClr, hSize);
     
        if (hParentExists) {
            nodeMap.updateParent(x, hParent);
            if (nodeMap.leftIsSet(hParent) 
                && nodeMap.getLeft(hParent) == h) {
                nodeMap.updateLeft(hParent, x);
            } else {
                assert(nodeMap.rightIsSet(hParent) 
                    && nodeMap.getRight(hParent) == h);
                nodeMap.updateRight(hParent, x);
            }
        } else {
            //h is root
        }
        
        if (hLeftExists) {
            nodeMap.updateParent(hLeft, x);
            nodeMap.updateLeft(x, hLeft);
        }
        if (hRightExists) {
            nodeMap.updateParent(hRight, x);
            nodeMap.updateRight(x, hRight);
        }
        
        nodeMap.remove(h); 
        
        if (h == root) {
            root = x;
        }
    }
    
}
