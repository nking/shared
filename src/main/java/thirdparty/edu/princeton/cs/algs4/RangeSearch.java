package thirdparty.edu.princeton.cs.algs4;

/******************************************************************************
adapted from http://algs4.cs.princeton.edu/92search/
    copyright for authors Robert Sedgewick and Kevin Wayne
    is GPLV3, http://algs4.cs.princeton.edu/faq/ 
*
* Compilation:  javac RangeSearch.java
 *  Execution:    java RangeSearch left-pipe words.txt
 *  
 *  Range search implemented using a randomized BST.
 
 
 *
 @param <Key> parameter type of key
 @param <Value>  parameter type of value
 ******************************************************************************/

public class RangeSearch<Key extends Comparable<Key>, Value>  {

    /**
     *
     */
    protected RangeSearchNode<Key, Value> root;   // root of the BST
    
    //BST helper node data type

    /**
     *
     @param <T> parameter type of key
     @param <S> parameter type of value
     */
    protected class RangeSearchNode<T, S> {
        T key;              // key
        S val;              // associated data
        RangeSearchNode<T, S> left, right;   // left and right subtrees
        int N;              // node count of descendents

        /**
         *
         @param key
         @param val
         */
        public RangeSearchNode(T key, S val) {
            this.key = key;
            this.val = val;
            this.N = 1;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder("key=");
            sb.append(key).append(" n=").append(N)
                .append(" left=").append(left)
                .append(" right=").append(right)
                .append(" val=").append(val);
            return sb.toString();
        }
    }
    
   /***************************************************************************
    *  BST search
     @param key
     @return 
    ***************************************************************************/

    public boolean contains(Key key) {
        return (get(key) != null);
    }
    
    /**
     * runtime complexity at worst is O(lg2(N)).
     * 
     @param key
     @return value associated with the given key
               if no such value, return null
    */
    public Value get(Key key) {
        return get(root, key);
    }

    /**
     * runtime complexity at worst is O(lg2(N)).
     * 
     @param x
     @param key
     @return 
     */
    private Value get(RangeSearchNode<Key, Value> x, Key key) {
        if (x == null) return null;
        int cmp = key.compareTo(x.key);
        if      (cmp == 0) return x.val;
        else if (cmp  < 0) return get(x.left,  key);
        else               return get(x.right, key);
    }

   /***************************************************************************
    *  randomized insertion
    ***************************************************************************/
    /**
     * insert interval, but if it intersects with an
     * interval already in tree, set the existing value
     * to the given val, and return the previous value
     * before the overwrite by val.  If null is returned
     * from this method, then the insert succeeded,
     * that is, there were no collisions.
     * 
     * runtime complexity at worst is O(lg2(N)).
     * 
     @param key
     @param val
     @return the value that was replaced with given val
     * if the Key intersected with another, preventing 
     * an insert of Key, but updating existing with val.
     * Note that the return is null when the insert
     * succeeded.
     */
    @SuppressWarnings({"unchecked"})
    public Value put(Key key, Value val) {
        
        //to return whether insert was successful
        Object[] replaced = new Object[1];
        
        root = put(root, key, val, replaced);
        
        return (replaced[0] == null) ? null : (Value)replaced[0];
        
        //System.out.println("<==root=" + root);
    }
    
    /**
     * put key in map, and if compareVal is greater than or equal to the
     * replaced value, re-insert the replaced value and return false, else
     * the insert succeeded and returns true.
     * 
     * runtime complexity at worst is O(lg2(N)).
     * 
     @param <Value2>
     @param key
     @param val
     @param compareVal
     @return 
     */
    public <Value2 extends Comparable<Value>> boolean 
        putIfLessThan(Key key, Value val, Value2 compareVal) {
          
        //if the insert replaced an object, this holds the value, then key
        Object[] replaced = new Object[1];
        boolean[] inserted = new boolean[1];
        
        root = putIfLessThan(root, key, val, compareVal, replaced, inserted);
        
        return inserted[0];
        
        //System.out.println("<==root=" + root);
    }
    
    /**
     * make new node the root with uniform probability.
     * 
     * runtime complexity at worst is O(lg2(N)).
     *
     @param <Value2>
     @param x
     @param key
     @param val
     @param compareVal
     @param replaced
     @param inserted
     @return 
     */
    @SuppressWarnings({"unchecked"})
    private <Value2 extends Comparable<Value>> 
        RangeSearchNode<Key, Value> 
        putIfLessThan(RangeSearchNode<Key, Value> x, 
        Key key, Value val, Value2 compareVal, Object[] replaced,
        boolean[] inserted) {
                    
        if (x == null) {
            inserted[0] = true;
            return new RangeSearchNode<Key, Value>(key, val);
        }
        
        int cmp = key.compareTo(x.key);
        if (cmp == 0) {
            int cmpV = compareVal.compareTo(val);
            if (cmpV < 0) {
                // continue w/ insert
                replaced[0] = x.val;
                x.val = val;
                inserted[0] = true;
            } else {
                inserted[0] = false;
            }
            return x;
        }
        if (StdRandom.bernoulli(1.0 / (size(x) + 1.0))) {
            return putRootIfLessThan(x, key, val, compareVal, replaced, inserted);
        }
        if (cmp < 0) {
            x.left  = putIfLessThan(x.left,  key, val, compareVal, replaced,
                inserted);
        } else {
            x.right = putIfLessThan(x.right, key, val, compareVal, replaced,
                inserted);
        }
        // (x.N)++;
        fix(x);
        return x;
    }
    
    /**
     * make new node the root with uniform probability.
     * 
     * runtime complexity at worst is O(lg2(N)).
     * 
     @param x
     @param key
     @param val
     @param replaced
     @return 
     */
    private RangeSearchNode<Key, Value> put(RangeSearchNode<Key, Value> x, 
        Key key, Value val, Object[] replaced) {
    
        if (x == null) return new RangeSearchNode<Key, Value>(key, val);
        
        int cmp = key.compareTo(x.key);
        if (cmp == 0) {
            replaced[0] = x.val;
            x.val = val;
            return x;
        }
        if (StdRandom.bernoulli(1.0 / (size(x) + 1.0))) {
            return putRoot(x, key, val, replaced);
        }
        if (cmp < 0) {
            x.left  = put(x.left,  key, val, replaced);
        } else {
            x.right = put(x.right, key, val, replaced);
        } 
        // (x.N)++;
        fix(x);
        return x;
    }

    private <Value2 extends Comparable<Value>> 
    RangeSearchNode<Key, Value> putRootIfLessThan(
        RangeSearchNode<Key, Value> x, Key key, Value val, Value2 compareVal,
        Object[] replaced, boolean[] inserted) {
        
        if (x == null) {
            inserted[0] = true;
            return new RangeSearchNode<Key, Value>(key, val);
        }
        
        int cmp = key.compareTo(x.key);
        if (cmp == 0) {
            int cmpV = compareVal.compareTo(val);
            if (cmpV < 0) {
                replaced[0] = x.val;
                x.val = val;
                inserted[0] = true;
            } else {
                inserted[0] = false;
            }
            return x; 
        } else if (cmp  < 0) { 
            x.left  = putRootIfLessThan(x.left,  key, val, compareVal, replaced,
                inserted); 
            x = rotR(x); 
        } else { 
            x.right = putRootIfLessThan(x.right, key, val, compareVal, replaced,
                inserted); 
            x = rotL(x); 
        }
        return x;
    }
    
    /**
     * runtime complexity at most is O(lg2(N)).
     * 
     @param x
     @param key
     @param val
     @param replaced
     @return 
     */
    private RangeSearchNode<Key, Value> putRoot(
        RangeSearchNode<Key, Value> x, Key key, Value val,
        Object[] replaced) {
        
        if (x == null) return new RangeSearchNode<Key, Value>(key, val);
        int cmp = key.compareTo(x.key);
        if (cmp == 0) {
            replaced[0] = x.val;
            x.val = val;
            return x; 
        } else if (cmp  < 0) { 
            x.left  = putRoot(x.left,  key, val, replaced); 
            x = rotR(x); 
        } else { 
            x.right = putRoot(x.right, key, val, replaced); 
            x = rotL(x); 
        }
        return x;
    }

   /***************************************************************************
    *  deletion
    
    * runtime complexity worse case is O(lg2(N)).
     @param a
     @param b
     @return 
    ***************************************************************************/
    private RangeSearchNode<Key, Value> joinLR(RangeSearchNode<Key, Value> a, 
        RangeSearchNode<Key, Value> b) { 
   
        if (a == null) return b;
        if (b == null) return a;

        // O(1), but statement requires use of random algorithm.
        if (StdRandom.bernoulli((double) size(a) / (size(a) + size(b))))  {
            a.right = joinLR(a.right, b);
            fix(a);
            return a;
        } else {
            b.left = joinLR(a, b.left);
            fix(b);
            return b;
        }
    }

    /**
     * runtime complexity at worst is O(lg2(N)).
     * 
     @param x
     @param key
     @return 
     */
    private RangeSearchNode<Key, Value> remove(RangeSearchNode<Key, Value> x, Key key) {
        if (x == null) return null; 
        int cmp = key.compareTo(x.key);
        if      (cmp == 0) x = joinLR(x.left, x.right);
        else if (cmp  < 0) x.left  = remove(x.left,  key);
        else               x.right = remove(x.right, key);
        fix(x);
        return x;
    }

    /**
     * runtime complexity at worse is O(lg2(N)).
     * 
     * remove and return value associated with given key; if no such key, 
     * return null
     * 
     @param key
     @return 
     */
    public Value remove(Key key) {
        Value val = get(key);
        root = remove(root, key);
        return val;
    }

   /***************************************************************************
    *  Range searching
    ***************************************************************************/

    /**
     * return all keys in given interval
     * 
     * runtime complexity is O(1) to O(lg2(N)).
     *
     @param min
     @param max
     @return 
     */
    public Iterable<Key> range(Key min, Key max) {
        return range(new Interval<Key>(min, max));
    }
    /**
     * return all keys in given interval
     * 
     * runtime complexity is O(1) to O(lg2(N)).
     * 
     @param interval
     @return 
     */
    public Iterable<Key> range(Interval<Key> interval) { 
        Queue<Key> list = new Queue<Key>();
        range(root, interval, list);
        return list;
    }

    /**
     * 
     * runtime complexity is O(1) to O(lg2(N)).
     * 
     @param x
     @param interval
     @param list 
     */
    private void range(RangeSearchNode<Key, Value> x, Interval<Key> interval, 
        Queue<Key> list) {
        if (x == null) return;
        if (!less(x.key, interval.min()))  range(x.left, interval, list);
        if (interval.contains(x.key))      list.enqueue(x.key);
        if (!less(interval.max(), x.key))  range(x.right, interval, list);
    }

    /***************************************************************************
    *  Utility functions
    ***************************************************************************/
    /**
     *
     * runtime complexity is O(1)
     * 
     @return the smallest key
     */
    public Key min() {
        Key key = null;
        for (RangeSearchNode<Key, Value> x = root; x != null; x = x.left)
            key = x.key;
        return key;
    }
    
    /**
     *
     * runtime complexity is O(1)
     * 
     @return the largest key
     */
    public Key max() {
        Key key = null;
        for (RangeSearchNode<Key, Value> x = root; x != null; x = x.right)
            key = x.key;
        return key;
    }

   /***************************************************************************
    *  useful binary tree functions
    ***************************************************************************/ 
    /**
     * runtime complexity is O(1)
     * 
     @return number of nodes in subtree rooted at x
     */
    public int size() { return size(root); }
    private int size(RangeSearchNode<Key, Value> x) { 
        if (x == null) return 0;
        else           return x.N;
    }

    /**
     * height of tree (empty tree height = 0)
     * 
     * runtime complexity is O(lg_2(N))
     @return 
     */
    public int height() { return height(root); }
    private int height(RangeSearchNode<Key, Value> x) {
        if (x == null) return 0;
        return 1 + Math.max(height(x.left), height(x.right));
    }


   /***************************************************************************
    *  helper BST functions
    ***************************************************************************/

    /**
     * fix subtree count field.
     * 
     * runtime complexity O(1)
     @param x 
     */
    private void fix(RangeSearchNode<Key, Value> x) {
        if (x == null) return;                 // check needed for remove
        x.N = 1 + size(x.left) + size(x.right);
    }

    /**
     * right rotate
     * 
     * runtime complexity O(1)
     @param h
     @return 
     */
    private RangeSearchNode<Key, Value> rotR(RangeSearchNode<Key, Value> h) {
        RangeSearchNode<Key, Value> x = h.left;
        h.left = x.right;
        x.right = h;
        fix(h);
        fix(x);
        return x;
        /*      X                       H
                0                      10
             1     2               11     12
        
            
                 11                     
             1        10  
                    2  `12
        */
    }

    /**
     * left rotate
     * 
     * runtime complexity is O(1).
     * 
     @param h
     @return 
     */
    private RangeSearchNode<Key, Value> rotL(RangeSearchNode<Key, Value> h) {
        
        RangeSearchNode<Key, Value> x = h.right;
        h.right = x.left;
        x.left = h;
        
        fix(h);
        fix(x);
        
        return x;
    }


   /***************************************************************************
    *  Debugging functions that test the integrity of the tree
     @return 
    ***************************************************************************/

    // check integrity of subtree count fields
    public boolean check() { return checkCount() && isBST(); }

    // check integrity of count fields
    private boolean checkCount() { return checkCount(root); }
    private boolean checkCount(RangeSearchNode<Key, Value> x) {
        if (x == null) return true;
        return checkCount(x.left) && checkCount(x.right) && (x.N == 1 + size(x.left) + size(x.right));
    }


    // does this tree satisfy the BST property?
    private boolean isBST() { return isBST(root, min(), max()); }

    // are all the values in the BST rooted at x between min and max, and recursively?
    //  runtime complexity worse case is O(lg_2(N))
    private boolean isBST(RangeSearchNode<Key, Value> x, Key min, Key max) {
        if (x == null) return true;
        if (less(x.key, min) || less(max, x.key)) return false;
        return isBST(x.left, min, x.key) && isBST(x.right, x.key, max);
    } 



   /***************************************************************************
    *  helper comparison functions
     @param k1
     @param k2
     @return 
    ***************************************************************************/

    private boolean less(Key k1, Key k2) {
        return k1.compareTo(k2) < 0;
    }


   /***************************************************************************
    *  test client
     @param args
    ***************************************************************************/
    public static void main(String[] args) {
        /*
        int N = 0;
        RangeSearch<String, Integer> st = new RangeSearch<String, Integer>();
        while (!StdIn.isEmpty()) {
            String s = StdIn.readString();
            st.put(s, N++);
        }

        StdOut.println("height:          " + st.height());
        StdOut.println("size:            " + st.size());
        StdOut.println("min key:         " + st.min());
        StdOut.println("max key:         " + st.max());
        StdOut.println("integrity check: " + st.check());
        StdOut.println();

        StdOut.println(new Interval<String>("kevin", "kfg"));
        Iterable<String> list = st.range(new Interval<String>("kevin", "kfg"));
        for (String s : list)
            StdOut.println(s + " " + st.get(s));
        StdOut.println();

        StdOut.println(new Interval<String>("paste", "pasty"));
        list = st.range(new Interval<String>("paste", "pasty"));
        for (String s : list)
            StdOut.println(s);
        StdOut.println();
        */
    }

}
