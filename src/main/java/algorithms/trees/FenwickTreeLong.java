package algorithms.trees;

import java.util.Arrays;

/**
 * Fenwick Tree a.k.a. Binary Indexed Tree (BIT).
 * an update-able datastructure which can be used to hold numbers, given an index for them
 * that starts at 1 (or 0 w/ flag).  It's an efficient holder of prefix sums that is update-able in a
 * runtime complexity of O(log(n)) where n is the number of values it is constructed with.
 *
 <pre>
 references:
 adapted from the YouTube lecture of William Fiset (add reference)
 and his source code:
 https://github.com/williamfiset/Algorithms/tree/master/src/main/java/com/williamfiset/algorithms/datastructures/fenwicktree
 The William Fiset code uses the MIT license:
 https://github.com/williamfiset/Algorithms/blob/master/LICENSE
 </pre>

 internally, the tree is built to store prefix sums over varying ranges of indexes of a.
 that is, the positions of partial sums of a are dependent upon the bits set in index i.
 NOTE: the actual tree positions are i+1 for a 0-based indexing, but the math is
 explained here for 1-based indexing:

 build tree:
    tree = copy of a.
    then loop over tree indices 1 thru n-1
       tree[1 + lsb(i)] += tree[i]
 can see this for the loop:
    i    LSB(i)  LSB(i)_bit    # of #s summed    j
    1     1        0             1               i+1
    2     2        1             2               i+2
    3     1        0             1               i+1
    4     4        2             4               i+4
    5     1        0             1               i+1
    6     2        1             2               i+2
    7     1        0             1               i+1
    8     8        3             8               i+8
 ...

 then reading from the tree:
    get(i,j) = prefixSum(j) - prefixSum(i - 1);
    where prefixSum(idx) = {
        long sum = 0L;
        while (idx != 0) {
           sum += tree[idx];
           idx -= lSB(idx);
        }
        return sum;
    }

 updating tree:
     update(idx, val) = {
         while (idx < tree.length) {
             tree[idx] += val;
             idx += lSB(_idx);
         }
     }

 <pre>
 for 2D Fenwick tree, see end of article:
 https://www.topcoder.com/thrive/articles/Binary%20Indexed%20Trees
 </pre>
 */
public class FenwickTreeLong {

    // Note, the internal arrays are represented as bit structures.
    // the first bit, 0, is not used.

    final protected long[] tree;
    final boolean use0Based;

    /**
     * constuct a fenwick tree to hold n values.
     * @param n maximum number of numbers to be placed in the tree.
     @param use0Based if true, you are inserting all values in the array and
      *                  access this FenwickTree using indexes
      *                  0 through n-1, inclusive where n is the length of the array,
      *                  else the values array is taken to have an unused 0 as the
      *                  first element and you will access this FewnwickTree using
      *                  indexes 1 through n, inclusive.
     */
    public FenwickTreeLong(int n, boolean use0Based) {
        this.tree = new long[n+1];
        this.use0Based = use0Based;
    }

    /**
     * construct a fenwick tree for the given values.
     * runtime complexity is O(n) where n is values.length.  the space complexity is O(n).
     * @param values the values to insert into the FenwickTree.
     * @param use0Based if true, you are inserting all values in the array and
     *                  access this FenwickTree using indexes
     *                  0 through n-1, inclusive where n is the length of the array,
     *                  else the values array is taken to have an unused 0 as the
     *                  first element and you will access this FewnwickTree using
     *                  indexes 1 through n, inclusive.
     */
    public FenwickTreeLong(long[] values, boolean use0Based) {

        if (!use0Based && (values[0] != 0L)) {
            throw new IllegalArgumentException("when use0Based is false, a 0 is expected " +
                    "as the first element in values");
        }
        this.use0Based = use0Based;
        if (use0Based) {
            this.tree = new long[values.length + 1];
            System.arraycopy(values, 0, tree, 1, values.length);
        } else {
            this.tree = Arrays.copyOf(values, values.length);
        }

        // add item to immediate cell above
        int j;
        for (int i = 1; i < tree.length; ++i) {
            j = i + lSB(i);
            if (j < tree.length) {
                tree[j] += tree[i];
            }
        }
    }

    protected int lSB(int idx) {
        return idx & -idx;
    }

    protected int getIndex(int idx) {
        if (use0Based) return idx+1;
        return idx;
    }

    private long prefixSum(final int i) {
        int idx = use0Based ? i+1: i;
        long sum = 0L;
        while (idx != 0) {
            sum += tree[idx];
            // these are all equiv.  have chosen the one complementary to build and add methods
            // drop the LSB  A = A & (A - 1)
            //idx &= (idx - 1);
            //idx &= ~lSB(idx);
            idx -= lSB(idx);
        }
        return sum;
    }


    // r.t.c. O(log(n))
    // Returns the sum of the interval [leftIdx, rightIdx],
    public long sum(int leftIdx, int rightIdx) {
        if (rightIdx < leftIdx) throw new IllegalArgumentException("Make sure rightIdx >= leftIdx");
        return prefixSum(rightIdx) - prefixSum(leftIdx - 1);
    }

    /**
     * get the value of the tree element at index idx
     * @param idx tree index
     * @return
     */
    public long get(int idx) {
        return sum(idx, idx);
    }

    /**
     * add v to the element at index idx in the tree.
     * r.t.c. is O(log(n))
     * @param idx tree index
     * @param v value to add to value at tree index idx
     */
    public void add(final int idx, long v) {
        int _idx = use0Based ? idx+1: idx;
        while (_idx < tree.length) {
            tree[_idx] += v;
            _idx += lSB(_idx);
        }
    }

    /**
     * et index idx to be equal to v.
     * r.t.c. O(log(n))
     * @param idx index of tree to set value v to
     * @param v value to set tree element to
     */
    public void set(int idx, long v) {
        add(idx, v - sum(idx, idx));
    }

    /*
    min in a range
    max in a range
    n elements in range >= a value
     */

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("use 0 based indexes = %b", use0Based));
        sb.append(String.format("/ntree=%s\n",java.util.Arrays.toString(tree)));
        return sb.toString();
    }

}
