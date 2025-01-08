package algorithms.trees;

import java.util.Arrays;

/**
 * an update-able datastructure which can be used to hold numbers, given an index for them
 * that starts at 1.  It's an efficient holder of prefix sums that is update-able in a
 * runtime complexity of O(log(n)) where n is the number of values it is constructed with.
 *
 * If one needs range minimum or maximum queries in an update-able structure for
 * r.t.c. of O(log(n)), should probably use a sortedtree like TreeSet which is a Red-Blac Tree.
 *
 <pre>
 references:
 adapted from the YouTube lecture of William Fiset (add reference)
 and his source code:
 https://github.com/williamfiset/Algorithms/tree/master/src/main/java/com/williamfiset/algorithms/datastructures/fenwicktree
 The William Fiset code uses the MIT license:
 https://github.com/williamfiset/Algorithms/blob/master/LICENSE
 </pre>
 */
public class FenwickTreeLong {

    // Note, the internal arrays are represented as bit structures.
    // the first bit, 0, is not used.

    final protected long[] tree;
    final boolean use0Based;

    /**
     * constuct a fenwick tree to hold n values.
     * @param n
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
            if (j < tree.length) tree[j] += tree[i];
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
            // drop the LSB  A = A & (A - 1)
            idx &= (idx - 1);
            //idx &= ~lSB(idx); // Equivalently, i -= lsb(i);
        }
        return sum;
    }


    // Returns the sum of the interval [leftIdx, rightIdx], O(log(n))
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
