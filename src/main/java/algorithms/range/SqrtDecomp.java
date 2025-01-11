package algorithms.range;

import java.util.Arrays;

/**
 * a data structure that takes an integer array of data and performs sum and update
 * operations on it.  
 * The update operation has r.t.c. O(sqrt(n)) which is better
 * performance than the O(n) update for a prefix sum array.
 * If no updates are ever performed, one should prefer a prefix sum array because of the O(1)
 * sum operation.
 * <pre>
 *     reference https://usaco.guide/plat/sqrt?lang=java#blocking
 * </pre>
 */
public class SqrtDecomp {
    
    protected int blockSize;//nBins
    protected int n;
    protected long[] blocks;
    protected long[] a;

    public SqrtDecomp(long[] a) {
        n = a.length;
        blockSize = (int)Math.sqrt(n) + 1;
        blocks = new long[blockSize];
        for (int i = 0; i < n; ++i) {
            blocks[i/blockSize] += a[i];
        }
        this.a = Arrays.copyOf(a, n);
    }

    /**
     * r.t.c. O(queries.length * sqrt(n)) where n = a.length;
     * @param queries array of elements containing ranges, e.g. [l, r] inclusive
     * @param indexesAre0Based whether the query range points are 0-based indexes
     *                         (else false is 1-based indexing).
     * @return
     */
    public long[] sum(int[][] queries, boolean indexesAre0Based) {
        long[] out = new long[queries.length];
        int r, l;
        for (int i = 0; i < queries.length; ++i) {
            out[i] = sum(queries[i][0], queries[i][1], indexesAre0Based);
        }
        return out;
    }

    /**
     * r.t.c. O(sqrt(n)) where n = a.length;
     * @param l start index of query range
     * @param r stop index of query range inclusive
     * @param indexesAre0Based whether the query range points are 0-based indexes
     *                         (else false is 1-based indexing).
     * @return
     */
    public long sum(int l, int r, boolean indexesAre0Based) {
        if (indexesAre0Based) {
            ++r;
            ++l;
        }
        long s1 = sum(r);
        long s2 = sum(l - 1);
        return s1 - s2;
        //return sum(r) - sum(l - 1);
    }

    //r.t.c. O(sqrt(n)) where n = a.length;
    protected long sum(int r) {
        long res = 0;
        for (int i = 0; i < r / blockSize; i++) {
            res += blocks[i];
        }
        for (int i = (r / blockSize) * blockSize; i < r; i++) {
            res += a[i];
        }
        return res;
    }

    /** O(1) update to set a[i] to v */
    public void set(int i, long v, boolean indexesAre0Based) {
        if (!indexesAre0Based) {
            --i;
        }
        int bIdx = i / blockSize;
        blocks[bIdx] -= a[i];
        a[i] = v;
        blocks[bIdx] += a[i];
    }

    /**
     * update the internal data for the given range and addition
     * @param updateArray array holding left, right, add
     * @param indexesAre0Based
     */
    public void updateAdd(int[] updateArray, boolean indexesAre0Based) {
        int l = updateArray[0];
        int r = updateArray[1];
        if (!indexesAre0Based) {
            --l;
            --r;
        }
        for (int j = l; j <= r; ++j) {
            set(j, a[j] + updateArray[2], indexesAre0Based);
        }
    }

    /**
     * r.t.c. is (updates.length * (query range))
     * update internal data to add updates[i][2]
     * @param updates array of [left, right, add] arrays.
     * @param indexesAre0Based
     */
    public void updateAdd(int[][] updates, boolean indexesAre0Based) {
        for (int i = 0; i < updates.length; ++i) {
            updateAdd(updates[i], indexesAre0Based);
        }
    }

}
