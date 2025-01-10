package algorithms.range;

import java.util.Arrays;

/**
 * a class to handle minimum range queries.
 * with r.t.c. O(1) after a pre-processing r.t.c. of O(n*log(n))
 */
public class MinSparseTable {

    // note that since the table is built for minimum and not aggregate operations,
    // one could store indicies in min and retain the original array a
    // instead if useful.

    int[][] min;

    public MinSparseTable(int[] a) {
        int n = a.length;
        int n2 = 1 + (int)Math.floor(Math.log(n)/Math.log(2));

        min = new int[n][n2];

        /*
        we calc the minimum for index ranges of powers of 2
        col min[0] holds min of index ranges (0,0}, (1,1)...
        col min[1] holds min of index ranges (0,1), (1,2)...
        col min[2] ...                       (0,3), (1,4)

        https://www.topcoder.com/thrive/articles/Range%20Minimum%20Query%20and%20Lowest%20Common%20Ancestor

        TODO: for better locality, rewrite for transposed min.
        */
        int i, j;
        for (i = 0; i < n; ++i) {
            min[i][0] = a[i];
        }
        for (j = 1; j < n2; j++) {
            for (i = 0; (i + (1 << j)) <= n; i++) {
                min[i][j] = Math.min(min[i][j - 1], min[i + (1 << (j - 1))][j - 1]);
            }
        }
    }

    /**
     * given an array of query ranges, return the min within them.
     * @param queries array of ranges [left index, right index] inclusive
     * @return
     */
    public int[] min(int[][] queries) {
        int[] out = new int[queries.length];
        for (int i = 0; i < queries.length; ++i) {
            out[i] = min(queries[i][0], queries[i][1]);
        }
        return out;
    }
    public int min(int leftIdx, int rightIdx) {
        if (rightIdx < leftIdx) {
            throw new IllegalArgumentException("leftIdx must be <= rightIdx");
        }

        // largest power of 2 <= r
        int n2 = (int)(Math.log(rightIdx - leftIdx + 1)/Math.log(2));

        // min of [range][n2] where range is leftIdx to rightIdx - n2 + 1
        int i2 = rightIdx - (1<<n2) + 1;
        return Math.min(min[leftIdx][n2], min[i2][n2]);
    }
}
