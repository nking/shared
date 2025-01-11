package algorithms.range;

import algorithms.util.FormatArray;

import java.util.Arrays;

public class PrefixSumArray {
    // length is a.length + 1
    // array[0] = 0.
    // a query for [i0,i1] is array[i1+1] - array[i0]
    // e.g. [1,2] = array[2+1] - array[1]
    protected long[] array;

    public PrefixSumArray(long[] a) {
        array = new long[a.length + 1];
        System.arraycopy(a, 0, array, 1, a.length);
        for (int i = 1; i < array.length; ++i) {
            array[i] += array[i-1];
        }
    }

    public long[] sum(int[][] queries, boolean queriesAre0Based) {
        long[] out = new long[queries.length];
        int i0, i1;
        for (int i = 0; i < queries.length; ++i) {
            int[] q = queries[i];
            i0 = q[0];
            i1 = q[1];
            if (queriesAre0Based) {
                ++i0;
                ++i1;
            }
            if (i0 < 0 || i1 >= array.length) {
                throw new IllegalArgumentException("i0 or i1 out of bounds");
            }
            if (i0-1 < 0) {
                i0 = 0;
            }
            out[i] = array[i1] - array[i0-1];
        }
        return out;
    }

    public long sum(int i0, int i1, boolean queriesAre0Based) {
        if (queriesAre0Based) {
            ++i0;
            ++i1;
        }
        if (i0 < 0 || i1 >= array.length) {
            throw new IllegalArgumentException("i0 or i1 out of bounds");
        }
        if (i0-1 < 0) {
            i0 = 0;
        }
        return array[i1] - array[i0-1];
    }

    /**
     * given an array of updates, update the internal data for it.
     * each update is a range start, end, and amount to add to the region.
     * The runtime complexity is O(updates.length + a.length).
     * @param updates an array of [start index, end index, add value].
     * @param rangesAre0Based if true, the range indexes in the updates are 0-based
     *                        indexes, else are 1-based indexes.
     */
    public void update(int[][] updates, boolean rangesAre0Based) {
        // create a difference of the queries, then apply it to the prefix sum
        int[] diffs = new int[array.length];
        int i0, i1;
        for (int[] update : updates) {
            i0 = update[0];
            i1 = update[1];
            if (rangesAre0Based) {
                ++i0;
                ++i1;
            }
            if (i0 < 0 || i1 >= array.length) {
                throw new IllegalArgumentException("i0 or i1 out of bounds");
            }
            diffs[i0] += update[2];
            if (i1+1 < array.length) {
                diffs[i1 + 1] -= update[2];
            }
        }

        // update prefix array
        int sum = 0;
        int curr = 0;
        for (int i = 0; i < array.length; ++i) {
            curr += diffs[i];
            sum += curr;
            array[i] += sum;
        }
    }
}
