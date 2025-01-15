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

    public static long[] createPrefixArray(int[] a) {
        int n = a.length;
        // create prefix array to use for range queries
        long[] psa = new long[n+1];
        for (int i = 0; i < n; ++i) {
            psa[i+1] = a[i];
        }
        for (int i = 1; i < psa.length; ++i) {
            psa[i] += psa[i-1];
        }
        return psa;
    }

    public static long[] createPrefixArray(long[] a) {
        int n = a.length;
        // create prefix array to use for range queries
        long[] psa = new long[n+1];
        System.arraycopy(a, 0, psa, 1, a.length);
        for (int i = 1; i < psa.length; ++i) {
            psa[i] += psa[i-1];
        }
        return psa;
    }

    public long[] sum(int[][] queries, boolean queriesAre0Based) {
        return sum(array, queries, queriesAre0Based);
    }

    protected static long[] sum(long[] prefixArray, int[][] queries, boolean queriesAre0Based) {
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
            if (i0 < 0 || i1 >= prefixArray.length) {
                throw new IllegalArgumentException("i0 or i1 out of bounds");
            }
            if (i0-1 < 0) {
                i0 = 0;
            }
            out[i] = prefixArray[i1] - prefixArray[i0-1];
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
    public void updateAdd(int[][] updates, boolean rangesAre0Based) {
         updateAdd(array, updates, rangesAre0Based);
    }

    protected static void updateAdd(long[] prefixArray, int[][] updates, boolean rangesAre0Based) {
        // create a difference array of the queries, then apply it to the prefix sum
        int[] diffs = new int[prefixArray.length];
        int i0, i1;
        for (int[] update : updates) {
            i0 = update[0];
            i1 = update[1];
            if (rangesAre0Based) {
                ++i0;
                ++i1;
            }
            if (i0 < 0 || i1 >= prefixArray.length) {
                throw new IllegalArgumentException("i0 or i1 out of bounds");
            }
            diffs[i0] += update[2];
            if (i1+1 < prefixArray.length) {
                diffs[i1 + 1] -= update[2];
            }
        }

        // update prefix array
        int sum = 0;
        int curr = 0;
        for (int i = 0; i < prefixArray.length; ++i) {
            curr += diffs[i];
            sum += curr;
            prefixArray[i] += sum;
        }
    }

    /**
     * given an array 'a' of data, and array updateAndQuery of update ranges,
     * update the array 'a' by the update add amounts, and return sums of the ranges for each update.
     * each updateAndQuery row is a range start, end, inclusive, and amount to add to the region.
     * The runtime complexity is O(updates.length + a.length).
     * @param updateAndQuery an array of [start index, end index, add value].
     * @param rangesAre0Based if true, the range indexes in the updates are 0-based
     *                        indexes, else are 1-based indexes.
     *                        Note that if they are 1-based indexes, the updateAndQuery ranges
     *                        are converted to 0-based.
     * @raturn the sums of the ranges
     */
    public static long[] updateAddQuery(int[] a, int[][] updateAndQuery, boolean rangesAre0Based) {
        int n = a.length;
        if (!rangesAre0Based) {
            for (int[] update : updateAndQuery) {
                --update[0];
                --update[1];
            }
        }

        Misc.updateAddUsingDifferenceArray(a, updateAndQuery);

        // create prefix array to use for range queries
        long[] psa = createPrefixArray(a);

        return sum(psa, updateAndQuery, true);
    }
}
