package algorithms.range;

import java.util.TreeSet;

public class Misc {

    /**
     * update array a with the updates.  r.t.c. is O(Q + N) where Q = updates.length and N = a.length;
     * see PrefixSumArray.java for query methods.
     * @param a array to update in place
     * @param updates array of index range inclusive and add amount as [indexLo, indexHi, addValue]
     *                indexes are expected to be 0-based.
     */
    public static void updateAddUsingDifferenceArray(int[] a, int[][] updates) {
        int n = a.length;

        int[] diffs = new int[n];
        for (int[] update : updates) {
            diffs[update[0]] += update[2];
            if (update[1] + 1 < n) {
                diffs[update[1] + 1] -= update[2];
            }
        }

        int sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += diffs[i];
            a[i] += sum;
        }

    }

}
