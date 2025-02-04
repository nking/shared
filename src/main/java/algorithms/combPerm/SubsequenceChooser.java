package algorithms.combPerm;

import algorithms.misc.MiscMath0;
import gnu.trove.set.hash.TIntHashSet;

import java.util.Arrays;

/**
 * class to calculate subsequences of an array
 */
public class SubsequenceChooser {

    private int outIdx = 0;

    public SubsequenceChooser() {
    }

    /**
     * calculate all subsequences of size k from size a.
     * the r.t.c. is O(n!/((n-k)!)).
     * @param a array
     * @param k subsequence size
     * @return the subsequences of size k of a
     */
    public int[][] calcSubSequences(int[] a, int k) {
        int n = a.length;
        if (n < 1) {
            throw new IllegalArgumentException("n must be larger than 0");
        }
        if (k < 1) {
            throw new IllegalArgumentException("k must be larger than 0");
        }
        if (k > n) {
            throw new IllegalArgumentException("k must be less than or equal to n");
        }
        // n!/(n-k)! number of subsequences
        long npk = MiscMath0.computeNDivNMinusK(n, k);
        if (npk > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("the number of combinations is larger than max length of an array," +
                    "so this algorithm needs to be adjusted to return one element at a time");
        }
        int np = (int)npk;
        int[][] out = new int[np][k];
        this.outIdx = 0;

        recurseSeq(a, new int[k], 0, out, new TIntHashSet());

        return out;
    }

    int nIter = 0;
    private void recurseSeq(int[] a, int[] s, int i, int[][] out, TIntHashSet drawn) {
        ++nIter;
        if (i >= s.length) {
            out[outIdx] = Arrays.copyOf(s, s.length);
            ++outIdx;
            return;
        }

        for (int j = 0; j < a.length; ++j) {
            if (drawn.contains(j)) continue;
            s[i] = a[j];
            drawn.add(j);
            recurseSeq(a, s, i+1, out, drawn);
            drawn.remove(j);
        }
    }
}
