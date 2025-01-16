package algorithms.combPerm;

import algorithms.misc.MiscMath0;
import gnu.trove.set.hash.TIntHashSet;

import java.util.Arrays;

public class SubsequenceChooser {

    private int outIdx = 0;

    public SubsequenceChooser() {
    }

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
        // n!/(n-k)!
        long nnk = MiscMath0.computeNDivNMinusK(n, k);
        if (nnk > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("the number of combinations is larger than max length of an array," +
                    "so this algorithm needs to be adjusted to return one element at a time");
        }
        int np = (int)nnk;
        int[][] out = new int[np][k];
        this.outIdx = 0;

        recurseSeq(a, new int[k], 0, out, new TIntHashSet());

        return out;
    }

    private void recurseSeq(int[] a, int[] s, int i, int[][] out, TIntHashSet drawn) {
        if (i >= s.length) {
            out[outIdx] = Arrays.copyOf(s, s.length);
            ++outIdx;
            return;
        }

        for (int j = 0; j < a.length; ++j) {
            if (drawn.contains(a[j])) continue;
            s[i] = a[j];
            drawn.add(a[j]);
            recurseSeq(a, s, i+1, out, drawn);
            drawn.remove(a[j]);
        }
    }
}
