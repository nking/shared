package algorithms.alignment;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;

/**
 * Knuth-Morris-Pratt string matcher.
 * runtime complexity is O(m) + O(n) where m is the length of the
 * pattern and n is the length of the text to search.
 *
 * algorithm is from Introduction to Algorithms, by Corman, Leiserson, Rivest, and Stein (CLRS)
 */
public class KMP {

    /**
     * find the indexes of t where the complete matches of p begin.
     *
     * the method implements KMP-Matcher from
     * Introduction to Algorithms, by Corman, Leiserson, Rivest, and Stein (CLRS)
     * @param p pattern of characters to find in t
     * @param t text to be searched for pattern p
     * @return the indexes of t where the complete matches of p begin.
     */
    public static int[] findPatternInText(char[] p, char[] t) {
        int m = p.length;
        int n = t.length;
        int[] pi = computePrefixFunction(p);
        TIntList idxs = new TIntArrayList();

        int q = 0;
        for (int i = 0; i < n; ++i) {
            while (q > 0 && p[q] != t[i]) {
                q = pi[q - 1];
            }
            if (p[q] == t[i]) {
                ++q;
            }
            if (q == m) {
                idxs.add(i - m + 1);
                q = pi[q - 1];
            }
        }
        return idxs.toArray();
    }

    /**
     * compute the prefix function which returns an offset character index of match of repeated pattern
     * @param p the pattern
     * @return prefix function
     */
    protected static int[] computePrefixFunction(char[] p) {
        int m = p.length;
        int[] pi = new int[m];
        int k = 0;
        for (int q = 1; q < m; ++q) {
            while (k > 0 && p[k] != p[q]) {
                k = pi[k - 1] - 1;
            }
            if (p[k] == p[q]) {
                ++k;
            }
            pi[q] = k;
        }
        return pi;
    }
}
