package algorithms.combPerm;

import java.util.Arrays;

/**
 * a class to select k objects from n such that the sum of the selections
 * is minimized.
 *
 * The problem is adapted from the Competitive Programmer's Handbook by Antti Laaksonen, Chap 10.5.
 * And the extra conditional logic needed is adapted from a post by Zeks on:
 *     https://stackoverflow.com/questions/71248436/optimal-selection-for-minimum-total-sum
 */
public class SelectKMinSum {

    /**
     * select exactly 1 element from each row of kByNChoices such that the total sum of
     * the choices is minimized.
     *
     * The problem uses dynamic programming and a power set iteration pattern to solve
     * with a runtime complexity of O(k*n*2^(k)).
     *
     * The use of dynamic programming solution quickly becomes faster than a
     * brute force C(n,k) comparison
     * of all possible combinations of k selections from n elements for n>10.
     *
     * @param kByNChoices a data array with k rows and n columns.
     * @return the minimum sum of selecting exactly one item from each row.
     */
    public static long selectKMinSum(int[][] kByNChoices) {
        int k = kByNChoices.length;
        int n = kByNChoices[0].length;

        long[][] tab = new long[1<<k][n];
        long sentinel = Long.MAX_VALUE;
        for (long[] t : tab) {
            Arrays.fill(t, sentinel);
        }

        // init for 1st column of kByNChoices
        for (int iK = 0; iK < k; iK++) {
            tab[1<<iK][0] = kByNChoices[iK][0];
        }

        for (int i = 1; i < n; i++) {

            for (int s = 0; s < 1<<k; s++) {

                tab[s][i] = tab[s][i-1];

                for (int iK = 0; iK < k; iK++) {

                    if ((s&(1<<iK)) != 0) {
                        int sInclude = s^(1<<iK);

                        // too large (too many bits) for sequential from 0 to i:
                        if (Integer.bitCount(sInclude) > i) continue;

                        // add only if we haven't selected it already in s, and that iK is the
                        // only addition to s
                        if (Integer.bitCount(s ^ sInclude) != 1) continue;

                        // include current selection by including previous and adding to it.
                        // exclude by not changing (above we already set to previous without
                        // adding current)
                        long includePlusCurr = tab[sInclude][i-1] != sentinel ?
                                tab[sInclude][i-1] + kByNChoices[iK][i] : kByNChoices[iK][i];

                        if (tab[s][i] == sentinel || includePlusCurr <= tab[s][i]) {
                            tab[s][i] = includePlusCurr;
                        }
                    }
                }
            }
        }

        return tab[(1<<k)-1][n-1];
    }
}
