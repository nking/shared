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
     * with a runtime complexity of O(k*n*(2^k)).
     *
     * The use of dynamic programming quickly becomes faster than a brute force C(n,k) comparison
     * of all possible combinations of k selections from n elements for n GT 10.
     *
     * @param kByNChoices a data array with k rows and n columns. one item from each
        row must be picked.
     * @return the minimum sum of selecting exactly one item from each row.
     */
    public static long selectKMinSum(int[][] kByNChoices) {
        int k = kByNChoices.length;
        int n = kByNChoices[0].length;

        long[] tabPrev = new long[1<<k];
        long sentinel = Long.MAX_VALUE;
        Arrays.fill(tabPrev, sentinel);

        // init for 1st column of kByNChoices
        for (int iK = 0; iK < k; iK++) {
            tabPrev[1<<iK] = kByNChoices[iK][0];
        }

        for (int i = 1; i < n; i++) {
            long[] tabCurr = new long[1<<k];
            for (int s = 0; s < 1<<k; s++) {
                tabCurr[s] = tabPrev[s];
                for (int iK = 0; iK < k; iK++) {

                    // use s only if it includes iK
                    if ((s&(1<<iK)) == 0) continue;

                    // difference in set bits between s and iK.  this is a candidate previous
                    // set to add to that does not include iK
                    int sPrev = s^(1<<iK);

                    // too large (too many bits) for sequential approach from 0 to i:
                    if (Integer.bitCount(sPrev) > i) continue;

                    // include current selection by including previous and adding to it.
                    // exclude by not changing (above we already set to previous without
                    // adding current)
                    long prevPlusCurr = tabPrev[sPrev] != sentinel ?
                            tabPrev[sPrev] + kByNChoices[iK][i] : kByNChoices[iK][i];

                    if (tabCurr[s] == sentinel || prevPlusCurr <= tabCurr[s]) {
                        tabCurr[s] = prevPlusCurr;
                    }
                    
                } // end loop over iK
            } // end loop over s
            tabPrev = tabCurr;
        }

        return tabPrev[(1<<k)-1];
    }
}
