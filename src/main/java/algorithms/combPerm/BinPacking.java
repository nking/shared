package algorithms.combPerm;

/**
 * given an array and a constraint that each bin sum must be less than
 * or equal to capacity, find the minimum number of bins to place each
 * element in the array.
	 */
	public class BinPacking {

    /**
      given an array and a constraint that each bin sum must be less that
      or equal to capacity, find the minimum number of bins to place each
      element in the array.
     The algorithm uses dynamic programming and a bitstring power set pattern to solve w/ a
     r.t.c. of O(n*(2^n)) which is much better than trying all permutations for n!.
     The s.c. is O(n).
     <pre>
     The initial problem is adapted from Competitve Programmer's Handbook by Antti Laaksonen
      chapt 10.5.
     corrections are made here for the initialization. Also changed update of table
     to preserve integrity of previous calculations.
     </pre>
     * @param weights array of weights.   each weight must be .LEQ. capacity.
     * @param capacity weight capacity for each bin
     * @return the minimum number of bins to place each element in the
     * array.
     */
    public static int minNumberOfBins(int[] weights, int capacity) {
        int n = weights.length;

        // tab[S]:  S enumerates all subsets of size n as a bitstring, with set bits being
        //     the indexes of weights which are present in the calculation.
        // tab[S][0]: column 0 of each row holds the number of bins used for tab[s].
        // tab[S][1]: column 1 of each row holds the summed weight of the last bin not
        //            yet filled for this set S.
        int[][] tab = new int[1<<n][2];
        tab[0][0] = n+1;//1; // init with infeasible value
        tab[0][1] = 0;

        // init for 1 bin with each weight
        for (int i = 0; i < n; ++i) {
            tab[1<<i][0] = 1;
            tab[1<<i][1] = weights[i];
        }

        int[] tmpIncl = new int[2];

        // we solve for each subset successively.  start at s=1 so that prev is already init s=0
        for (int s = 1; s < (1<<n); ++s) {

            // if not a power of 2, hasn't been initialized yet:
            if ((s & (s-1)) != 0) {
                tab[s][0] = n + 1; // init with infeasible value
                tab[s][1] = 0;
            }

            for (int i = 0; i < n; ++i) {

                if ((s & (1<<i)) == 0) continue;// skip if S doesn't includes i

                int incl = s ^ (1 << i); // prev set is s minus i

                if (i == incl) continue; // cannot add item to its own set which already has item

                System.arraycopy(tab[incl], 0, tmpIncl, 0, 2);

                if (tmpIncl[1] + weights[i] <= capacity) {// add to existing latest bin
                    tmpIncl[1] += weights[i];
                } else { // new bin
                    ++tmpIncl[0];
                    tmpIncl[1] = weights[i];
                }
                // assign to S the best of current and candidate + current weight
                if ((tmpIncl[0] == tab[s][0] && tmpIncl[1] < tab[s][1]) || (tmpIncl[0] < tab[s][0])) {
                    tab[s][0] = tmpIncl[0];
                    tab[s][1] = tmpIncl[1];
                }
            }
        }
        return tab[(1<<n)-1][0];
    }
}
