package algorithms.bipartite;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.Map;

public class GaleShapley {

    /**
     * given two groups a and b and ranked preferences of each element for elements
     * in the opposite group, match each element of a with and element in b
     * in a "stable" way (at least one prefers their assigned partner).
     *
     * The results tend towards group 'a' receiving their most prefered pairings while
     * group 'b' receives their least-preferred pairings.
     *
     <pre>
     CLRS chap 25.2, Gale-Shapley algorithm
     </pre>

     Note that both aPrefs and bPrefs have to rank every member of the other list.

     * @param aPrefs each element is a ranked list of b indexes ordered by preference.
     *               note that the algorithm modifies the list.
     * @param bPrefs each element is a ranked list of b preferences where the
     *               Map has key=element in a, value = rank (preference).
     * @return
     */
    public static int[] match(LinkedList<Integer>[] aPrefs, Map<Integer, Integer>[] bPrefs) {
        int n = aPrefs.length;
        if (bPrefs.length != n) {
            throw new IllegalArgumentException("aPrefs.length must == bPrefs.length");
        }
        int[] mAB = new int[n];
        Arrays.fill(mAB, -1);

        int[] mBA = new int[n];
        Arrays.fill(mBA, -1);

        int nAMatched = 0;
        while (nAMatched != n) {

            for (int a0 = 0; a0 < n; ++a0) {
                if (mAB[a0] == -1) {
                    assert(aPrefs[a0] != null && !aPrefs[a0].isEmpty());
                    int b = aPrefs[a0].pollFirst();
                    if (mBA[b] == -1) {
                        mAB[a0] = b;
                        mBA[b] = a0;
                        ++nAMatched;
                    } else {
                        int a1 = mBA[b];
                        assert(mAB[a1] == b);
                        if (bPrefs[b].get(a1) > bPrefs[b].get(a0)) {
                            mAB[a1] = -1;
                            mAB[a0] = b;
                            mBA[b] = a0;
                        }
                        //  else a0 remains free
                    }
                }
            }
        }
        return mAB;
    }

}
