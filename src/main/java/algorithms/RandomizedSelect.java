package algorithms;

import algorithms.util.FormatArray;

import java.util.Arrays;
import java.util.Random;

public class RandomizedSelect {

    /**
     *  find the value with rank i in array a with average runtime complexity O(n) where n = a.length
     *  and i as a rank is 0-based.
     * The worst case runtime complexity is O(n^2).
     *
     <pre>
     references :
     Randomized Select, Section 9.2
     "Introduction to Algorithms" by
     Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and Clifford Stein
     </pre>
     * @param a an unsorted array
     * @param idxLo smallest index to of range
     * @param idxHi largest index of range, inclusive
     * @param i the rank of the item to select in 0-based numbering.
     * @return the value of a's rank i item where i is 0-based
     */
    static double select(double[] a, int idxLo, int idxHi, int i, Random rand) {
        if (idxLo < 0 || idxLo >= a.length) {
            throw new IllegalArgumentException("idxLo is out of bounds");
        }
        if (idxHi < 0 || idxHi >= a.length) {
            throw new IllegalArgumentException("idxHi is out of bounds");
        }
        if (i < 0 || i > idxHi) {
            throw new IllegalArgumentException("i is out of bounds");
        }
        if (idxHi < idxLo) {
            throw new IllegalArgumentException("idxHi < idxLo");
        }

        if (idxLo == idxHi) {
            return a[idxLo];
        }

        int q = RandomizedQuickSort.partition(a, idxLo, idxHi, rand);
        // q is pivotIndex w.r.t 0
        // k is its rank
        int k = q - idxLo;

        if (i == k) {
            return a[q];
        } else if (i < k) {
            return select(a, idxLo, q-1, i, rand);
        } else {
            return select(a, q+1, idxHi, i - k - 1, rand);
        }
    }

}
