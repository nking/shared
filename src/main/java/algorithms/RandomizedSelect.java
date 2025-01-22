package algorithms;

import algorithms.util.FormatArray;

import java.util.Arrays;
import java.util.Random;

public class RandomizedSelect {

     //TODO: consider improving this with the Chapter notes for chapter 9 which
     // uses Floyd & Rivest 1975, "Expected Time Bounds for Selection: improvements
     // see MedianOfMediansSelect for the changes to 0-based indexing
    /**
     *  find the value with rank rank in array a with average runtime complexity O(n) where n = a.length
     *  and rank as a rank is 0-based.
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
     * @param rank the rank of the item to select in 0-based numbering.
     * @return the value of a's rank rank item where rank is 0-based
     */
    static double select(double[] a, int idxLo, int idxHi, int rank, Random rand) {
        if (idxLo < 0 || idxLo >= a.length) {
            throw new IllegalArgumentException("idxLo is out of bounds");
        }
        if (idxHi < 0 || idxHi >= a.length) {
            throw new IllegalArgumentException("idxHi is out of bounds");
        }
        if (rank < 0 || rank > idxHi) {
            throw new IllegalArgumentException("rank is out of bounds");
        }
        if (idxHi < idxLo) {
            throw new IllegalArgumentException("idxHi < idxLo");
        }

        if (idxLo == idxHi) {
            return a[idxLo];
        }

        int idxPivot = RandomizedQuickSort.partition(a, idxLo, idxHi, rand);
        // idxPivot is pivotIndex w.r.t 0
        // k is its rank w.r.t. range [idxLo, idxHi]
        int k = idxPivot - idxLo;

        if (rank == k) {
            return a[idxPivot];
        } else if (rank < k) {
            return select(a, idxLo, idxPivot-1, rank, rand);
        } else {
            return select(a, idxPivot+1, idxHi, rank - k - 1, rand);
        }
    }

}
