package algorithms.optimization;

import algorithms.search.MiscBisectingSearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Longest increasing subsequence (LIS) is general term for algorithms that find
 * the largest length subsequence of increasing values
 * within an array, where a subsequence maintains the order of the given list, but
 * can exclude members that do not fit the ordering.
 *
 * The implementations here use "patience sort" to improve runtime complexities
 * from O(n^2) to O(n*log(n)) where n is the length of the given array.

 <em> Patience sort</em>
 * Patience sort is named after the game of solitaire where one makes piles of
 * sequences of sorted numbers .
 *
 * As one iterates over the input list, the number  is compared to the
 * leftmost pile for which the pile's smallest card is LEQ (for ascending sort, or GEQ for descending sort)
 * that number,
 * and the card is added to that pile or forms a new pile to right of it if cannot add to a left pile.
 *
 * Note that there are versions for strictly increasing, or non-decreasing.

 <pre>
 for length of longest incr subseq, ascending sort example:
 a = [2,1,7,3,9]
 we would have piles:
 pile 0   pile 1    pile 2
   2      7         9
   1      3

 Every LIS is (2,7,9), (2,3,9), (1,7,9), (1,3,9).
 The length of the longest increasing subsequence is 3.
 There are 4 of LIS (2*2*1=4)

 Another example:

 a = [3,2,1,7,6,9,8]
 we would have piles:
 pile 0   pile 1    pile 2
 3        7         9
 2        6         8
 1

 Every LIS is (3,7,9), (3,7,8), (3,6,9),(3,6,8)
 (2,7,9), (2,7,8), (2,6,9),(2,6,8)
 (1,7,9), (1,7,8), (1,6,9),(1,6,8)

 The length of the longest increasing subsequence is 3.

 There are 12 of LIS (3*2*2=12)

 </pre>
 If we are interested in only the length of the longest increasing sequence, we
 only need to store the comparison number (i.e. the representation) for each pile.
 That means we can store each pile representative in a a single list and
 the length of that list is the length of the longest increasing aubsequence
 and we replace the representative with the new number that would be added to the pile.
 We find the leftmost pile that a number belongs to by  using a pileary search
 of type floor if the resulting series should be monontonic (increasing but equal values
 are in same piles) or type ceiling if
 the resulting series should be strictly monotonic (equal values are in different piles).

 If we are interested in only the number of LISes and the length of them,
 we can store the pile representation (see paragraph above this) in a list,
 and we can store the number of items in a pile in a list.

 Patience sort can be used on multi-column data too.
 In that case, we first sort by the first dimension,
 and if the sort is ascending of a[i][0], the tie breakers must be descending on a[i][1].

 <pre>
 a = [ [4,1], [2,3], [2,2], [5,3] ]

 the longest sequence wherein the first and second column values of a[i] are GEQ the
 values in a[i+1]

 for GT and ascending sort we would have presorted list:
 sort ascending of a[i][0], the tie breakers must be descending on a[i][1]:
     results in a = [ [2,3], [2,2], [4,1], [6,3] ]
 then we use a pileary search of type ceiling to find were a[i][1] can be placed
 in patience piles of decreasing values

 pile 0                        pile 1
 [2,3] where 3 is repr         [6,3]
 [2,2] where 2 is new repr
 [4,1] where 1 is new repr

 The longest sequences are
    [[2,3], [5,3]],
    [[2,2], [5,3]],
    [[4,1], [6,3]].
 There are 3 longest sequences.
 The length of the longest sequence is 2.

 ------------

 for an ascending sort in which each column is GEQ the same column for the previous element,
 we use patience sort similarly, but have an additional dimension.
 If we store all possible sequences, the problem becomes compuationally large, but still
 O(n^2 * log(n)) at worst.

 example a = [ [1,1], [2,2], [2,3], [4,1], [5,1], [5,2], [6,3] ]

 [1,1] through [2,3] placement follow patience rules.
 [4,1] does not fit after [2,3].
       ceiling search in pile 0 returns index=1.
       so we copy pile 0 up to index 1 into a new pile and add [4,1] to it
       (or keep reference to previous pile last index to include)
 [5,1] needs to be compared to pile 1 and pile 0
       for pile 1 we can add it to end.
       for pile 0, a ceiling search return index 1 and we see that in pile 1 already
                   so we do nothing wile pile 0 results.
 [5,2] needs to be compared to pile 1, and pile 0
       for pile 1, we can add to end and note that pile 1 ref to pile 0 is (0)
       for pile 0, ceil search returns index = 2, and that is greater than the pile 0, 0
           reference in pile 1, so we create a new pile 2, add the pile 0 items (or reference)
           and then add [5,2]] to pile 2
 [6,3] needs to be compared to pile 2, pile 1, and pile 0
       for pile 2, can add to end.
       for pile 1, can add to end.
       for pile 0, can add to end.

     pile 0       pile 1                       pile 2
     [1,1]        [1,1] ref (pile0, idx 0)     [1,1]
     [2,2]        [4,1]                        [2,2] ref (pile0, idx 1)
     [2,3]        [5,1]                        [5,2]
     [6,3]        [5,2]                        [6,3]
                  [6,3]

 The longest sequence which is increasing, but not strictly increasing for both columns of array 'a'
 is then in pile 1 for this example.
 The length of the LIS is 5.
 There is only 1 way to compose the LIS for this example.

 The worse case runtime complexity would be O(n^2 * log(n))
 for all items in their own piles, e.g. [1,4], [2, 2], [3,1].

 </pre>

 */
public class LongestIncreasingSubsequence {

    /**
     * given a list of integers, return the indexes of the longest increasing subsequence
     * within the list, where item a[i+1] is GT item a[i].  If there is more than one sequence
     * with same max length, return one of them.
     <pre>
     a = new int[]{2,2,1,7,6};

     2      7
     2
     1      6
     by value produces LIS:
     [2,7], [2,6], [2,7], [2,6], [1,7], [1,6]
     by index:
     [0,3], [0,4], [1,3], [1,4], [2,3], [2,4]
     </pre>
     * @param a
     * @return
     */
    public static List<Integer> findAnyStrictlyIncreasing(int[] a) {

        List<Integer> pileReps = new ArrayList<>();
        List<Integer> pileRepIdxs = new ArrayList<>();
        int v;
        for (int i = 0; i < a.length; ++i) {
            v = a[i];
            if (pileReps.isEmpty() || v > pileReps.get(pileReps.size() - 1)) {
                // create new pile for v which is GT last element
                pileReps.add(v);
                pileRepIdxs.add(i);
            } else {
                int idx = MiscBisectingSearch.ceiling(pileReps, v);
                assert(idx >-1 && idx < pileReps.size());
                pileReps.set(idx, v);
                pileRepIdxs.set(idx, i);
            }
        }

        return pileRepIdxs;
    }

    /**
     * given a list of integers, find the length of the longest increasing sequence, and return all
     * sequences of that length (specifically, the indexes of the items).
     * The list values should be GEQ previous item values, that is
     * item a[i+1] is GEQ item a[i].
     * @param a
     * @return
     */
    public static List<int[]> findAllStrictlyIncreasing(int[] a) {

        List<Integer> pileReps = new ArrayList<>();
        List<List<Integer>> pileIdxs = new ArrayList<>();
        int v;
        for (int i = 0; i < a.length; ++i) {
            v = a[i];
            if (pileReps.isEmpty() || v > pileReps.get(pileReps.size() - 1)) {
                // create new pile if v > last element
                pileReps.add(v);
                pileIdxs.add(new ArrayList<>());
                pileIdxs.get(pileIdxs.size() - 1).add(i);
            } else {
                int idx = MiscBisectingSearch.ceiling(pileReps, v);
                assert(idx >-1 && idx < pileReps.size());
                pileReps.set(idx, v);
                pileIdxs.get(idx).add(i);
            }
        }

        // return every combination of pileIdxs

        List<int[]> out = new ArrayList<>();
        int[] curr = new int[pileIdxs.size()];
        combineRecursively(0, curr, out, pileIdxs);

        int nSize = 1;
        for (int i = 0; i < pileIdxs.size(); ++i) {
            nSize *= pileIdxs.get(i).size();
        }
        assert(out.size() == nSize);

        return out;
    }

    private static void combineRecursively(int i, int[] curr, List<int[]> out, List<List<Integer>> pileIdxs) {
        if (i == pileIdxs.size()) {
            out.add(Arrays.copyOf(curr, curr.length));
            return;
        }
        List<Integer> currPile = pileIdxs.get(i);
        for (int ii = 0; ii < currPile.size(); ++ii) {
            int v = currPile.get(ii);
            curr[i] = v;
            combineRecursively(i+1, curr, out, pileIdxs);
        }
    }

    /**
     * given a list of integers, return the indexes of the longest increasing subsequence
     * within the list, where increasing is defined as a[i+1] GEQ a[i].
     * @param a
     * @return an array of length 2 holding as first item, the max size of longest increasing subequence
     * and as second item, the number of max size LISes.
     */
    public static int[] findSizeAndNumberStrictlyIncreasing(int[] a) {

         /*
         a = [3,2,1,7,6,9,8]
         we would have piles:
         pile 0   pile 1    pile 2
         3        7         9
         2        6         8
         1

         The length of the longest increasing subsequence is 3.

         Every LIS is (3,7,9), (3,7,8), (3,6,9), (3,6,8)
         (2,7,9), (2,7,8), (2,6,9), (2,6,8)
         (1,7,9), (1,7,8), (1,6,9), (1,6,8)

         There are 12 of LIS (3*2*2=12)
         */
        List<Integer> pileReps = new ArrayList<>();
        List<Integer> pileSizes = new ArrayList<>();
        int v;
        for (int i = 0; i < a.length; ++i) {
            v = a[i];
            if (pileReps.isEmpty() || v > pileReps.get(pileReps.size() - 1)) {
                pileReps.add(v);
                pileSizes.add(1);
            } else {
                int idx = MiscBisectingSearch.ceiling(pileReps, v);
                assert(idx >-1 && idx < pileReps.size());
                pileReps.set(idx, v);
                pileSizes.set(idx, pileSizes.get(idx) + 1);
            }
        }

        int nSize = 1;
        for (int i = 0; i < pileSizes.size(); ++i) {
            nSize *= pileSizes.get(i);
        }

        return new int[]{pileReps.size(), nSize};
    }

    /**
     * given a list of integers, return the size of the longest increasing sequence within the list,
     * where increasing is defined as a[i+1] GEQ a[i].
     * @param a
     * @return an array of length 2 holding as first item, the max size of longest increasing subequence
     * and as second item, the number of max size LISes.
     */
    public static int findSizeStrictlyIncreasing(int[] a) {

         /*
         a = [3,2,1,7,6,9,8]
         we would have piles:
         pile 0   pile 1    pile 2
         3        7         9
         2        6         8
         1

         The length of the longest increasing subsequence is 3.

         Every LIS is (3,7,9), (3,7,8), (3,6,9), (3,6,8)
         (2,7,9), (2,7,8), (2,6,9), (2,6,8)
         (1,7,9), (1,7,8), (1,6,9), (1,6,8)

         There are 12 of LIS (3*2*2=12)
         */
        List<Integer> pileReps = new ArrayList<>();
        int v;
        for (int i = 0; i < a.length; ++i) {
            v = a[i];
            if (pileReps.isEmpty() || v > pileReps.get(pileReps.size() - 1)) {
                pileReps.add(v);
            } else {
                int idx = MiscBisectingSearch.successor(pileReps, v);
                assert(idx >-1 && idx < pileReps.size());
                pileReps.set(idx, v);
            }
        }

        return pileReps.size();
    }

    /**
     * given a list of integers, return the indexes of the longest increasing subsequence
     * within the list, items are non-decreasing. 
     <pre>

     example: a = [7,6, 2,2,1,7,6]

     expected sequences:
      7,7
      6,7
      6,6
      2,2,7   <== LIS
      2,2,6   <== LIS
      2,7
      2,6
      1,7
      1,6

     needs multiple lists of patience piles and to remember top and bottom of pile

     list 0
     pile 0   line 1
     7         7
     6         cannot add 2nd 6 here as it is not >= 7, so it gets a new list
     2
     1

     start new list 1
     pile 0    pile 1   pile 2
     2          2        7
                         6

     start new list 2
     pile 0    pile 1
     6          6
     2
     1

     </pre>
     * @param a
     * @return
     */
    public static List<int[]> _findAll(int[] a) {

        throw new UnsupportedOperationException("not ready for use");

    }

}
