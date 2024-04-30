package algorithms.optimization;

import algorithms.search.MiscBisectingSearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.function.ToIntFunction;
import java.util.stream.IntStream;

/**
 * Longest increasing subsequence (LIS) is general term for algorithms that find
 * the largest length subsequence of increasing values
 * within an array, where a subsequence maintains the order of the given list, but
 * can exclude members that do not fit the ordering.
 *
 * The implementations here use "patience sort" on comparisons to improve runtime complexities
 * from O(n^2) to O(n*log(n)) where n is the length of the given array.

 <em> Patience sort</em>
 * Patience sort is named after the game of solitaire where one makes piles of
 * sequences of sorted numbers .
 *
 * As one iterates over the input list, the number  is compared to the
 * leftmost pile for which the pile's smallest card is LEQ
 * that number,
 * and the card is added to that pile or forms a new pile to right of it if cannot add to a left pile.
 *
 * Note that there are methods here for strictly increasing, or non-decreasing subsequences.

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

 <pre>
 --------------------------------
 </pre>
 Patience sort can be used on multi-column data too.
 In that case, we first sort by the first dimension,
 and if the sort is ascending of a[i][0], the tie breakers must be descending on a[i][1].

 <pre>
 a = [ [4,1], [2,3], [2,2], [5,3] ]

 the longest sequence wherein the first and second column values of a[i] are GEQ the
 values in a[i+1]

 for strictly increasing sequences:

 (1) sort ascending of a[i][0], the tie breakers must be descending on a[i][1]:
     results in a = [ [2,3], [2,2], [4,1], [6,3] ]
 (2) then binary search of type ceiling to find where a[i][1] can be placed
 in patience piles of decreasing values

 pile 0                        pile 1
 [2,3] where 3 is repr         [6,3]
 [2,2] where 2 is new repr
 [4,1] where 1 is new repr

 The longest sequences are
    [[2,3], [6,3]],
    [[2,2], [6,3]],
    [[4,1], [6,3]].
 There are 3 longest sequences.
 The length of the longest sequence is 2.

 <pre>
 --------------------------------
 </pre>

 For multi-column data and non-decreasing subsequences:

 TODO: rewrite this one

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
     * where increasing is defind as stricly increasing.  If there is more than one sequence
     * with same max length, returns one of them.
     <pre>
     <pre>
     a = new int[]{2,2,1,7,6};

     pile 0    pile 1
     2         7
     2
     1          6

     The result is 6 LIS of length 2.

     by value the sequences are:
     [2,7], [2,6], [2,7], [2,6], [1,7], [1,6]

     by index the sequences are:
     [0,3], [0,4], [1,3], [1,4], [2,3], [2,4]

     returns one of the sequences of LIS indexes.
     </pre>
     The runtime complexity is O(n * log(n)) where n = a.length.
     The space complexity is O(n).
     * @param a an array of integers
     * @return one of the possible many index sequences of longest increasing subsequence.
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
     * Increasing is defined as strictly increasing.
     <pre>
     a = new int[]{2,2,1,7,6};

     pile 0    pile 1
     2         7
     2
     1          6

     The result is 6 LIS of length 2.

     by value the sequences are:
     [2,7], [2,6], [2,7], [2,6], [1,7], [1,6]

     by index the sequences are:
     [0,3], [0,4], [1,3], [1,4], [2,3], [2,4]

     returns indexes of the LIS
     </pre>
     The runtime complexity is O(k * n * log(n)) where n = a.length and k
     is the product of the number of sequences in each pile.
     The space complexity is O(k*n).
     * @param a an array of integers
     * @return list of indexes of each LIS of maximum length
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
     * within the list, where increasing is strictly increasing.
     <pre>
     a = new int[]{2,2,1,7,6};

     pile 0    pile 1
     2         7
     2
     1          6

     The result is 6 LIS of length 2.

     by value the sequences are:
     [2,7], [2,6], [2,7], [2,6], [1,7], [1,6]

     by index the sequences are:
     [0,3], [0,4], [1,3], [1,4], [2,3], [2,4]

     returns new int[]{2, 6}

     </pre>
     The runtime complexity is O(n * log(n)) where n = a.length.
     The space complexity is O(n).
     * @param a
     * @return an array of length 2 holding as first item, the max size of longest increasing subequence
     * and as second item, the number of max size LISes.
     */
    public static int[] findSizeAndNumberStrictlyIncreasing(int[] a) {

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

    public class AscDescComparator implements Comparator<int[]> {
        @Override
        public int compare(int[] a1, int[] a2) {
            if (a1[0] == a2[0]) return 0;
            return Integer.compare(a2[1], a1[1]);
        }
    }

    /**
     find the longest strictly increasing subsequence of a where columns 0
     and columns 1 of row i are greater than columns 0 and 1, respectively of
     row (i-1).
     * The runtime complexity is O(k * n * log(n)) where n = a.length and k
     *      is the product of the number of sequences in each pile.
     *      The space complexity is O(k*n).
     * If only 1 LIS is needed not the entire enumeration,
     * or only the size of LIS is needed, can use methods
     * findAnyStrictlyIncreasing or findSizeAndNumberStrictlyIncreasing
     * for smaller runtime complexities of O(n * log(n)).
     @param a 2 dimensional array of length n X 2.  Though, rows of a can be longer than 2,
      *          only the first 2 columns of a row are used in comparisons.
     @return list of indexes of each longest LIS of maximum length
     */
    public static List<int[]> findAllStrictlyIncreasing(int[][] a) {

        int n = a.length;

        if (a[0].length < 2) {
            throw new IllegalStateException("each row length must be at least 2");
        }

        // created indexes sorted on ascending a[i][0] with tie breaking descending
        // sort on a[i][1]

        /*
        int[] sortedIdxs = IntStream.range(0, n)
                .boxed() // produces Stream<Integer> needed for Comparator
                .sorted(
                        //type argument of Comparator cannot be primitive
                       new Comparator<Integer>() {
                           @Override
                           public int compare(Integer i, Integer j) {
                                // a is external to this anonymous inner class
                                int c = Integer.compare(a[i][0], a[j][0]);
                                if (c != 0) {
                                    return c;
                                }
                                return Integer.compare(a[j][1], a[i][1]);
                           }
                       }
                )
                .mapToInt(new ToIntFunction<Integer>() {
                              @Override
                              public int applyAsInt(Integer value) {
                                  return value;
                              }
                          })
                .toArray();
        */
        // nice and compact, but produces Integer Objects in between, so this can be improved by writing a custom sort method
        int[] sortedIdxs = IntStream.range(0, n)
                .boxed() // produces Stream<Integer> needed for Comparator
                .sorted((i, j) -> {
                    int c = Integer.compare(a[i][0], a[j][0]);
                    if (c != 0) {
                        return c;
                    }
                    return Integer.compare(a[j][1], a[i][1]);
                })
                .mapToInt(element -> element).toArray();

        List<Integer> pileReps = new ArrayList<>();
        List<List<Integer>> pileIdxs = new ArrayList<>();
        int v2;
        int idx;
        int[] aI;
        for (int i0 = 0; i0 < n; ++i0) {
            idx = sortedIdxs[i0];
            aI = a[idx];

            // sorted already by a[i][0], now find where a[i][1] fits in patience piles
            v2 = aI[1];

            if (pileReps.isEmpty() || v2 > pileReps.get(pileReps.size() - 1)) {
                // create new pile if v > last element
                pileReps.add(v2);
                pileIdxs.add(new ArrayList<>());
                pileIdxs.get(pileIdxs.size() - 1).add(idx);
            } else {
                int _idx = MiscBisectingSearch.ceiling(pileReps, v2);
                assert(_idx >-1 && _idx < pileReps.size());
                pileReps.set(_idx, v2);
                pileIdxs.get(_idx).add(idx);
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

    /**
     find the longest strictly increasing subsequence of a where
     columns 0 and 1 of row a[i] are greater than columns 0 and 1, respectively,
     of row[i-1].
     If more than one sequence of indexes exists for the maximum length,
     only the last is returned.
     * The runtime complexity is O(n * log(n)) where n = a.length.
     * The space complexity is O(n).
     @param a 2 dimensional array of length n X 2.  Though, rows of a can be longer than 2,
      *          only the first 2 columns of a row are used in comparisons.
     @return indexes of one of the longest increasing subsequences
     */
    public static List<Integer> findAnyStrictlyIncreasing(int[][] a) {

        int n = a.length;

        if (a[0].length < 2) {
            throw new IllegalStateException("each row length must be at least 2");
        }

        // created indexes sorted on ascending a[i][0] with tie breaking descending
        // sort on a[i][1]

        // nice and compact, but produces Integer Objects in between, so this can be improved by writing a custom sort method
        int[] sortedIdxs = IntStream.range(0, n)
                .boxed() // produces Stream<Integer> needed for Comparator
                .sorted((i, j) -> {
                    int c = Integer.compare(a[i][0], a[j][0]);
                    if (c != 0) {
                        return c;
                    }
                    return Integer.compare(a[j][1], a[i][1]);
                })
                .mapToInt(element -> element).toArray();

        List<Integer> pileReps = new ArrayList<>();
        List<Integer> pileRepIdxs = new ArrayList<>();
        int v2;
        int idx;
        int[] aI;
        for (int i0 = 0; i0 < n; ++i0) {
            idx = sortedIdxs[i0];
            aI = a[idx];

            // sorted already by a[i][0], now find where a[i][1] fits in patience piles
            v2 = aI[1];

            if (pileReps.isEmpty() || v2 > pileReps.get(pileReps.size() - 1)) {
                // create new pile if v > last element
                pileReps.add(v2);
                pileRepIdxs.add(idx);
            } else {
                int _idx = MiscBisectingSearch.ceiling(pileReps, v2);
                assert(_idx >-1 && _idx < pileReps.size());
                pileReps.set(_idx, v2);
                pileRepIdxs.set(_idx, idx);
            }
        }

        return pileRepIdxs;
    }

    /**
     find the length and number of longest strictly increasing subsequence of a where
     columns 0 and 1 of row a[i] are greater than columns 0 and 1, respectively,
     of row[i-1].
     The runtime complexity is O(n * log(n)) where n = a.length.
     The space complexity is O(n).
     * @param a 2 dimensional array of length n X 2.  Though, rows of a can be longer than 2,
     *          only the first 2 columns of a row are used in comparisons.
     * @return the size of the LIS and the number of them as an integer array.
     */
    public static int[] findSizeAndNumberStrictlyIncreasing(int[][] a) {

        int n = a.length;

        if (a[0].length < 2) {
            throw new IllegalStateException("each row length must be at least 2");
        }

        // created indexes sorted on ascending a[i][0] with tie breaking descending
        // sort on a[i][1]

        // nice and compact, but produces Integer Objects in between, so this can be improved by writing a custom sort method
        int[] sortedIdxs = IntStream.range(0, n)
                .boxed() // produces Stream<Integer> needed for Comparator
                .sorted((i, j) -> {
                    int c = Integer.compare(a[i][0], a[j][0]);
                    if (c != 0) {
                        return c;
                    }
                    return Integer.compare(a[j][1], a[i][1]);
                })
                .mapToInt(element -> element).toArray();

        List<Integer> pileReps = new ArrayList<>();
        List<List<Integer>> pileIdxs = new ArrayList<>();
        int v2;
        int idx;
        int[] aI;
        for (int i0 = 0; i0 < n; ++i0) {
            idx = sortedIdxs[i0];
            aI = a[idx];

            // sorted already by a[i][0], now find where a[i][1] fits in patience piles
            v2 = aI[1];

            if (pileReps.isEmpty() || v2 > pileReps.get(pileReps.size() - 1)) {
                // create new pile if v > last element
                pileReps.add(v2);
                pileIdxs.add(new ArrayList<>());
                pileIdxs.get(pileIdxs.size() - 1).add(idx);
            } else {
                int _idx = MiscBisectingSearch.ceiling(pileReps, v2);
                assert(_idx >-1 && _idx < pileReps.size());
                pileReps.set(_idx, v2);
                pileIdxs.get(_idx).add(idx);
            }
        }

        int nSize = 1;
        for (int i = 0; i < pileIdxs.size(); ++i) {
            nSize *= pileIdxs.get(i).size();
        }

        return new int[]{pileReps.size(), nSize};
    }

    /**
     given a list of integers, return the indexes of the longest increasing subsequence
     within the list, items are non-decreasing.

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
     6
     2         cannot add 2nd 2 as it is not >= top of pile 0 = 7, so it gets a new list
     1         cannot add 2nd 6 as it is not >= top of pile 0 = 7, so it gets a new list

     list 1 for 2nd 2 that can't be added to list 0
     pile 0    pile 1   pile 2
     2          2
                cannot add 1 because it is not >= top of pile 0 = 2 so it gets a new list
                          7
                          6

     list 2 for 1 that cant be added to list 1
     pile 0    pile 1
     1          7
                6

     list 3 for 6 that cant be added to list 0
     pile 0    pile 1
     6          6
     2
     1

     </pre>
     * @param a
     * @return
     */
    public static List<int[]> _findAllNonDecreasing(int[] a) {

        throw new UnsupportedOperationException("not ready for use");

    }

}
