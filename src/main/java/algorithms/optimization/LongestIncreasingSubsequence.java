package algorithms.optimization;

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
 * intervals of sorted numbers .
 *
 * As one iterates over the input list, the number  is compared to the
 * leftmost pile for which the pile's smallest card is LEQ (for ascending sort, or GEQ for descending sort)
 * that number,
 * and the card is added to that pile or formd a new pile to right of it if cannot add to a left.

 <pre>
 for length of longest incr subseq, ascending sort example:
 a = [2,1,7,3,9]
 we would have piles:
 pile 0   pile 1    pile 2
   2      7         9
   1      3

 The length of the longest increasing subsequence is 3.
 Every LIS is (2,7,9), (2,3,9), (1,7,9), (1,3,9).
 There are 4 of LIS (2*2*1=4)

 a = [3,2,1,7,6,9,8]
 we would have piles:
 pile 0   pile 1    pile 2
 3        7         9
 2        6         8
 1

 The length of the longest increasing subsequence is 3.

 Every LIS is (3,7,9), (3,7,8), (3,6,9),(3,6,8)
 (2,7,9), (2,7,8), (2,6,9),(2,6,8)
 (1,7,9), (1,7,8), (1,6,9),(1,6,8)

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
 If we store all possible sequences, the problem becomes exponential.

 example a = [ [1,1], [2,2], [2,3], [4,1], [5,1], [5,2], [6,3] ]

 [1,1] through [2,3] placement follow patience rules.
 [4,1] does not fit after [2,3].
       successor search in pile 0 returns index=1.
       so we copy pile 0 up to index 1 into a new pile and add [4,1] to it
       (or keep reference to previous pile last index to include)
 [5,1] needs to be compared to pile 1 and pile 0
       for pile 1 we can add it to end.
       for pile 0, a succ search return index 1 and we see that in pile 1 already
                   so we do nothing wile pile 0 results.
 [5,2] needs to be compared to pile 1, and pile 0
       for pile 1, we can add to end and note that pile 1 ref to pile 0 is (0)
       for pile 0, succ search returns index = 2, and that is greater than the pile 0, 0
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

 The longest sequence which is increasing, but not strictly increasing for both columns of a
 is then in pile 1.
 The length of the LIS is 5.
 There is only 1 way to compose it.


 </pre>

 */
public class LongestIncreasingSubsequence {


    /**
     * given a list of integers, return the indexes of the longest increasing subsequence
     * within the list.
     * @param a
     * @return
     */
    public List<List<Integer>> findIndexes(int[] a) {
        throw new UnsupportedOperationException("not yet implemented");
    }

    /**
     * given a list of integers, return the indexes of the longest increasing subsequence
     * within the list.
     * @param a
     * @return an array of length 2 holding as first item, the max size of longest increasing subequence
     * and as second item, the number of max size LISes.
     */
    public int[] findSizeAndNumber(int[] a) {
        throw new UnsupportedOperationException("not yet implemented");
    }

    /**
     * given a list of integers, return the indexes of the longest increasing subsequence
     * within the list.
     * @param a
     * @return
     */
    public int findSize(int[] a) {
        throw new UnsupportedOperationException("not yet implemented");
    }
}
