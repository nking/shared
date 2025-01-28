package algorithms.range;

import java.util.*;

/**
 * an algorithm used to process many queries of a static array (no modifications).
 *
 */
public class MosAlgorithm {

    /**
     * calculates the sum of all range queries of array 'a'.
     * r.t.c. is O((N+Q)*sqrt(N)) where N is a.length and Q is q.length.
     * Note that in contrast, a prefix array method is O(N + Q).
     *
     * Mo's algorithm is more helpful for problems when any use of a prefix array
     * would be multidimensional.  see other example methods in this class.
     <pre>
     adapted from
     author Ajay https://www.geeksforgeeks.org/mos-algorithm-query-square-root-decomposition-set-1-introduction/
     and adapted from
     https://cp-algorithms.com/data_structures/sqrt_decomposition.html#mos-algorithm
     </pre>
     * @param a array of integers.
     * @param q array of queries that are 0-based ranges of indexes of array 'a', inclusive.
     * @return sum of each query range.
     */
    public static long[] querySums(int[] a, int[][] q) {

        int blockSize = (int)Math.sqrt(a.length);

        // sort so that those within same block are grouped together, and those
        // outside are ascending or descending ordered by right range.
        // this sorting tries to reduce the number of times that currL and currR are changed
        // in the code below it.
        // r.t.c. O(Q*log(Q))
        Arrays.sort(q, new Comparator<int[]>() {
            public int compare(int[] x, int[] y) {
                if (x[0] / blockSize != y[0] / blockSize) {
                    return x[0] - y[0];
                }
                /*
                improvement from
                In odd blocks sort the right index in ascending order and
                in even blocks sort it in descending order.
                This will minimize the movement of right pointer, as the normal sorting will move the right pointer
                from the end back to the beginning at the start of every block. With the improved version this
                resetting is no more necessary.
                NLK: for the while loops below, reduction of nIter is when asc sort right for even,
                else descending sort right for odd.
                */
                if (((x[0]/blockSize) &1)==0) {
                    return x[1] - y[1];
                } else {
                    return y[1] - x[1];
                }
                //return x[1] - y[1];// default
            }
        });

        int currL = 0;
        int currR = 0;
        long currSum = 0;
        int qL, qR;
        long[] out = new long[q.length];

        for (int i = 0; i < q.length; i++) {
            // qL and qR values of current range
            qL = q[i][0];
            qR = q[i][1];

            while (currL < qL) {
                currSum -= a[currL];
                currL++;
            }
            while (currL > qL) {
                currSum += a[currL - 1];
                currL--;
            }
            while (currR <= qR) {
                currSum += a[currR];
                currR++;
            }
            while (currR > qR + 1) {
                currSum -= a[currR - 1];
                currR--;
            }
            out[i] = currSum;
        }
        return out;
    }

    /**
     * given integer array a, count the occurences of every number within query range for each query.
     * r.t.c. is O((N+Q)*sqrt(N)) where N is a.length and Q is q.length.
     <pre>
     adapted from
         author Ajay https://www.geeksforgeeks.org/mos-algorithm-query-square-root-decomposition-set-1-introduction/
     and adapted from
         https://cp-algorithms.com/data_structures/sqrt_decomposition.html#mos-algorithm
     </pre>
     * @param a array of integers
     * @param q array of queries of indexes into array 'a' as inclusive query ranges.
     * @return frequency maps for each query range are returned.
     */
    public static List<Map<Integer, Integer>> queryFrequencies(int[] a, int[][] q) {
        //if store as long[][][] we have ragged arrays of result[i] for q[i] = array of [number, count]
        // which is a ragged array, so we lose some of the benefits of efficient memory load/store
        //
        // List<Map<Integer, Integer>> takes more space but is more usable.

        List<Map<Integer, Integer>> out = new ArrayList<>();

        int n = a.length;

        int blockSize = (int)Math.ceil(Math.sqrt(n));

        Arrays.sort(q, new Comparator<int[]>(){
            public int compare(int[]x, int[] y) {
                if ((x[0] / blockSize) != (y[0] / blockSize)) {
                    return x[0] - y[0];
                }
                if (((x[0] / blockSize) & 1) == 0) {
                    return y[1] - x[1];
                } else {
                    return x[1] - y[1];
                }
            }
        });

        int currL = 0;
        int currR = -1;
        long currSum = 0;
        int qL, qR;

        Map<Integer, Integer> freqMap = new HashMap<>();

        for (int i = 0; i < q.length; i++) {
            // qL and qR values of current range
            qL = q[i][0];
            qR = q[i][1];
            while (currL > qL) {
                currL--;
                freqMap.put(a[currL], freqMap.getOrDefault(a[currL], 0) + 1);
            }
            while (currR < qR) {
                currR++;
                freqMap.put(a[currR], freqMap.getOrDefault(a[currR], 0) + 1);
            }
            while (currR > qR) {
                freqMap.put(a[currR], freqMap.getOrDefault(a[currR], 0) - 1);
                currR--;
            }
            while (currL < qL) {
                freqMap.put(a[currL], freqMap.getOrDefault(a[currL], 0) - 1);
                currL++;
            }
            //O(nUnique from a[minQ] to a[maxQ])
            out.add(writeMap(freqMap));
        }
        return out;
    }

    private static Map<Integer, Integer> writeMap(Map<Integer, Integer> freqMap) {
        Map<Integer, Integer> out = new HashMap<>();
        for (Map.Entry<Integer, Integer> entry : freqMap.entrySet()){
            if (entry.getValue() != 0) {
                out.put(entry.getKey(), entry.getValue());
            }
        }
        return out;
    }
}
