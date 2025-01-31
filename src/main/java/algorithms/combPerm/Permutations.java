package algorithms.combPerm;

import algorithms.misc.MiscMath0;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
/**
 *
 * @author nichole
 */
public class Permutations {

    // not counting array copies, just loop iterations
    protected static int nIter0 = 0;
    protected static int nIter1 = 0;

    /**
     permute the given set of numbers and store each result in a row of outPermutations.
     
        from https://en.wikipedia.org/wiki/Heap%27s_algorithm
        who reference:
        Sedgewick, Robert. "a talk on Permutation Generation Algorithms
        http://www.cs.princeton.edu/~rs/talks/perms.pdf
     @param set the set of numbers to generate all permuations for.
     * NOTE that set.length at most can be 12 due to the limit in the 
     * length of an array in java (which is outPermutations.length).
     * other methods can be created upon need.
     @param outPermutations output variable to fill with the permutations
    */
    public static void permute(int[] set, int[][] outPermutations) {
        
        long np = MiscMath0.factorial(set.length);
        int n = set.length;
        if (n > 12) {
            throw new IllegalArgumentException("set.length must be 12 or less "
            + " so that factorial(s.length) can fit in a java array");
        }
        if (outPermutations.length != np) {
            throw new IllegalArgumentException("outPermutations length must be set.length!");
        }
        if (outPermutations[0].length != n) {
            throw new IllegalArgumentException("outPermutations[0].length must be set.length");
        }
        
        set = Arrays.copyOf(set, n);
        
        int[] c = new int[n];
        
        int oc = 0;
        //output(A)
        outPermutations[oc] = Arrays.copyOf(set, n);
        oc++;

        /*
        procedure recursive(k : integer, A : array of any):
            if k = 0 then {
                output(A)
                return;
            }
            // Recursively call once for each k
            for i := 0; i < k; i += 1 do
                recursive(k - 1, A)
                // avoid swap when i==k-1
                if (i < k - 1)
                    // swap choice dependent on parity of k
                    if k is even then
                        swap(A[i], A[k-1])
                    else
                        swap(A[0], A[k-1])
                    end if
                end if
            end for
        */

        nIter0 = 0;

        int i = 0;
        int swap;
        while (i < n) {
            if (c[i] < i) {
                if ((i & 1) != 1) {
                    // i is even number
                    swap = set[0];
                    set[0] = set[i];
                    set[i] = swap;
                } else {
                    swap = set[c[i]];
                    set[c[i]] = set[i];
                    set[i] = swap;
                }
                outPermutations[oc] = Arrays.copyOf(set, n);
                oc++;

                ++nIter0;
                
                //Swap has occurred ending the for-loop. Simulate the increment 
                //of the for-loop counter
                c[i] += 1;
                //Simulate recursive call reaching the base case by bringing the 
                //pointer to the base case analog in the array
                i = 0;
            } else {
                //Calling generate(i+1, A) has ended as the for-loop terminated. 
                //Reset the state and simulate popping the stack by incrementing the pointer.
                c[i] = 0;
                i++;
            }
        }
    }

    static int findSuccessor(int[] a, int srch, int lo, int hi) {
        int sIdx = lo;
        for (int i = lo + 1; i <= hi; i++){
            if (a[i] > srch && a[i] < a[sIdx]) {
                sIdx = i;
            }
        }
        return sIdx;
    }

    /**
     * permute the array, lexicographically.
     <pre>
     adapted from
     Swarn Pallav Bhaskar, https://www.geeksforgeeks.org/lexicographic-permutations-of-string/
     </pre>
     * @param a array of numbers to permuta
     * @param ignoreDuplicates if true, will not print duplicate sequences.
     *                         e.g. output will contain [1,1,2] once instead of twice.
     * @return list of permutations, in sorted order
     */
    public static List<int[]> permuteLexicographically(int[] a, boolean ignoreDuplicates) {

        int n = a.length;

        Arrays.sort(a);

        List<int[]> out = new ArrayList<>();

        boolean isFinished = false;
        int[] prev = null;

        nIter1 = 0;

        while (!isFinished) {
            int[] cp = Arrays.copyOf(a, n);
            //when there is more than one of same number, we can avoid printing
            if (!ignoreDuplicates) {
                out.add(cp);
            } else if (prev == null || !Arrays.equals(prev, cp)) {
                out.add(cp);
            }
            prev = cp;

            ++nIter1;

            // find largest index i for which there is a larger value at a[i+1]
            int i = -1;
            for (i = n - 2; i >= 0; --i) {
                ++nIter1;
                if (a[i] < a[i + 1]) {
                    break;
                }
            }
            if (i == -1) {
                isFinished = true;
                continue;
            }

            // given a[i], find min value larger than it with idx >= i+1
            int sIdx = i + 1;
            for (int j = i + 2; j <= n-1; j++) {
                if (a[j] > a[i] && a[j] < a[sIdx]) {
                    sIdx = j;
                }
            }

            if (i != sIdx) { // swap
                a[i] ^= a[sIdx];
                a[sIdx] ^= a[i];
                a[i] ^= a[sIdx];
            }

            // reverse from i+1 to n-1
            for (int j = i+1, k = n-1; j <k; ++j, --k) {
                a[j] ^= a[k];
                a[k] ^= a[j];
                a[j] ^= a[k];
            }
        }
        return out;
    }

    public static boolean findNextLexicographically(int[] a) {
        int n = a.length;
        // find largest index i for which there is a larger value at a[i+1]
        int i = -1;
        for (i = n - 2; i >= 0; --i) {
            ++nIter1;
            if (a[i] < a[i + 1]) {
                break;
            }
        }
        if (i == -1) {
            return false;
        }

        // given a[i], find min value larger than it with idx >= i+1
        int sIdx = i + 1;
        for (int j = i + 2; j <= n-1; j++) {
            if (a[j] > a[i] && a[j] < a[sIdx]) {
                sIdx = j;
            }
        }

        if (i != sIdx) { // swap
            a[i] ^= a[sIdx];
            a[sIdx] ^= a[i];
            a[i] ^= a[sIdx];
        }

        // reverse from i+1 to n-1
        for (int j = i+1, k = n-1; j <k; ++j, --k) {
            a[j] ^= a[k];
            a[k] ^= a[j];
            a[j] ^= a[k];
        }
        return true;
    }

    public static boolean findPrevLexicographically(int[] a) {
        int n = a.length;
        int i = n-1;
        while (i > 0 && a[i] >= a[i-1]) {
            --i;
        }
        if (i == 0) {
            return false;
        }

        int sIdx = n-1;
        while (sIdx > -1 && a[sIdx] >= a[i-1]) {
           --sIdx;
        }

        if (i-1 >= 0 && (i-1) != sIdx) { // swap
            a[i-1] ^= a[sIdx];
            a[sIdx] ^= a[i-1];
            a[i-1] ^= a[sIdx];
        }

        // reverse a[i] thru n
        for (int j = i, k = n-1; j <k; ++j, --k) {
            a[j] ^= a[k];
            a[k] ^= a[j];
            a[j] ^= a[k];
        }

        return true;
    }

    public static List<int[]> recursivePermute(int[] a) {
        int n = a.length;
        List<int[]> out = new ArrayList<>();
        r(n, a, out);
        return out;
    }

    private static void r(int k, int[] a, List<int[]> out) {
        if (k == 0) {
            out.add(Arrays.copyOf(a, a.length));
            return;
        }
        int i = k - 1;
        for (int c = 0; c <= i; ++c) {
            r(i, a, out);
            if (c == i) break;
            if ((i&1) == 0) {
                // swap 0, i
                if (i != 0) {
                    a[0] ^= a[i];
                    a[i] ^= a[0];
                    a[0] ^= a[i];
                }
            } else {
                // swap c, i
                if (c != i) {
                    a[c] ^= a[i];
                    a[i] ^= a[c];
                    a[c] ^= a[i];
                }
            }
        }
    }
}
