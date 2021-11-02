package algorithms;

import algorithms.misc.MiscMath0;
import java.util.Arrays;

/**
 *
 */
public class Permutations {
    
    /**
     permute the given set of numbers and store each result in a row of outPermutations.
     
        from https://en.wikipedia.org/wiki/Heap%27s_algorithm
        who reference:
        Sedgewick, Robert. "a talk on Permutation Generation Algorithms
        http://www.cs.princeton.edu/~rs/talks/perms.pdf
     * @param set the set of numbers to generate all permuations for.
     * NOTE that set.length at most can be 12 due to the limit in the 
     * length of an array in java (which is outPermutations.length).
     * other methods can be created upon need.
     * @param outPermutations        
    */
    public static void permute(int[] set, int[][] outPermutations) {
        
        long np = MiscMath0.factorial(set.length);
        int n = set.length;
        if (n > 12) {
            throw new IllegalArgumentException("set.length must be 12 or less "
            + " so that factorial(s.ength) can fit in a java array");
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
}
