package algorithms;

import algorithms.misc.MiscMath0;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.stream.IntStream;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PermutationsTest extends TestCase {
    
    public PermutationsTest(String testName) {
        super(testName);
    }

    protected void swap(int[] a, int i, int j) {
        if (a[i] == a[j]) return;
        a[i] ^= a[j];
        a[j] ^= a[i];
        a[i] ^= a[j];
    }
    public void testPermute() {

        int i, j;

        /*
        Arrays.sort(out, (o1,o2)-> {
            for (int k  = 0; k < n; ++k) {
                int c = Integer.compare(o1[k], o2[k]);
                if (c != 0) return c;
            }
            // arrive here only for identical matrices
            return Integer.compare(o1[0], o2[0]);
        });
        */

        // public void permute(int[] set, int[][] outPermutations) {
    
        int[] set = new int[]{1, 2, 3};
        int[][] outPermutations = new int[6][];
        for (i = 0; i < outPermutations.length; ++i) {
            outPermutations[i] = new int[set.length];
        }
        
        int[][] expected = new int[6][];
        expected[0] = new int[]{1, 2, 3};
        expected[1] = new int[]{1, 3, 2};
        expected[2] = new int[]{2, 1, 3};
        expected[3] = new int[]{2, 3, 1};
        expected[4] = new int[]{3, 1, 2};
        expected[5] = new int[]{3, 2, 1};
        
        Permutations.permute(set, outPermutations);
        
        TIntSet found = new TIntHashSet();
        int[] ej, pi;
        for (i = 0; i < expected.length; ++i) {
            pi = outPermutations[i];
            for (j = 0; j < expected.length; ++j) {
                ej = expected[j];
                if (Arrays.equals(ej, pi)) {
                    found.add(j);
                }
            }
        }
        assertEquals(6, found.size());
    }
    
}
