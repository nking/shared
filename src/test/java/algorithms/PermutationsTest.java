package algorithms;

import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PermutationsTest extends TestCase {
    
    public PermutationsTest(String testName) {
        super(testName);
    }
    
    public void testPermute() {
        // public void permute(int[] set, int[][] outPermutations) {
    
        int[] set = new int[]{1, 2, 3};
        int[][] outPermutations = new int[6][];
        for (int i = 0; i < outPermutations.length; ++i) {
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
        int i, j;
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
