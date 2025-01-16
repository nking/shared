package algorithms.combPerm;

import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.util.*;

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

    public void test3() {
        int[] m = new int[]{1,2,3};

        List<int[]> expected = new ArrayList<>();
        expected.add(new int[]{1,2,3});
        expected.add(new int[]{1,3,2});
        expected.add(new int[]{2,1,3});
        expected.add(new int[]{2,3,1});
        expected.add(new int[]{3,1,2});
        expected.add(new int[]{3,2,1});

        List<int[]> ans1 = Permutations.permuteLexicographically(Arrays.copyOf(m, m.length), true);
        assertEquals(expected, ans1);
        List<int[]> ans2 = Permutations.recursivePermute(Arrays.copyOf(m, m.length));
        assertEquals(expected, ans2);
    }

    public void test4() {
        int[] a = new int[] {1,2,4,1,3};
        int[] expAns = new int[]{1,2,4,3,1};

        int[] ans = Arrays.copyOf(a, a.length);
        Permutations.findNextLexicographically(ans);
        assertTrue(Arrays.equals(expAns, ans));

        Permutations.findPrevLexicographically(ans);
        assertTrue(Arrays.equals(a, ans));
    }

    protected void assertEquals(List<int[]> expected, List<int[]> ans) {
        assertEquals(expected.size(), ans.size());
        expected = new ArrayList<>(expected);
        for (int[] a : ans) {
            for (int i = 0; i < expected.size(); ++i) {
                int[] e = expected.get(i);
                if (Arrays.equals(a, e)) {
                    expected.remove(i);
                    break;
                }
            }
        }
        assertEquals(0, expected.size());
    }

    public void _testR() {
        int[] m = new int[]{1,2,3};
        r(m.length, m);
    }

    protected void r(int mIdx, int[] m) {
        if (mIdx == 0) {
            //System.out.printf("%s\n", Arrays.toString(m));
            return;
        }
        for (int i = 0; i < mIdx; ++i) {
            r(mIdx-1, m);
            int j = mIdx - 1;
            if (i < j) {
                if ((mIdx&1)!=1) { // even
                    m[i]      ^= m[j];
                    m[j] ^= m[i];
                    m[i]      ^= m[j];
                } else {//odd
                    if (0 != j) {
                        m[0]      ^= m[j];
                        m[j] ^= m[0];
                        m[0]      ^= m[j];
                    }
                }
            }
        }
    }
    
}
