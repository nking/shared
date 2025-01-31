package algorithms.combPerm;

import algorithms.misc.MiscMath0;
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

    public void testRecursiveRandom() {
        long seed = System.nanoTime();
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        int nTests = 3;
        for (int nTest = 0; nTest < nTests; ++nTest) {
            int n = 1 + rand.nextInt(7);
            int[] a = new int[n];
            for (int i = 0; i < n; ++i) {
                a[i] = rand.nextInt(1_000_000);
            }
            long nExp = MiscMath0.factorial(n);
            List<int[]> p = Permutations.recursivePermute(a);
            assertEquals(nExp, p.size());

            // redraw a uniquely
            Set<Integer> in = new HashSet<>();
            int i = 0;
            while (i < n) {
                int j = rand.nextInt(1_000_000);
                if (!in.contains(j)) {
                    a[i++] = j;
                    in.add(j);
                }
            }
            int[] b = Arrays.copyOf(a, n);
            Arrays.sort(b);
            p = Permutations.recursivePermute(a);
            assertEquals(nExp, p.size());
            Set<String> unique = new HashSet<>();
            for (int[] pi : p) {
                unique.add(Arrays.toString(pi));
                Arrays.sort(pi);
                assertTrue(Arrays.equals(b, pi));
            }
            assertEquals(nExp, unique.size());
        }
    }

    public void testNonRecursiveRandom() {
        long seed = System.nanoTime();
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        int nTests = 3;
        for (int nTest = 0; nTest < nTests; ++nTest) {
            int n = 1 + rand.nextInt(7);
            int[] a = new int[n];
            //draw a uniquely
            Set<Integer> in = new HashSet<>();
            int i = 0;
            while (i < n) {
                int j = rand.nextInt(1_000_000);
                if (!in.contains(j)) {
                    a[i++] = j;
                    in.add(j);
                }
            }
            int[] b = Arrays.copyOf(a, n);
            Arrays.sort(b);
            long nExp = MiscMath0.factorial(n);
            int[][] out = new int[(int)nExp][n];
            Permutations.permute(a, out);

            Set<String> unique = new HashSet<>();
            for (int[] pi : out) {
                unique.add(Arrays.toString(pi));
                Arrays.sort(pi);
                assertTrue(Arrays.equals(b, pi));
            }
            assertEquals(nExp, unique.size());
        }
    }
}
