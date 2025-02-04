package algorithms.combPerm;

import algorithms.misc.MiscMath0;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

import java.util.*;

import java.util.Random;

public class SubsequenceChooserTest extends TestCase {

    public void testRand() {
        long seed = System.nanoTime();
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);
        int nTests = 10;
        for (int nTest = 0; nTest < nTests; ++nTest) {
            int n = 4 + rand.nextInt(4);
            int k = 2 + rand.nextInt(n-2);
            int[] a = new int[n];
            TIntSet aSet = new TIntHashSet(); // might be duplicates.  only using this to check output
            for (int i = 0; i < n; ++i) {
                a[i] = rand.nextInt();
                while (a[i] == 0) { // avoid 0 to check output for null entry
                    a[i] = rand.nextInt();
                }
                aSet.add(a[i]);
            }
            int np = (int)MiscMath0.computeNDivNMinusK(n, k);
            SubsequenceChooser sc = new SubsequenceChooser();
            int[][] out = sc.calcSubSequences(a, k);
            //System.out.printf("n=%d, k=%d, nIter=%d, npk=%d\n", a.length, k, sc.nIter, np);
            assertEquals(np, out.length);
            // assert all unique sequences
            Set<String> set = new HashSet<>();
            for (int[] o : out) {
                set.add(Arrays.toString(o));
                // assert all members are in a[i]
                for (int _o : o) {
                    assertTrue(aSet.contains(_o));
                }
            }
            assertEquals(np, set.size());
        }
    }

    public void test0() {
        int[] a = new int[]{1,2,3,4};
        int k = 2;
        int np = (int)MiscMath0.computeNDivNMinusK(a.length, k);

        TIntSet expected = new TIntHashSet();
        expected.add(12);
        expected.add(13);
        expected.add(14);
        expected.add(21);
        expected.add(23);
        expected.add(24);
        expected.add(31);
        expected.add(32);
        expected.add(34);
        expected.add(41);
        expected.add(42);
        expected.add(43);

        SubsequenceChooser sc = new SubsequenceChooser();
        int[][] out = sc.calcSubSequences(a, k);

        //System.out.printf("n=%d, k=%d, nIter=%d\n", a.length, k, sc.nIter);

        assertEquals(expected.size(), out.length);

        for (int i = 0; i < np; ++i) {
            int s = convertToInteger(out[i]);
            assertTrue(expected.remove(s));
        }
        assertTrue(expected.isEmpty());

    }

    int convertToInteger(int[] a) {
        int total = 0;
        int n = a.length;
        for (int i = n-1; i > -1; --i) {
            total += a[i] * Math.pow(10, n-i-1);
        }
        return total;
    }

}
