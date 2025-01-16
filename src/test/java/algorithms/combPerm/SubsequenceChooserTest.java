package algorithms.combPerm;

import algorithms.misc.MiscMath0;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

public class SubsequenceChooserTest extends TestCase {

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

        assertEquals(np, out.length);

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
