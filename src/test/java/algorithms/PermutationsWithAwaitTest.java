package algorithms;

import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PermutationsWithAwaitTest extends TestCase {
    
    public PermutationsWithAwaitTest(String testName) {
        super(testName);
    }
    
    public void test0() throws InterruptedException {
        int[] set = new int[]{1,2,3};
        
        Set<TIntList> expected = new HashSet<TIntList>();
        expected.add(new TIntArrayList(new int[]{1, 2, 3}));
        expected.add(new TIntArrayList(new int[]{2, 1, 3}));
        expected.add(new TIntArrayList(new int[]{3, 1, 2}));
        expected.add(new TIntArrayList(new int[]{1, 3, 2}));
        expected.add(new TIntArrayList(new int[]{2, 3, 1}));
        expected.add(new TIntArrayList(new int[]{3, 2, 1}));
        
        long i;
        long np = MiscMath0.factorial(set.length);
        assertEquals(6L, np);
        
        int ii;
        
        int[][] results = new int[(int)np][];
        for (ii = 0; ii < results.length; ++ii) {
            results[ii] = new int[set.length];
            Arrays.fill(results[ii], -1);
        }
        
        PermutationsWithAwait p = new PermutationsWithAwait(set);
        assertTrue(p.hasNext());

        for (i = 0; i < np; ++i) {
            p.getNext(results[(int)i]);
        }
        boolean done = !p.hasNext();
        if (!done) {
            System.out.printf("results=\n%s\n", FormatArray.toString(results, "%d"));
            System.out.flush();
            fail("expecting getNext == false");
        }
        assertTrue(done);

        int[] r, r2;
        TIntList re2;
        for (ii = 0; ii < results.length; ++ii) {
            r = results[ii];
            re2 = null;
            for (TIntList re : expected) {
                r2 = re.toArray();
                if (Arrays.equals(r2, r)) {
                    re2 = re;
                    break;
                }
            }
            assertNotNull(re2);
            expected.remove(re2);
        }
        assertTrue(expected.isEmpty());

        // test that invocations past np do not deadlock
        p = new PermutationsWithAwait(set);
        assertTrue(p.hasNext());
        int n2 = (int)np + 10;
        results = new int[n2][];
        for (ii = 0; ii < results.length; ++ii) {
            results[ii] = new int[set.length];
            Arrays.fill(results[ii], -1);
        }
        for (i = 0; i < n2; ++i) {
            p.getNext(results[(int)i]);
        }
        done = !p.hasNext();
        if (!done) {
            System.out.printf("results=\n%s\n", FormatArray.toString(results, "%d"));
            System.out.flush();
            fail("expecting getNext == false");
        }
        assertTrue(done);
    }
}
