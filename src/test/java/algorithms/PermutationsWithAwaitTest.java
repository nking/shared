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
        assertEquals(6, (int)np);
        
        int ii;
        
        int[][] results = new int[(int)np][];
        for (ii = 0; ii < results.length; ++ii) {
            results[ii] = new int[set.length];
        }
        
        PermutationsWithAwait p = new PermutationsWithAwait(set);
        
        for (i = 0; i < np; ++i) {
            p.getNext(results[(int)i]);
        }
        
        //System.out.printf("results=\n%s\n", FormatArray.toString(results, "%d"));
        
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
    }
}
