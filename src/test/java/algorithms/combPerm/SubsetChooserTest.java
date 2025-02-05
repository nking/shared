package algorithms.combPerm;

import algorithms.misc.MiscMath0;
import algorithms.util.PairInt;

import java.util.*;

import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SubsetChooserTest extends TestCase {
    
    public SubsetChooserTest() {
    }
    
    public void testGetNextSubset() throws Exception {
        
        int n = 3;
        int k = 2;
        
        SubsetChooser chooser = new SubsetChooser(n, k);
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(0, 1));
        expected.add(new PairInt(0, 2));
        expected.add(new PairInt(1, 2));
        
        int[] selectedIndexes = new int[k];
        
        while (chooser.getNextSubset(selectedIndexes) != -1) {
            int idx0 = selectedIndexes[0];
            int idx1 = selectedIndexes[1];
            if (idx0 > idx1) {
                int swap = idx0;
                idx0 = idx1;
                idx1 = swap;
            }
            PairInt p = new PairInt(idx0, idx1);
            assertTrue(expected.remove(p));
        }
        
        assertTrue(expected.isEmpty());
    }
    
    public void testGetNextSubset2() throws Exception {
        
        // k=2
        
        int n = 68;
        long np = (MiscMath0.computeNDivNMinusK(n, 2)) >> 1;
        
        Set<PairInt> expected = new HashSet<PairInt>();
        for (int i0 = 0; i0 < n; ++i0) {
            for (int i1 = (i0 + 1); i1 < n; ++i1) {
                PairInt p = new PairInt(i0, i1);
                expected.add(p);
            }
        }
                
        assertTrue(expected.size() == np);
        
        int[] selectedIndexes = new int[2];
        
        SubsetChooser chooser = new SubsetChooser(n, 2);
        
        while (chooser.getNextSubset(selectedIndexes) != -1) {
            int idx0 = selectedIndexes[0];
            int idx1 = selectedIndexes[1];
            if (idx0 > idx1) {
                int swap = idx0;
                idx0 = idx1;
                idx1 = swap;
            }
            PairInt p = new PairInt(idx0, idx1);
            assertTrue(expected.remove(p));
        }
        
        assertTrue(expected.isEmpty());
    }

    public void testRand() {
        long seed = System.nanoTime();
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);
        int nTests = 10;
        for (int nTest = 0; nTest < nTests; ++nTest) {
            int n = 4 + rand.nextInt(6);
            int k = 2 + rand.nextInt(n-2);
            int[] a = new int[n];
            TIntSet aSet = new TIntHashSet(); // might be duplicates.  only using this to check output
            for (int i = 0; i < n; ++i) {
                a[i] = rand.nextInt();
                while (a[i] == 0 || aSet.contains(a[i])) { // avoid 0 to check output for null entry
                    a[i] = rand.nextInt();
                }
                aSet.add(a[i]);
            }
            int nck = (int)MiscMath0.computeNDivKTimesNMinusK(n, k);
            List<int[]> out = SubsetChooser.calcSubSets(a, k);
            //System.out.printf("n=%d, k=%d, nIter=%d, nck=%d, nck*n=%d, nck*k=%d\n",
            //        a.length, k, SubsetChooser.nIter, nck, nck*n, nck*k);
            assertEquals(nck, out.size());
            // assert all unique sequences
            Set<String> set = new HashSet<>();
            for (int[] o : out) {
                set.add(Arrays.toString(o));
                // assert all members are in a[i]
                for (int _o : o) {
                    assertTrue(aSet.contains(_o));
                }
            }
            assertEquals(nck, set.size());
        }
    }
}
