package algorithms.graphs;

import algorithms.misc.MiscMath0;
import algorithms.util.PairInt;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CuthillMcKeeTest extends TestCase {
    
    public CuthillMcKeeTest(String testName) {
        super(testName);
    }
    
    public void testRCM() {
        Set<PairInt> gE = new HashSet<PairInt>();
        gE.add(new PairInt(0, 1));
        gE.add(new PairInt(1, 0));
        gE.add(new PairInt(0, 2));
        gE.add(new PairInt(2, 0));
        gE.add(new PairInt(0, 3));
        gE.add(new PairInt(3, 0));
        gE.add(new PairInt(2, 4));
        gE.add(new PairInt(4, 2));
        gE.add(new PairInt(2, 5));
        gE.add(new PairInt(5, 2));
        gE.add(new PairInt(4, 5));
        gE.add(new PairInt(5, 4));
        
        int[] expected0 = new int[]{2-1, 1-1, 4-1, 3-1, 5-1, 6-1};
        int[] expected1 = new int[]{2-1, 1-1, 4-1, 3-1, 6-1, 5-1};
        int[] expected2 = new int[]{4-1, 1-1, 2-1, 3-1, 6-1, 5-1};
        int[] expected3 = new int[]{4-1, 1-1, 2-1, 3-1, 5-1, 6-1};
        
        int[] expRCM0 = Arrays.copyOf(expected0, expected0.length);
        int[] expRCM1 = Arrays.copyOf(expected1, expected1.length);
        int[] expRCM2 = Arrays.copyOf(expected2, expected2.length);
        int[] expRCM3 = Arrays.copyOf(expected3, expected3.length);
        MiscMath0.reverse(expRCM0);
        MiscMath0.reverse(expRCM1);
        MiscMath0.reverse(expRCM2);
        MiscMath0.reverse(expRCM3);
        
        int[] rcmIdxs = CuthillMcKee.rcm(gE);

        assertTrue(Arrays.equals(expRCM0, rcmIdxs) ||
            Arrays.equals(expRCM1, rcmIdxs) || Arrays.equals(expRCM2, rcmIdxs)
            || Arrays.equals(expRCM3, rcmIdxs));
        
        int[][] a = new int[6][6];
        int i;
        for (i = 0; i < a.length; ++i) {
            a[i] = new int[6];
        }
        a[0][1] = 1;
        a[1][0] = 1;
        a[0][2] = 1;
        a[2][0] = 1;
        a[0][3] = 1;
        a[3][0] = 1;
        a[2][4] = 1;
        a[4][2] = 1;
        a[2][5] = 1;
        a[5][2] = 1;
        a[4][5] = 1;
        a[5][4] = 1;
        rcmIdxs = CuthillMcKee.rcm(gE);
        
        assertTrue(Arrays.equals(expRCM0, rcmIdxs) ||
            Arrays.equals(expRCM1, rcmIdxs) || Arrays.equals(expRCM2, rcmIdxs)
            || Arrays.equals(expRCM3, rcmIdxs));
    }
}
