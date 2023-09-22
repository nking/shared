package algorithms.alignment;

import algorithms.util.PairIntArray;
import junit.framework.TestCase;

import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class KMPTest extends TestCase {

    public KMPTest() {
    }
    
    public void test0() {

        char[] p = "ababaca".toCharArray();
        char[] t = "ababababacababaca".toCharArray();
        int[] expectedPi = new int[]{0, 0, 1, 2, 3, 0, 1};
        int[] expectedM = new int[]{4, 10};

        int[] pi = KMP.computePrefixFunction(p);
        assertTrue(Arrays.equals(expectedPi, pi));

        int[] m = KMP.findPatternInText(p, t);
        assertTrue(Arrays.equals(expectedM, m));

        /*
        * m = 7
              0   1   2   3   4   5   6   7   8   9  10  11  13
         P    a   b   a   b   a   c   a
         pi   0   0   1   2   3   0   1
         T    a   b   a   b   a   b   a   b   a   c   a   b   a   b   a   c   a
         *                    *4                      *10
                  X               X
                                 pi[4]
                                 q=3                  q=1
                                 i=5                  i=10
         */

    }

    public void test1() {

        char[] p = "ababaca".toCharArray();
        char[] t = "babababacababaca".toCharArray();
        int[] expectedPi = new int[]{0, 0, 1, 2, 3, 0, 1};
        int[] expectedM = new int[]{3, 9};

        int[] pi = KMP.computePrefixFunction(p);
        assertTrue(Arrays.equals(expectedPi, pi));

        int[] m = KMP.findPatternInText(p, t);
        assertTrue(Arrays.equals(expectedM, m));

        /*
        * m = 7
              0   1   2   3   4   5   6   7   8   9  10  11  13
         P    a   b   a   b   a   c   a
         pi   0   0   1   2   3   0   1
         T    b   a   b   a   b   a   b   a   c   a   b   a   b   a   c   a
         *                *3                     *9
                                      X
                                     pi[4]
                                     q=3
         */

    }


}
