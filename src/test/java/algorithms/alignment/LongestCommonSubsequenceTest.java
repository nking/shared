package algorithms.alignment;

import java.util.Arrays;

import algorithms.util.PairIntArray;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LongestCommonSubsequenceTest extends TestCase {
    
    public LongestCommonSubsequenceTest() {
    }
    
    public void testCalculateWithCormenEtAl() {
        System.out.println("testCalculateWithCormenEtAl");
        
        // Fig. 15.6 from Cormen et al.
        
        /*
         *  C  j  .  0  1   2   3   4   5   6
         *        .
         *  i  B  . y_i B   D   C   A   B   A
         *  ......---------------------------------
         *  0 x_i | 0   0   0   0   0   0   0
         *  1  A  | 0
         *  2  B  | 0
         *  3  C  | 0
         *  4  B  | 0
         *  5  D  | 0
         *  6  A  | 0
         *  7  B  | 0
         *
         * for X=ABCBDAB
         * and Y=BDCABA
         *
         *  C  j  .  0   1    2    3    4    5    6
         *        .
         *  i  B  . y_i  B    D    C    A    B    A
         *  ......---------------------------------
         *  0 x_i | 0    0    0    0    0    0    0
         *  1  A  |*0   U0   U0   U0  UL1   L1   UL1
         *  2  B  | 0 *UL1  *L1   L1   U1  UL2    L2
         *  3  C  | 0   U1   U1 *UL2  *L2   U2    U2
         *  4  B  | 0  UL1   U1   U2   U2 *UL3    L3
         *  5  D  | 0   U1  UL2   U2   U2  *U3    U3
         *  6  A  | 0   U1   U2   U2  UL3   U3  *UL4
         *  7  B  | 0  UL1   U2   U2   U3  UL4   *U4
        */
        
        char[] x = (new String("ABCBDAB")).toCharArray();
        char[] y = (new String("BDCABA")).toCharArray();
        
        char[] ans = (new String("BCBA")).toCharArray();
        char[] result = LongestCommonSubsequence.calculateWithCormenEtAl(x, y, true);
        
        System.out.println("result=" + new String(result));
        
        System.out.println("\nCompare to answer:");
        System.out.println("   0   0   0   0   0   0   0");
        System.out.println("   0  U0  U0  U0 UL1  L1 UL1");
        System.out.println("   0 UL1  L1  L1  U1 UL2  L2");
        System.out.println("   0  U1  U1 UL2  L2  U2  U2");
        System.out.println("   0 UL1  U1  U2  U2 UL3  L3");
        System.out.println("   0  U1 UL2  U2  U2  U3  U3");
        System.out.println("   0  U1  U2  U2 UL3  U3 UL4");
        System.out.println("   0 UL1  U2  U2  U3 UL4  U4");
        
        assertTrue(Arrays.equals(ans, result));

        PairIntArray out = new PairIntArray();
        StringEditDistance sed = new StringEditDistance();
        int r = sed.calculateWithWagnerFischer(new String("ABCBDAB"),
                new String("BDCABA"), out, false);
        System.out.println("r=" + r + " out=" + out.toString());
        PairIntArray expected0 = new PairIntArray();
        expected0.add(5,5);
        expected0.add(3,4);
        expected0.add(2,2);
        expected0.add(1,0);
        assertEquals(expected0.toString(), out.toString());
    }


}
