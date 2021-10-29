package algorithms.tsp;

import algorithms.SubsetChooser;
import algorithms.misc.Distances;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TSPDynamicOutlineTest extends TestCase {
    
    public TSPDynamicOutlineTest() {
    }
    
    public void test0() {
        
        //nP = C(n-1,k) for 1 fixed node
        int n = 10;
        int k = 3;
        SubsetChooser chooser = new SubsetChooser(n-1, k);
        int c = 0;
        long s, sInv;
        StringBuilder sb;
        //https://www.geeksforgeeks.org/program-to-invert-bits-of-a-number-efficiently/
        long all1s = 1 << (n-1);
        all1s = all1s|all1s-1;
        while (true) {
            s = chooser.getNextSubset64Bitstring();
            if (s == -1) {
                break;
            }
            
            s <<= 1;
            sb = new StringBuilder(Long.toBinaryString(s));
            
            while (sb.length() < n) {
                sb = sb.insert(0, "0");
            }
            
            sInv = s^all1s;
            sInv &= ~1; // clear bit 0
            
            System.out.printf("%d (%7s)  (%7s)\n", s, sb, Long.toBinaryString(sInv));
            c++;
        }
        System.out.printf("n=%d count=%d\n", n, c);
    }
}
