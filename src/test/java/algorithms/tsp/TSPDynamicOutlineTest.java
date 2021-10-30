package algorithms.tsp;

import algorithms.SubsetChooser;
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
        int nSetBits;
        long tZeros;
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
            
            nSetBits = Long.bitCount(sInv);
            tZeros = Long.numberOfTrailingZeros(sInv);
            
            System.out.printf("%d (%7s)  (%7s) nSetBits=%d tz=%d\n", 
                s, sb, Long.toBinaryString(sInv), nSetBits, tZeros);
            
            // this permutation of nSetBits can be used for all 84 of the first
            //   3-node paths.
            //   each use of the permutation can be left shifted by tZeros
            //   then assigned positions relative to the set bits in sInv
            //   to get the permuation for this specific si
            SubsetChooser chooseri = new SubsetChooser(nSetBits,  k);
            long si;
            while (true) {
                si = chooseri.getNextSubset64Bitstring();
                if (si == -1) {
                    break;
                }
                long s2 = si << tZeros;
                
                // the nSetBits in si are to be shifted by tZeros
                // then converted to the si set bit positions
                //
                
            }
            
            c++;
        }
        System.out.printf("n=%d count=%d\n", n, c);
    }
}
