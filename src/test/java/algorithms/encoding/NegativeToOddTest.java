package algorithms.encoding;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class NegativeToOddTest extends TestCase {
    
    public NegativeToOddTest(String testName) {
        super(testName);
    }

    /**
     * Test of negativeToOddEncode method, of class NegativeToOdd.
     */
    public void testNegativeToOddEncode() {
        System.out.println("negativeToOddEncode");
        int[] a0 = new int[]{0, -1, 4};
        int[] a = Arrays.copyOf(a0, a0.length);
        int[] expected = new int[]{0, 3, 8};
        
        NegativeToOdd.negativeToOddEncode(a);
        assertTrue(Arrays.equals(expected, a));
      
        NegativeToOdd.negativeToOddDecode(a);
        
        assertTrue(Arrays.equals(a0, a));
    }
    
}
