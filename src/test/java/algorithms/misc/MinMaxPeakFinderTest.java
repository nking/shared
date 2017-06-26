package algorithms.misc;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MinMaxPeakFinderTest extends TestCase {
    
    public MinMaxPeakFinderTest(String testName) {
        super(testName);
    }
    
    public void testFindPeaks() {
        
        MinMaxPeakFinder mmpf = new MinMaxPeakFinder();
        
        float[] values = new float[] {
            0.03f, 0.05f, 0.04f, 2.4f, 0.04f, 0.02f, 6.0f 
        };
        
        int[] indexes = mmpf.findPeaks(values);
        System.out.println(Arrays.toString(indexes));
        assertEquals(2, indexes.length);
        assertEquals(3, indexes[0]);
        assertEquals(6, indexes[1]);
        
    }
}
