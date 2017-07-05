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
    
    public void testCalculateMeanOfSmallest() {
                
        int n = 100;
        int[] a = new int[n];
        Arrays.fill(a, 100);
        int n0 = (int)(0.03f*n);
        for (int i = 0; i < n0; ++i) {
            a[i] = 10;
        }
        
        MinMaxPeakFinder finder = new MinMaxPeakFinder();
        float avgMin = finder.calculateMeanOfSmallest(a, 0.03f);
        assertEquals(10.f, avgMin);
    }
    
    public void testCalculateMeanOfSmallest_float() {
                
        int n = 100;
        float[] a = new float[n];
        Arrays.fill(a, 100);
        int n0 = (int)(0.03f*n);
        for (int i = 0; i < n0; ++i) {
            a[i] = 10;
        }
        
        MinMaxPeakFinder finder = new MinMaxPeakFinder();
        float avgMin = finder.calculateMeanOfSmallest(a, 0.03f);
        assertEquals(10.f, avgMin);
    }
    
}
