package algorithms.signalProcessing;

import algorithms.misc.MiscMath0;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TriangleFunction1DTest extends TestCase {
    
    public TriangleFunction1DTest(String testName) {
        super(testName);
    }
    
    public void testCalculateNextLevel() {
        
        //public float[] calculateNextLevel(float[] input, int j) {
        
        //c_(j + 1,k) = (1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) + (1/4)*c_(j,k+(2^j))
        
        float[] a = new float[] {
            2, 2, 4, 6, 8, 10, 12, 10, 8, 6, 4, 2, 2
        };
        
        TriangleFunction1D sf = new TriangleFunction1D();
        
        float[] b = sf.calculateNextLevel(a, 0);
        //System.out.println("a=" + Arrays.toString(a));
        //System.out.println("b=" + Arrays.toString(b));
        
        int maxIdx0 = MiscMath0.findYMaxIndex(a);
        int maxIdx1 = MiscMath0.findYMaxIndex(b);

        assertEquals(maxIdx0, maxIdx1);        
        
        for (int i = 1; i < maxIdx0 - 1; ++i) {
            float v0 = 0.25f*a[i-1] + 0.5f*a[i] + 0.25f*a[i+1];
            float v1 = b[i];
            float d = Math.abs(v1 - v0);
            assertEquals(0.f, d);
            //System.out.println("d=" + d);
        }
    }
    
}
