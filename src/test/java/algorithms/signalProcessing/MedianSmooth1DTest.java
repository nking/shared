package algorithms.signalProcessing;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MedianSmooth1DTest extends TestCase {
    
    public MedianSmooth1DTest(String testName) {
        super(testName);
    }
        
    public void testCalculate0() throws Exception {
        
        /*
        k=3
        curve=10 4's
        expectedN = 10 - kPoints + 1
        4 4 4 4 4 4 4 4 4 4
            4 4 4 4 4 4 4 4
            0 1 2 3 4 5 6 7
        */
        MedianSmooth1D interp = new MedianSmooth1D();
        
        float[] curveY = new float[]{4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
        int kPoints = 3;
        float[] result = interp.calculate(curveY, kPoints);
                
        float[] expected = new float[]{4, 4, 4, 4, 4, 4, 4, 4};
        
        assertTrue(result.length == expected.length);
        
        assertTrue(Arrays.equals(result, expected));
        
    } 
      
    public void testCalculate1() throws Exception {
        
        /*
        k=3
        curve=10 4's
        expectedN = 10 - kPoints + 1
        2 3 4 3 2 4 4 4 4 4
            3 3 3 3 4 4 4 4
            0 1 2 3 4 5 6 7
        */
        MedianSmooth1D interp = new MedianSmooth1D();
        
        float[] curveY = new float[]{2, 3, 4, 3, 2, 4, 4, 4, 4, 4};
        int kPoints = 3;
        float[] result = interp.calculate(curveY, kPoints);
                
        float[] expected = new float[]{3, 3, 3, 3, 4, 4, 4, 4};
        
        assertTrue(result.length == expected.length);
        
        assertTrue(Arrays.equals(result, expected));
        
        
        curveY = null;
        boolean caughtException = false;
        try {
            result = interp.calculate(curveY, kPoints);
        } catch (Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
    }
    
}
