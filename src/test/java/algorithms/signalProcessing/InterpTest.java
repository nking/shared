package algorithms.signalProcessing;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class InterpTest extends TestCase {
    
    public InterpTest(String testName) {
        super(testName);
    }

    public void testBin() {
        
        System.out.println("bin");
                
        int binFactor = 2;
        
        float[] values = new float[]{2, 2, 4, 4, 2, 2};
        float[] expected = new float[]{2, 4, 2};
        
        Interp interp = new Interp();
        float[] binned = interp.bin(values, binFactor);
        
        assertTrue(Arrays.equals(expected, binned));
        
        float[] unbinned = interp.unbin(binned, binFactor);
        
        assertTrue(Arrays.equals(unbinned, values));
    }

    public void testLinearInterp() {
        
        float scale = 3.333f;
    
        Interp interp = new Interp();
        
        float[] values = new float[10];
        for (int i = 0; i < values.length; ++i) {
            values[i] = i;
        }
        
        int outLength = (int)Math.ceil(scale * values.length);
        
        float[] r = interp.linearInterp(
            values, outLength, 0.f, 255.f);
        
        assertEquals(outLength, r.length);
        
        System.out.println(Arrays.toString(values));
        System.out.println(Arrays.toString(r));
        
    }
    
}
