package algorithms.signalProcessing;

import algorithms.misc.MiscMath0;
import static junit.framework.Assert.assertEquals;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class B3SplineFunction1DTest extends TestCase {
    
    public B3SplineFunction1DTest(String testName) {
        super(testName);
    }

    public void testCalculate() {
        
        /*
        1/16, 1/4, 3/8, 1/4, 1/16
        */
        
        float[] a = new float[] {
            2, 2, 4, 6, 8, 10, 33, 10, 8, 6, 4, 2, 2
        };
        
        B3SplineFunction1D sf = new B3SplineFunction1D();
        
        float[] b = sf.calculate(a);
        //System.out.println("a=" + Arrays.toString(a));
        //System.out.println("b=" + Arrays.toString(b));
        
        int maxIdx0 = MiscMath0.findYMaxIndex(a);
        int maxIdx1 = MiscMath0.findYMaxIndex(b);

        assertEquals(maxIdx0, maxIdx1);        
        
        for (int i = 2; i < a.length - 2; ++i) {
            float v0 = 
                (1.f/16.f)*a[i-2] + (1.f/4.f)*a[i-1]
                + (3.f/8.f)*a[i] + (1.f/4.f)*a[i+1] + (1.f/16.f)*a[i+2];
            float v1 = b[i];
            float d = Math.abs(v1 - v0);
            assertEquals(0.f, d);
            //System.out.println("d=" + d);
        }
    }
    
}
