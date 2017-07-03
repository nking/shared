package algorithms.misc;

import algorithms.util.PixelHelper;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Random;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class Misc0Test extends TestCase {
    
    public Misc0Test(String testName) {
        super(testName);
    }

    public void testGetSecureRandom() {
        
        Random sr = Misc0.getSecureRandom();
        
        assertNotNull(sr);
    }
    
    public void testConvertBinary() {
        
        int w = 10;
        int h = 10;
        
        TIntSet pixIdxs = new TIntHashSet();
        
        PixelHelper ph = new PixelHelper();
        
        for (int i = 0; i < w; i+=2) {
            for (int j = 1; j < h; j+=2) {
                int pixIdx = ph.toPixelIndex(i, j, w);
                pixIdxs.add(pixIdx);
            }
        }
        
        double[][] d = Misc0.convertToBinary(pixIdxs, w, h);
        
        assertNotNull(d);
        
        assertEquals(w, d.length);
        assertEquals(h, d[0].length);
        
        for (int i = 0; i < w; ++i) {
            for (int j = 1; j < h; ++j) {
                int pixIdx = ph.toPixelIndex(i, j, w);
                if (d[i][j] == 1) {
                    assertTrue(pixIdxs.contains(pixIdx));
                } else {
                    assertFalse(pixIdxs.contains(pixIdx));
                }
            }
        }
    }
}
