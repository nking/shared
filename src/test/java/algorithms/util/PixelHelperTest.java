package algorithms.util;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PixelHelperTest extends TestCase {
    
    public PixelHelperTest(String testName) {
        super(testName);
    }

    public void testToPixelIndex_PairInt_int() {
        System.out.println("toPixelIndex");
        
        int width = 10;
        int height = 4;
        PixelHelper ph = new PixelHelper();
        int[] xy = new int[2];
        
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                
                int pixIdx = ph.toPixelIndex(i, j, width);
                
                PairInt p = new PairInt(i, j);
                
                int pixIdx2 = ph.toPixelIndex(p, width);
                
                ph.toPixelCoords(pixIdx, width, xy);
                assertEquals(i, xy[0]);
                assertEquals(j, xy[1]);
                
                ph.toPixelCoords(pixIdx2, width, xy);
                assertEquals(i, xy[0]);
                assertEquals(j, xy[1]);
            }
        }
    }

    public void testConvert_Set_int() {
        System.out.println("convert");
        
        int width = 10;
        int height = 4;
        PixelHelper ph = new PixelHelper();
        int[] xy = new int[2];
        
        Set<PairInt> points = new HashSet<PairInt>();
        TIntSet pixIdxs = new TIntHashSet();
        
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                
                int pixIdx = ph.toPixelIndex(i, j, width);
                pixIdxs.add(pixIdx);
                
                PairInt p = new PairInt(i, j);
                points.add(p);
                
            }
        }
        
        TIntSet pixIdxs2 = ph.convert(points, width);
        
        Set<PairInt> points2 = ph.convert(pixIdxs, width);
        
        for (PairInt p : points2) {
            assertTrue(points.contains(p));
            assertTrue(points.remove(p));
        }
        assertTrue(points.isEmpty());
        
        TIntIterator iter2 = pixIdxs2.iterator();
        while (iter2.hasNext()) {
            int pixIdx = iter2.next();
            assertTrue(pixIdxs.contains(pixIdx));
            assertTrue(pixIdxs.remove(pixIdx));
        }
        assertTrue(pixIdxs.isEmpty());
    }
    
}
