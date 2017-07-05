package algorithms.imageProcessing;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FiltersTest extends TestCase {
    
    public FiltersTest(String testName) {
        super(testName);
    }

    public void testMaximumFilter_floatArrArr_int() {
        //public static float[][] maximumFilter(float[][] img, int size)
    }

    public void testMaximumFilter_floatArr_int() {
       //public static float[] maximumFilter(float[] img, int size)
    }

    public void testPeakLocalMax_5args() {
       /*
        public void peakLocalMax(float[][] img, int minDistance,
            float thresholdRel,
            TIntList outputKeypoints0, TIntList outputKeypoints1)
        */
    }

    public void testPeakLocalMax_4args() {
        
        //public void peakLocalMax(float[] img, int minDistance, float thresholdRel,
        //    TIntList outputKeypoints)
        
        float[] a = new float[]{10, 10, 10, 100, 20, 25, 14, 10, 15, 100, 13, 13, 13};
        
        Filters filters = new Filters();
        
        TIntList output = new TIntArrayList();
        
        filters.peakLocalMax(a, 0, 0.85f, output);
        
        assertEquals(2, output.size());
        assertEquals(3, output.get(0));
        assertEquals(9, output.get(1));
        
    }
    
}
