package algorithms.misc;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MiscSorterTest extends TestCase {
    
    public void testSortBy1stArgThen2nd(){

        /*   10
         *
         *                *p12
         *    9
         *
         *
         *    8
         *
         *
         *    7
         *
         *
         *    6                             *p10 *p11
         *
         *                   *p9       *p8
         *    5         *p7
         *
         *                        *p6
         *    4                                            *p5        *p4
         *
         *     p3*
         *    3                                                 *p2
         *
         *
         *    2                                                   *p1
         *
         *
         *    1            *p0
         *
         *
         *    0    1    2    3    4    5    6    7    8    9   10   11   12
    	 */
    	float[] x = new float[13];
    	float[] y = new float[13];
        x[0] =  2.5f;
        y[0] =  1.0f;
        x[1] = 10.5f;
        y[1] =  2.0f;
        x[2] = 10.0f;
        y[2] =  3.0f;
        x[3] =  0.5f;
        y[3] =  3.2f;
        x[5] =  9.0f;
        y[5] =  4.0f;

        x[4] = 11.2f;
        y[4] =  4.0f;
        x[6] =  4.0f;
        y[6] =  4.2f;
        x[7] =  2.0f;
        y[7] =  5.0f;
        x[9] =  3.0f;
        y[9] =  5.3f;
        x[8] =  5.0f;
        y[8] =  5.3f;

        x[10] =  6.0f;
        y[10] =  6.0f;
        x[11] =  7.0f;
        y[11] =  6.0f;
        x[12] =  2.3f;
        y[12] =  9.3f;


        float[] ex = Arrays.copyOf(x, x.length);
    	float[] ey = Arrays.copyOf(y, y.length);
        ex[4] =  x[5];
        ey[4] =  y[5];
        ex[5] =  x[4];
        ey[5] =  y[4];
        
        ex[8] =  x[9];
        ey[8] =  y[9];
        ex[9] =  x[8];
        ey[9] =  y[8];

    	MiscSorter.sortBy1stArgThen2nd(y, x);
    	assertTrue(x.length == ex.length);

    	for (int i=0; i < ex.length; i++) {
            assertTrue( Math.abs(ex[i] - x[i]) < 0.01);
            assertTrue( Math.abs(ey[i] - y[i]) < 0.01);
        }
    	
    	boolean caughtException = true;
        try {
            float[] x2 = null;
            float[] y2 = null;
            MiscSorter.sortBy1stArgThen2nd(x2, y2);
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            float[] y2 = null;
            MiscSorter.sortBy1stArgThen2nd(new float[]{1f,2f}, y2);
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            MiscSorter.sortBy1stArgThen2nd(new float[]{1f,2f}, new float[]{1f,2f,3f});
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
    }
    
     public void testSortBy1stArgThen2nd_int(){

        /*   10
         *
         *                *p12
         *    9
         *
         *
         *    8
         *
         *
         *    7
         *
         *
         *    6                             *p10 *p11
         *
         *                   *p9       *p8
         *    5         *p7
         *
         *                        *p6
         *    4                                            *p5        *p4
         *
         *     p3*
         *    3                                                 *p2
         *
         *
         *    2                                                   *p1
         *
         *
         *    1            *p0
         *
         *
         *    0    1    2    3    4    5    6    7    8    9   10   11   12
    	 */
    	int[] x = new int[13];
    	int[] y = new int[13];
        x[0] =  25;
        y[0] =  10;
        x[1] = 105;
        y[1] =  20;
        x[2] = 100;
        y[2] =  30;
        x[3] =  5;
        y[3] =  32;
        x[5] =  90;
        y[5] =  40;

        x[4] = 112;
        y[4] =  40;
        x[6] =  40;
        y[6] =  42;
        x[7] =  20;
        y[7] =  50;
        x[9] =  30;
        y[9] =  53;
        x[8] =  50;
        y[8] =  53;

        x[10] =  60;
        y[10] =  60;
        x[11] =  70;
        y[11] =  60;
        x[12] =  23;
        y[12] =  93;


        int[] ex = Arrays.copyOf(x, x.length);
    	int[] ey = Arrays.copyOf(y, y.length);
        ex[4] =  x[5];
        ey[4] =  y[5];
        ex[5] =  x[4];
        ey[5] =  y[4];
        
        ex[8] =  x[9];
        ey[8] =  y[9];
        ex[9] =  x[8];
        ey[9] =  y[8];

    	MiscSorter.sortBy1stArgThen2nd(y, x);
    	assertTrue(x.length == ex.length);

    	for (int i=0; i < ex.length; i++) {
            assertTrue( Math.abs(ex[i] - x[i]) < 0.01);
            assertTrue( Math.abs(ey[i] - y[i]) < 0.01);
        }
    	
    	boolean caughtException = true;
        try {
            int[] a0 = null;
            int[] a1 = null;
            MiscSorter.sortBy1stArgThen2nd(a0, a1);
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            int[] a1 = null;
            MiscSorter.sortBy1stArgThen2nd(new int[]{1, 2}, a1);
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            MiscSorter.sortBy1stArgThen2nd(new int[]{1, 2}, new int[]{1, 2, 3});
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
    }
    
    public void testSortByFirstArgument_1() {

        int[] a = new int[]{3, 7, 1, 5};
        int[] b = new int[]{0, 1, 2, 3};
        
        MiscSorter.sortBy1stArg(a, b);
        
        assertTrue(a[0] == 1);
        assertTrue(b[0] == 2);
        
        assertTrue(a[1] == 3);
        assertTrue(b[1] == 0);
        
        assertTrue(a[2] == 5);
        assertTrue(b[2] == 3);
        
        assertTrue(a[3] == 7);
        assertTrue(b[3] == 1); 
    }
 
    public void testSortByFirstArgument_1_mergesort() {

        int[] a = new int[]{3, 7, 1, 5};
        int[] b = new int[]{0, 1, 2, 3};
        
        MiscSorter.sortBy1stArg2(a, b);
        
        assertTrue(a[0] == 1);
        assertTrue(b[0] == 2);
        
        assertTrue(a[1] == 3);
        assertTrue(b[1] == 0);
        
        assertTrue(a[2] == 5);
        assertTrue(b[2] == 3);
        
        assertTrue(a[3] == 7);
        assertTrue(b[3] == 1); 
    }
 
    public void testSortByFirstArgument_int_float() {

        int[] a = new int[]{3, 7, 1, 5};
        float[] b = new float[]{0, 1, 2, 3};
        
        MiscSorter.sortBy1stArg(a, b);
        
        assertTrue(a[0] == 1);
        assertTrue(b[0] == 2.f);
        
        assertTrue(a[1] == 3);
        assertTrue(b[1] == 0.f);
        
        assertTrue(a[2] == 5);
        assertTrue(b[2] == 3.f);
        
        assertTrue(a[3] == 7);
        assertTrue(b[3] == 1.f); 
    }
    
    public void testSortByDecr() throws Exception {
        
        int[] a = new int[]{1, 2, 3, 4, 5, 6};
    	int[] b = new int[]{0, 1, 2, 3, 4, 5};

    	MiscSorter.sortByDecr(a, b);
    	assertTrue(a.length == b.length);

    	int[] expectedA = new int[]{6, 5, 4, 3, 2, 1};
        int[] expectedB = new int[]{5, 4, 3, 2, 1, 0};
        
        assertTrue(Arrays.equals(expectedA, a));
        assertTrue(Arrays.equals(expectedB, b));
    }
    
    public void testSortByDecr_float_int() throws Exception {
        
        float[] a = new float[]{1, 2, 3, 4, 5, 6};
    	int[] b = new int[]{0, 1, 2, 3, 4, 5};

    	MiscSorter.sortByDecr(a, b);
    	assertTrue(a.length == b.length);

    	float[] expectedA = new float[]{6, 5, 4, 3, 2, 1};
        int[] expectedB = new int[]{5, 4, 3, 2, 1, 0};
        
        assertTrue(Arrays.equals(expectedA, a));
        assertTrue(Arrays.equals(expectedB, b));
    }
    
    public void testSortBy1stArgDecrThen2ndIncr() {
        
        int[] a = new int[] {
            1, 2, 2, 1 
        };
        int[] b = new int[] {
            1, 2, 3, 1
        };
        
        int[] ea = new int[] {
            2, 2, 1, 1
        };
        int[] eb = new int[] {
            2, 3, 1, 1
        };
        
        MiscSorter.sortBy1stArgDecrThen2ndIncr(a, b);
        
        assertTrue(Arrays.equals(ea, a));
        assertTrue(Arrays.equals(eb, b));
        
    }
    
    public static void testMergeBy1stArgThen2nd() {
        
        int i = 0;
        double eps = 1.e-17;
        double dx, dy;
        
        double[] x = new double[]{4,   5,   4,   1,   2, 8};
        double[] y = new double[]{3.1, 1.1, 3, 1.1, 1.0, 1.1};
        int[] indexes = new int[x.length];
        for (i = 0; i < x.length; ++i) {
            indexes[i] = i;
        }
        
        double[] eX = new double[]{1,     2,  4,   4,   5,   8};
        double[] eY = new double[]{1.1, 1.0,  3, 3.1, 1.1, 1.1};
        int[] eIndexes = new int[]{3, 4, 2, 0, 1, 5};
                
        int[] sortedIndexes = MiscSorter.mergeBy1stArgThen2nd(x, y);
        
        assertEquals(eIndexes.length, sortedIndexes.length);
                
        for (i = 0; i < indexes.length; ++i) {
            assertEquals(eIndexes[i], sortedIndexes[i]);
            assertEquals(eIndexes[i], indexes[ eIndexes[i] ] );
            
            dx = Math.abs(eX[i] - x[i]);
            dy = Math.abs(eY[i] - y[i]);
            
            assertTrue(dx < eps);
            assertTrue(dy < eps);
        }
    }
}
