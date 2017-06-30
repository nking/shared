package algorithms.misc;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FrequencyTest extends TestCase {
    
   public void test0() {
       
       Frequency f = new Frequency();
       
       int[][] d = getData0();
       
       TIntList values = new TIntArrayList();
       TIntList counts = new TIntArrayList();
       f.calcFrequency(d, values, counts, true);
       
       int[] expectedV = new int[]{1, 2, 4, 8};
       int[] expectedC = new int[]{12, 4, 4, 1};
       
       assertEquals(expectedV.length, values.size());
       assertEquals(expectedC.length, counts.size());
       
       for (int i = 0; i < expectedV.length; ++i) {
           assertEquals(expectedV[i], values.get(i));
           assertEquals(expectedC[i], counts.get(i));
       }
   }
   
   private int[][] getData0() {
       
       /*
        row 0:  0 0 0 0 0 0 0 0 0 
        row 1:  0 0 0 0 0 0 0 0 0 
        row 2:  0 0 0 1 1 1 0 0 0 
        row 3:  0 0 1 2 4 2 1 0 0 
        row 4:  0 0 1 4 8 4 1 0 0 
        row 5:  0 0 1 2 4 2 1 0 0 
        row 6:  0 0 0 1 1 1 0 0 0 
        row 7:  0 0 0 0 0 0 0 0 0 
        row 8:  0 0 0 0 0 0 0 0 0
       */
       
       int[][] data = new int[9][];
       data[0] = new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0};
       data[1] = new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0};
       data[2] = new int[]{0, 0, 0, 1, 1, 1, 0, 0, 0};
       data[3] = new int[]{0, 0, 1, 2, 4, 2, 1, 0, 0};
       data[4] = new int[]{0, 0, 1, 4, 8, 4, 1, 0, 0};
       data[5] = new int[]{0, 0, 1, 2, 4, 2, 1, 0, 0};
       data[6] = new int[]{0, 0, 0, 1, 1, 1, 0, 0, 0};
       data[7] = new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0};
       data[8] = new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0};
   
       return data;
   }
}
