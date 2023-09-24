package algorithms.sort;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import junit.framework.TestCase;

public class BucketSortTest extends TestCase {

    public void testsort() throws Exception {
        
        int[] a = new int[]{2, 5, 3, 0, 2, 3, 0, 3};
        
        int[] indexes = BucketSort.sortAndReturnIndexes(a);
        
        int[] expected = new int[]{0, 0, 2, 2, 3, 3, 3, 5};
        
        assertTrue(Arrays.equals(expected, a));

        double[] a2 = new double[]{2, 5, 3, 0, 2, 3, 0, 3};
        int[] indexes2 = BucketSort.sortAndReturnIndexes(a2);
        assertTrue(Arrays.equals(indexes, indexes2));
        for (int i = 0; i < a.length; ++i){
            assertTrue(Math.abs(a2[i] - expected[i]) < 1E-7);
        }
    }

    public void testsort3() throws Exception {
        
        // use more than 46340 random numbers whose value is higher than
        // 46340 to show that the internal long summations are safely
        // reduced back to the integer values
        
        List<Integer> list = new ArrayList<Integer>();
        
        SecureRandom sr = new SecureRandom();
        long seed = System.nanoTime();
        sr.setSeed(seed);
        
        int n = 46340*2;
        
        int[] a = new int[n];
        double[] a2 = new double[n];
        
        int max = Integer.MIN_VALUE;
        
        for (int i = 0; i < n; i++) {
            int r = sr.nextInt(67000);
            while (r < 0) {
                r = sr.nextInt();
            }
            list.add(Integer.valueOf(r));
            a[i] = r;
            a2[i] = r;
            
            if (r > max) {
                max = r;
            }
        }

        Collections.sort(list);
        
        int[] indexes = BucketSort.sortAndReturnIndexes(a);
        for (int i = 0; i < n; i++) {
            assertTrue(list.get(i).intValue() == a[i]);
        }

        int[] indexes2 = BucketSort.sortAndReturnIndexes(a2);
        for (int i = 0; i < n; i++) {
            assertTrue(Math.abs(list.get(i).intValue() - a2[i]) < 1E-7);
        }

        assertTrue(Arrays.equals(indexes, indexes2));
    }

}
