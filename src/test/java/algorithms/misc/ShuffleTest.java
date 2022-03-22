package algorithms.misc;

import algorithms.util.FormatArray;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ShuffleTest extends TestCase {
    
    public ShuffleTest(String testName) {
        super(testName);
    }

    public void testFisherYates_doubleArr() throws Exception {
        
        // make a sequence of 7 numbers.
        //
        // run Shuffle 7! times and count the number of 5040 permutations.
       
        int[] a = new int[]{0, 1, 2, 3, 4, 5, 6};
        
        // each permutation p=1./5040
        
        long n = MiscMath0.factorial(7);
        int factor = 10;
        
        // sorting by key to more easily look at permutations
        TreeMap<Integer, Integer> map = new TreeMap<Integer, Integer>();

        int[] ai;
        Integer v;
        int j;
        int sum;
                
        for (int i = 0; i <= n * factor; ++i) {
            
            ai = Arrays.copyOf(a, a.length);
            
            Shuffle.fisherYates(ai);
            
            // store each number in decimal digit place
            sum = 0;
            for (j = 0; j < ai.length; ++j) {
                sum += (ai[j] * Math.pow(10, ai.length - j - 1));
            }
            
            v = map.get(sum);
            if (v == null) {
                map.put(sum, 1);
            } else {
                map.put(sum, v + 1);
            }
        }
        
        int nUniqueGenerated = map.size();
        
        // calculate the median and IQR or mean and st.dev or use MAD
        //median of absolute deviation of x, median, min, and max.
        double[] count = new double[nUniqueGenerated];
        
        Set<Integer> keySet = map.keySet();
        j = 0;
        for (int key : keySet) {
            count[j] = map.get(key);
            ++j;
        }
        Arrays.sort(count);
        
        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(count);
        double kMAD = 1.4826;
        double s = kMAD*mADMinMax[0];
        double r0 = mADMinMax[1] - 3*s;
        double r1 = mADMinMax[1] + 3*s;
        
        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(count);
        double[] avgAndStDev = MiscMath0.getAvgAndStDev(count);
        
        System.out.printf("shuffled %d * %d!\n", factor, a.length);
        
        System.out.printf("median of absolute deviation of x, the median, the min, and the max = \n  %s with r0=%.3f r1=%.3f\n", 
            FormatArray.toString(mADMinMax, "%.3f"), r0, r1);
        
        System.out.printf("medianAndIQR = %s\n", 
            FormatArray.toString(medianAndIQR, "%.3f"));
        
        System.out.printf("avgAndStDev = %s\n", 
            FormatArray.toString(avgAndStDev, "%.3f"));
       
        /*System.out.printf("number of permutations in shuffle = %d, expected %d (%.3f)\n", 
            map.size(), n, (float)map.size()/(float)n);
        keySet = map.keySet();
        for (int key : keySet) {
            v = map.get(key);
            System.out.printf("%7d (%d)\n", key, v);
        }*/
        
        System.out.printf("number of permutations in shuffle = %d, expected %d (%.3f)\n", 
            map.size(), n, (float)map.size()/(float)n);
        System.out.flush();
    }

}
