package algorithms.misc;

import algorithms.util.FormatArray;
import java.util.Arrays;
import java.util.Set;
import java.util.TreeMap;
import junit.framework.TestCase;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 *
 * @author nichole
 */
public class ShuffleTest extends TestCase {
    
    public ShuffleTest(String testName) {
        super(testName);
    }
    
    public void est0() throws NoSuchAlgorithmException, IOException {
        
        // random number generator method nextInt(bound) is not producing a uniform distribution
        //     between 0 and bound
        
        Random rand = new Random();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        int n = 7;
        
        int factor = 1000;
        int nP = (int)MiscMath0.factorial(n);
        int r;
        Map<Integer, Integer> map = new HashMap<Integer, Integer>();
        for (int i = 0; i < factor*nP; ++i) {
            r = rand.nextInt(nP);
            //r = (int)Math.round(nP*rand.nextDouble());
            //r = rand.nextInt();
            if (map.containsKey(r)) {
                map.put(r, map.get(r) + 1);
            } else {
                map.put(r, 1);
            }
        }
        System.out.printf("number of unique random keys = %d, number of draws = %d, expected %d (%.3f)\n", 
            map.size(), factor*nP, nP,  (float)map.size()/(float)nP);
        System.out.flush();
        
        int nUniqueGenerated = map.size();
        
        // calculate the median and IQR or mean and st.dev or use MAD
        //median of absolute deviation of x, median, min, and max.
        double[] count = new double[nUniqueGenerated];
        
        Set<Integer> keySet = map.keySet();
        int j = 0;
        for (int key : keySet) {
            count[j] = map.get(key);
            ++j;
        }
        
        // histogram of number of repeated permutations
        TIntIntMap countMap = new TIntIntHashMap();
        int c;
        for (j = 0; j < count.length; ++j) {
            c = (int)count[j];
            if (countMap.containsKey(c)) {
                countMap.put(c, countMap.get(c) + 1);
            } else {
                countMap.put(c, 1);
            }
        }
        int[] k2 = countMap.keys();
        int[] c2 = new int[k2.length];
        Arrays.sort(k2);
        for (j = 0; j < k2.length; ++j) {
            c2[j] = countMap.get(k2[j]);
        }
        // shows it's a poisson distribution, not a uniform distribution        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        //addPlot(double[] xPoints, double[] yPoints, double[] xPolygon, double[] yPolygon, String plotLabel)
        plotter.addPlot(k2, c2, null, null, "f=" + factor + " rand");
        plotter.writeFile("rand_hist_" + factor);
        
        Arrays.sort(count);
        
        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(count);
        double kMAD = 1.4826;
        double s = kMAD*mADMinMax[0];
        double r0 = mADMinMax[1] - 3*s;
        double r1 = mADMinMax[1] + 3*s;
        
        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(count);
        double[] avgAndStDev = MiscMath0.getAvgAndStDev(count);
                     
        System.out.printf("median of absolute deviation of x, the median, the min, and the max = \n  %s with r0=%.3f r1=%.3f\n", 
            FormatArray.toString(mADMinMax, "%.3f"), r0, r1);
        
        System.out.printf("medianAndIQR = %s\n", 
            FormatArray.toString(medianAndIQR, "%.3f"));
        
        System.out.printf("avgAndStDev = %s\n", 
            FormatArray.toString(avgAndStDev, "%.3f"));
       
    }

    public void testFisherYates_doubleArr() throws Exception {
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
                
        // not a uniform distribution
        
        // make a sequence of 7 numbers.
        //
        // run Shuffle 7! times and count the number of 5040 permutations.
       
        int[] a = new int[]{0, 1, 2, 3, 4, 5, 6};
        
        // each permutation p=1./5040
        
        long n = MiscMath0.factorial(7);
        // to see the Poisson distribution change factor to 1, 10, 100, etc or fit the profiles...
        int factor = 100;
        
        /*
        for factor=10:
        [junit] shuffled 10 * 7!
        [junit] median of absolute deviation of x, the median, the min, and the max =
        [junit]   2.000, 10.000, 1.000, 24.000  with r0=1.104 r1=18.896
        [junit] medianAndIQR = 10.000, 2.000 
        [junit] avgAndStDev = 10.002, 3.215 
        [junit] number of permutations in shuffle = 5039, expected 5040 (1.000)
        2./10 = 0.2
        
        for factor=100:
        [junit] shuffled 100 * 7!
        [junit] median of absolute deviation of x, the median, the min, and the max = 
        [junit]   7.000, 100.000, 64.000, 137.000  with r0=68.865 r1=131.135
        [junit] medianAndIQR = 100.000, 7.000 
        [junit] avgAndStDev = 100.000, 9.842 
        [junit] number of permutations in shuffle = 5040, expected 5040 (1.000)
        7.0/100 = 0.07
        
        for factor=1000:
        [junit] shuffled 1000 * 7!
        [junit] median of absolute deviation of x, the median, the min, and the max = 
        [junit]   22.000, 1000.000, 869.000, 1116.000  with r0=902.148 r1=1097.852
        [junit] medianAndIQR = 1000.000, 22.000 
        [junit] avgAndStDev = 1000.000, 31.803 
        [junit] number of permutations in shuffle = 5040, expected 5040 (1.000)
        22./1000 = 0.022
        */
        
        // sorting by key to more easily look at permutations
        TreeMap<Integer, Integer> map = new TreeMap<Integer, Integer>();

        int[] ai;
        Integer v;
        int j;
        int sum;
                
        for (int i = 0; i <= n * factor; ++i) {
            
            ai = Arrays.copyOf(a, a.length);
            
            Shuffle.fisherYates(ai, rand);
            
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
        
        // histogram of number of repeated permutations
        TIntIntMap countMap = new TIntIntHashMap();
        int c;
        for (j = 0; j < count.length; ++j) {
            c = (int)count[j];
            if (countMap.containsKey(c)) {
                countMap.put(c, countMap.get(c) + 1);
            } else {
                countMap.put(c, 1);
            }
        }
        int[] k2 = countMap.keys();
        int[] c2 = new int[k2.length];
        Arrays.sort(k2);
        for (j = 0; j < k2.length; ++j) {
            c2[j] = countMap.get(k2[j]);
        }
        // shows it's a poisson distribution, not a uniform distribution        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        //addPlot(double[] xPoints, double[] yPoints, double[] xPolygon, double[] yPolygon, String plotLabel)
        plotter.addPlot(k2, c2, null, null, "f=" + factor + " shuffle");
        plotter.writeFile("shuffle_perm_hist_" + factor);
        
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
