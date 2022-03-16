package algorithms.scheduling;

import algorithms.util.FormatArray;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MiscTest extends TestCase {
    
    public MiscTest(String testName) {
        super(testName);
    }
    
    public void testUnweightedIntervalMinimizeLateGreedy() {
        double[] d = new double[]{1, 2, 3, 4, 5};
        double[] t = new double[]{0.5, 0.75, 1.8, 1.8, 0.3};

        double[] outputStart = new double[t.length];
        double[] outputLate = new double[t.length];
        Misc misc = new Misc();
        int[] indexes = misc.unweightedIntervalMinimizeLateGreedy(t, d, outputStart, outputLate);
        
        /*
        System.out.printf("lateness=\n%s\n", FormatArray.toString(outputLate, "%.3f"));
        System.out.printf("s=\n%s\n", FormatArray.toString(outputStart, "%.3f"));
        System.out.printf("indexes=\n%s\n", Arrays.toString(indexes));
        */
        
        int[] expectedOrder = new int[]{0, 1, 2, 3, 4};
        double[] expectedLate = new double[]{0.000, 0.000, 0.050, 0.850, 0.150};
        double[] expectedStart = new double[]{0.000, 0.500, 1.250, 3.050, 4.850};
        
        assertTrue(Arrays.equals(expectedOrder, indexes));
        double eps = 1e-7;
        double diff;
        int i;
        for (i = 0; i < expectedLate.length; ++i) {
            diff = Math.abs(expectedLate[i] - outputLate[i]);
            assertTrue(diff < eps);
            diff = Math.abs(expectedStart[i] - outputStart[i]);
            assertTrue(diff < eps);
        }
    }
    
     /**
     * test follows the example in the lecture notes of David Mount for CMSC 451 
     * Design and Analysis of Computer Algorithms.
     * https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect10-dp-intv-sched.pdf
     * 
     */
    public void testWeightedIntervalBottomUp() {
        
        //TODO: add more tests and include changing the order of items in this one
        
        //                              *               *
        double[] s = new double[]{0.5,  1, 2.75, 1.25, 5.25,   6};
        double[] f = new double[]{2.5,  3,    5,    7,  7.5, 8.5};
        double[] v = new double[]{  2,  6,  3.5,    7,    8, 1.1};

        Misc misc = new Misc();
        int[] indexes = misc.weightedIntervalBottomUp(s, f, v);
        //System.out.println("scheduled intervals = " + Arrays.toString(indexes));
        int i;
        double sum = 0;
        for (i = 0; i < indexes.length; ++i) {
            //System.out.printf("%d [%.2f : %.2f]  sum=%.2f\n", indexes[i], s[indexes[i]], f[indexes[i]], sum);
            sum += v[indexes[i]];
        }
        //System.out.println("sum of values=" + sum);
        assertEquals(14.0, sum);
     
        int[] expected = new int[]{4, 1};
        assertTrue(Arrays.equals(expected, indexes));

        
        // change the order  
        s = new double[]{1.25, 0.5,  5.25,  1, 2.75,   6};
        f = new double[]{   7, 2.5,   7.5,  3,    5, 8.5};
        v = new double[]{   7,   2,     8,  6,  3.5, 1.1};
        
        misc = new Misc();
        indexes = misc.weightedIntervalBottomUp(s, f, v);
        //System.out.println("scheduled intervals = " + Arrays.toString(indexes));
        sum = 0;
        for (i = 0; i < indexes.length; ++i) {
            //System.out.printf("%d [%.2f : %.2f]  sum=%.2f\n", indexes[i], s[indexes[i]], f[indexes[i]], sum);
            sum += v[indexes[i]];
        }
        //System.out.println("sum of values=" + sum);
        assertEquals(14.0, sum);
     
        expected = new int[]{2, 3};
        assertTrue(Arrays.equals(expected, indexes));
        
        // add a redundant interval with lower weight  
        s = new double[]{1.25, 0.5,  5.25,  1, 2.75,   6, 1};
        f = new double[]{   7, 2.5,   7.5,  3,    5, 8.5, 3};
        v = new double[]{   7,   2,     8,  6,  3.5, 1.1, 1};
        
        misc = new Misc();
        indexes = misc.weightedIntervalBottomUp(s, f, v);
        //System.out.println("scheduled intervals = " + Arrays.toString(indexes));
        sum = 0;
        for (i = 0; i < indexes.length; ++i) {
            //System.out.printf("%d [%.2f : %.2f]  sum=%.2f\n", indexes[i], s[indexes[i]], f[indexes[i]], sum);
            sum += v[indexes[i]];
        }
        //System.out.println("sum of values=" + sum);
        assertEquals(14.0, sum);
     
        expected = new int[]{2, 3};
        assertTrue(Arrays.equals(expected, indexes));
    }
    
    public void testWeightedGreedy() {
        System.out.println("unit task with weighted deadline scheduling");
        
        // testing for task scheduling with penalty for lateness
        
        // from cormen et al. Fig 16.7
        /*
         * Fig 16.7: 
         *                      0   1   2   3   4   5   6
         *                     a1  a2  a3  a4  a5  a6  a7
           int[] d = new int[]{4,  2,  4,   3, 1,  4,   6};
           int[] w = new int[]{70, 60, 50, 40, 30, 20, 10};

         *     greedy algorithm selects a1, a2, a3, a4   => (4,70),(2,60),(4,50),(3,40)
         *                      rejects a5, a6           => (1,30),(4,20)
         *                      accepts a7               => (6,10)
         *     final optimal schedule is a2, a4, a1, a3, a7, a5, a6
         *       w/ total penalty of w5+w6=50
         *      (2,60), (3,40), (4,70), (4,50), (6,10), (1,30), (4,20)
         * 
        */
        int[] d = new int[]{4,  2,  4,   3, 1,  4,   6};
        int[] w = new int[]{70, 60, 50, 40, 30, 20, 10};
       
        Misc misc =new Misc();
        int[] indexes = misc.weightedGreedy(d, w);
        
        /*
        System.out.println("scheduled tasks = " + Arrays.toString(indexes));
        int t, idx, fi;
        int sum = 0;
        // tasks start at time=0, are 1 unit in duration, and follow one another without gaps in time.
        for (t = 0; t < d.length; ++t) {
            // finish time of task i
            fi = t + 1;
            idx = indexes[t];
            if (d[idx] < fi) {
                sum += w[idx];
            }
            System.out.printf("a%d (%d,%d)\n", idx+1, d[idx], w[idx]);
        }
        System.out.printf("\nsum of late tasks=%d\n", sum);
        */
        int[] expected = new int[]{1, 3, 0, 2, 6, 5, 4};
        assertTrue(Arrays.equals(expected, indexes));
    }
   
    public void testUnweightedIntervalNoConflicts() {
        
        double[] s = new double[]{1,4,0,7,2, 8,11,3};
        double[] t = new double[]{5, 6, 9, 10, 11, 13, 14, 15};
        int[] expected = new int[]{1-1, 4-1, 7-1};
        Misc misc = new Misc();
        int[] indexes = misc.unweightedIntervalNoConflicts(s, t);
        assertTrue(Arrays.equals(expected, indexes));
    }
    
    public void testIntervalPartitionGreedy() {
        
        /* Fig 4. from
        lecture 7 notes of David Mount for CMSC 451       
        Design and Analysis of Computer Algorithms (with some corrections for pseudocode indexes).
        https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect07-greedy-sched.pdf
        */
        //                        0   1   2    3   4    5   6   7   8   9   10   11
        double[] s = new double[]{0,  1,  2,  6.5, 7.5, 10, 13, 19, 23, 26, 26,  31};
        double[] f = new double[]{4.5,10, 6.5,11,  22,  24,17.5,25, 32, 29.5, 33, 34};
        
        Misc misc = new Misc();
        int[] resources = misc.intervalPartitionGreedy(s, f);
        System.out.println("resources=" + Arrays.toString(resources));
        
        //                        0   1  2  3  4  5  6  7  8  9  10 11
        int[] expected = new int[]{0, 1, 2, 0, 2, 1, 0, 0, 2, 0, 1, 0};
        int[] expected2 = new int[]{0, 1, 2, 0, 2, 1, 0, 0, 2, 1, 0, 1};
        assertTrue(Arrays.equals(expected, resources) || 
            Arrays.equals(expected2, resources));
    }
}
