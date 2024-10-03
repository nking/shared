package algorithms.scheduling;

import algorithms.misc.MiscMath0;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

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
     */
    public void testweightedIntervalBottomUp() {

        //TODO: add more tests and include changing the order of items in this one

        //                              *               *
        double[] s = new double[]{0.5, 1, 2.75, 1.25, 5.25, 6};
        double[] f = new double[]{2.5, 3, 5, 7, 7.5, 8.5};
        double[] v = new double[]{2, 6, 3.5, 7, 8, 1.1};

        Misc misc = new Misc();
        int[] indexes = misc.weightedIntervalBottomUp2(s, f, v);
        //System.out.println("scheduled intervals = " + Arrays.toString(indexes));
        int i;
        double sum = 0;
        for (i = 0; i < indexes.length; ++i) {
            //System.out.printf("%d [%.2f : %.2f]  sum=%.2f\n", indexes[i], s[indexes[i]], f[indexes[i]], sum);
            sum += v[indexes[i]];
        }
        //System.out.println("sum of values=" + sum);
        assertEquals(14.0, sum);

        int[] expected = new int[]{1,4};
        assertTrue(Arrays.equals(expected, indexes));


        // change the order  
        s = new double[]{1.25, 0.5, 5.25, 1, 2.75, 6};
        f = new double[]{7, 2.5, 7.5, 3, 5, 8.5};
        v = new double[]{7, 2, 8, 6, 3.5, 1.1};

        misc = new Misc();
        indexes = misc.weightedIntervalBottomUp2(s, f, v);
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
        s = new double[]{1.25, 0.5, 5.25, 1, 2.75, 6, 1};
        f = new double[]{7, 2.5, 7.5, 3, 5, 8.5, 3};
        v = new double[]{7, 2, 8, 6, 3.5, 1.1, 1};

        misc = new Misc();
        indexes = misc.weightedIntervalBottomUp2(s, f, v);
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

        // --- test the other algorithm for same scheduling:
        s = new double[]{1.25, 0.5, 5.25, 1, 2.75, 6};
        f = new double[]{7, 2.5, 7.5, 3, 5, 8.5};
        v = new double[]{7, 2, 8, 6, 3.5, 1.1};

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

        /*
              0  1  2  3  4  5  6  7  8  9
           0  ---------| *
           1     ------------|
           2           --------| *
           3        ---------------|
           4                   ------|*
           5                   ---------|
         */
        s = new double[]{0, 1, 3, 2,   6, 6};
        f = new double[]{3, 5, 6, 7.5, 8, 9};
        v = new double[]{2, 2, 2, 2,   5, 2};

        expected = new int[]{0, 2, 4};
        indexes = misc.weightedIntervalBottomUp2(s, f, v);
        sum = 0;
        for (i = 0; i < indexes.length; ++i) {
            //System.out.printf("%d [%.2f : %.2f]  sum=%.2f\n", indexes[i], s[indexes[i]], f[indexes[i]], sum);
            sum += v[indexes[i]];
        }
        //System.out.println("scheduled intervals = " + Arrays.toString(indexes));
        assertEquals(9.0, sum);
        assertTrue(Arrays.equals(expected, indexes));
    }

    public void testWeightedGreedy() {
        System.out.println("unit task with weighted deadline scheduling");

        // testing for task scheduling with penalty for lateness

        // from Cormen, Leiserson, Rivest, and Stein Fig 16.7
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
        int[] d = new int[]{4, 2, 4, 3, 1, 4, 6};
        int[] w = new int[]{70, 60, 50, 40, 30, 20, 10};

        Misc misc = new Misc();
        int[] indexes = misc.weightedGreedySingleResource(d, w);

        int[] expected = new int[]{1, 3, 0, 2};//, 6, 5, 4}; // last 3 in any order
        Set<Integer> expectedLate = new HashSet<>();
        expectedLate.add(6); expectedLate.add(5); expectedLate.add(4);
        assertEquals(expected.length + expectedLate.size(), indexes.length);

        for (int i = 0; i < 4; ++i) {
            assertEquals(expected[i], indexes[i]);
        }
        for (int i = 4; i < indexes.length; ++i) {
            assertTrue(expectedLate.remove(indexes[i]));
        }
        assertTrue(expectedLate.isEmpty());
    }

    public void testUnweightedIntervalNoConflicts() {

        double[] s = new double[]{1, 4, 0, 7, 2, 8, 11, 3};
        double[] t = new double[]{5, 6, 9, 10, 11, 13, 14, 15};
        int[] expected = new int[]{1 - 1, 4 - 1, 7 - 1};
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
        double[] s = new double[]{0, 1, 2, 6.5, 7.5, 10, 13, 19, 23, 26, 26, 31};
        double[] f = new double[]{4.5, 10, 6.5, 11, 22, 24, 17.5, 25, 32, 29.5, 33, 34};

        Misc misc = new Misc();
        int[] resources;
        for (int method = 0; method < 2; ++method) {
            if (method == 0) {
                resources = misc.intervalPartitionGreedy2(s, f);
            } else {
                resources = misc.intervalPartitionGreedy(s, f);
            }
            //System.out.println("resources=" + Arrays.toString(resources));
            int cMax = MiscMath0.findMax(resources);
            assertEquals(2, cMax);
            double sPrev, fPrev;
            for (int c = 0; c < cMax; ++c) {
                fPrev = Double.NEGATIVE_INFINITY;
                sPrev = Double.NEGATIVE_INFINITY;
                for (int i = 0; i < resources.length; ++i) {
                    if (resources[i] == c) {
                        assertTrue(s[i] > sPrev);
                        assertTrue(s[i] >= fPrev);
                        assertTrue(f[i] > fPrev);
                        //System.out.printf("color=%d) %.2f -- %.2f\n", c, s[i], f[i]);
                        fPrev = f[i];
                        sPrev = s[i];
                    }
                }
            }
        }

    }

    public void testWeightedOptimal() {
        System.out.println("testWeightedOptimal");

        int[] duration;
        double[] deadline;
        double[] v;
        int[] outputSchedule;
        int[] outLastOnTimeIdx;
        double p;
        double pExp;
        int[] sExp;
        int n;

        /*
            t   d   p
            2   3   2
            2   2   3 <==
            1   3   4 <==
            5   5   5
         */
        duration = new int[]{2, 2, 1, 5};
        n = duration.length;
        deadline = new double[]{3, 2, 3, 5};
        v = new double[]{2, 3, 4, 5};
        outputSchedule = new int[n];
        outLastOnTimeIdx = new int[1];
        pExp = 7.;
        sExp = new int[]{1, 2};

        p = Misc.weightedDynamicSingleResource(duration, deadline, v, outputSchedule, outLastOnTimeIdx);

        System.out.printf("dynamic: p=%.3f, sched=%s\n", p, Arrays.toString(outputSchedule));

        assertEquals(1, outLastOnTimeIdx[0]);
        int i;
        for (i = 0; i <= outLastOnTimeIdx[0]; ++i) {
            assertEquals(sExp[i], outputSchedule[i]);
        }
        assertTrue(Math.abs(pExp - p) < 0.001);
    }

    public void testIntervalPartitionGreedySingleResource() {

        double[] s = new double[]{1, 3, 0, 5, 3, 5, 6, 7, 8, 2, 12};
        double[] f = new double[]{4, 5, 6, 7, 9, 9, 10, 11, 12, 14, 16};
        boolean isSortedByF = false;

        int[] expected = new int[]{1 - 1, 4 - 1, 8 - 1, 11 - 1};

        int[] schedule = Misc.intervalPartitionGreedySingleResource(s, f, isSortedByF);

        assertTrue(Arrays.equals(expected, schedule));
        //System.out.println(Arrays.toString(schedule));
    }

    public void test100() {
        int n = 100;
        int[] duration = new int[n];
        Arrays.fill(duration, 2);

        double[] p = new double[n];
        Arrays.fill(p, 10);

        int i;
        double[] deadline = new double[n];
        deadline[0] = duration[0];
        for (i = 1; i < n; ++i) {
            deadline[i] = duration[i] + deadline[i - 1];
        }
        int[] outputSchedule = new int[n];
        int[] outLastOnTimeIdx = new int[1];

        double sumP = Misc.weightedDynamicSingleResource(
                duration, deadline, p, outputSchedule, outLastOnTimeIdx);
        System.out.printf("dynamic: p=%.3f, sched=%s\n", sumP, Arrays.toString(outputSchedule));
    }
}