package algorithms.optimization;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class Knapsack01Test extends TestCase {
    
    public Knapsack01Test(String testName) {
        super(testName);
    }

    public void testMaxValueForTarget() {
        int[] w, v;
        int target, ans, expAns;

        w = new int[]{5, 10, 30};
        v = new int[]{11, 23, 30};
        target = 15;
        expAns = 34;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.maxValueForTarget(v, w, target);
        assertEquals(expAns, ans);
        ans = Knapsack01.maxValueForTarget2(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 23, 30};
        target = 16;
        expAns = 0;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.maxValueForTarget(v, w, target);
        assertEquals(expAns, ans);
        ans = Knapsack01.maxValueForTarget2(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30, 40};
        v = new int[]{11, 23, 30, 40};
        target = 45;
        expAns = 64;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.maxValueForTarget(v, w, target);
        assertEquals(expAns, ans);
        ans = Knapsack01.maxValueForTarget2(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30, 40};
        v = new int[]{11, 23, 30, 40};
        target = 10;
        expAns = 23;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.maxValueForTarget(v, w, target);
        assertEquals(expAns, ans);
        ans = Knapsack01.maxValueForTarget2(v, w, target);
        assertEquals(expAns, ans);
    }

    public void testMaxValueForCapacity() {
        int[] w, v;
        int capacity, ans, expAns;

        w = new int[]{5, 10, 30};
        v = new int[]{11, 23, 30};
        capacity = 16;
        expAns = 34;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.maxValueForCapacity(v, w, capacity);
        assertEquals(expAns, ans);
        ans = Knapsack01.maxValueForCapacity2(v, w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 23, 30};
        capacity = 15;
        expAns = 34;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.maxValueForCapacity(v, w, capacity);
        assertEquals(expAns, ans);
        ans = Knapsack01.maxValueForCapacity2(v, w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 23, 30};
        capacity = 4;
        expAns = 0;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.maxValueForCapacity(v, w, capacity);
        assertEquals(expAns, ans);
        ans = Knapsack01.maxValueForCapacity2(v, w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 23, 30};
        capacity = 8;
        expAns = 11;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.maxValueForCapacity(v, w, capacity);
        assertEquals(expAns, ans);
        ans = Knapsack01.maxValueForCapacity2(v, w, capacity);
        assertEquals(expAns, ans);
    }

    public void testNumberOfWaysForTarget() {
        int[] w;
        int target, ans, expAns;

        w = new int[]{5, 10, 30};
        target = 15;
        expAns = 1;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.numberOfWaysForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        target = 16;
        expAns = 0;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.numberOfWaysForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30, 40};
        target = 45;
        expAns = 2;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.numberOfWaysForTarget(w, target);
        assertEquals(expAns, ans);
    }

    public void testMinNumberOfItemsForTarget() {
        int[] w;
        int target, ans, expAns;

        w = new int[]{5, 10, 30};
        target = 15;
        expAns = 2;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.minNumberOfItemsForTarget(w, target);
        assertEquals(expAns, ans);
        ans = Knapsack01.minNumberOfItemsForTarget2(w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        target = 16;
        expAns = 0;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.minNumberOfItemsForTarget(w, target);
        assertEquals(expAns, ans);
        ans = Knapsack01.minNumberOfItemsForTarget2(w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30, 40};
        target = 45;
        expAns = 2;
        //System.out.printf("\target=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.minNumberOfItemsForTarget(w, target);
        assertEquals(expAns, ans);
        ans = Knapsack01.minNumberOfItemsForTarget2(w, target);
        assertEquals(expAns, ans);
    }

    public void testMinNumberOfItemsForCapacity() {
        int[] w;
        int capacity, ans, expAns;

        w = new int[]{5, 10, 30};
        capacity = 15;
        expAns = 2;
        //System.out.printf("\capacity=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = Knapsack01.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        capacity = 16;
        expAns = 2;
        //System.out.printf("\capacity=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = Knapsack01.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30, 40};
        capacity = 45;
        expAns = 2;
        //System.out.printf("\capacity=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = Knapsack01.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30, 40};
        capacity = 4;
        expAns = 0;
        //System.out.printf("\capacity=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = Knapsack01.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = Knapsack01.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);
    }
    
    public void test011() {

        // values
        int[] ratingsOfObjectOfInterest = new int[]{30, 14, 16, 9};

        // weights
        int[] distances = new int[]{6,  3,  4, 2};

        int capacity = 10;

        int expSumValue = 46;

        int result = Knapsack01.maxValueForCapacity(ratingsOfObjectOfInterest, distances, capacity);

        assertTrue(result == expSumValue);
        
        result = Knapsack01.approxDynamically(ratingsOfObjectOfInterest, distances, capacity);

        assertTrue(result == expSumValue);
        
        /*
        memo[5][11]
              0     1     2     3     4     5     6     7     8     9     10
        0     0     0     0     0     0     0     0     0     0     0     0
        1     0     0     0     0     0     0    30    30    30    30    30
        2     0     0     0    14    14    14    30    30    30    44    44
        3     0     0     0    14    16    16    30    30    30    44    46
        4     0     0     9    14    16    23    30    30    39    44    46
        
        memo[5][6] 
            0     1     2     3     4     5
        0   0     0     0     0     0     0
        1   0     0     0    30    30    30
        2   0     0    14    30    30    44
        3   0     0    16    30    30    46
        4   0     9    16    30    39    46
        */
    }
    
    public void test012() {

        int[] weights = new int[] {
            9, 13, 153, 50, 15, 68, 27, 39, 23, 52,
            11,32,24,48,73,42,43,22,7,18,4,30 
        };
        int[] values = new int[] {
            150,35,200,160,60,45,60,40,30,10,70,30,
            15,10,40,70,75,80,20,12,50,10
        };

        int result;
        int expSumValue = 1030;

        result = Knapsack01.maxValueForCapacity(values, weights, 400);
        
        //System.out.println("result: " +result);

        assertTrue(result == expSumValue);
        
        result = Knapsack01.approxDynamically(values, weights, 400);

       // System.out.println("result approx: " +result);
                
    }
        
}
