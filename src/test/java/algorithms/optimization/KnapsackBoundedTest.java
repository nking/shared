package algorithms.optimization;

import junit.framework.TestCase;

import java.util.Arrays;

public class KnapsackBoundedTest extends TestCase {

    public void testMaxValueForCapacity() {

        int[] w, v, q;
        int capacity, ans, expAns;

        w = new int[]{5, 10, 30};
        v = new int[]{11, 23, 30};
        q = new int[]{2, 2, 2};
        capacity = 16;
        expAns = 34;  // 1 5 and 1 10 = 11+23=34
        //System.out.printf("\ncapacity=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForCapacity2(v, w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{1, 2, 2};
        capacity = 8;
        expAns = 11;// 5
        //System.out.printf("\ncapacity=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForCapacity2(v, w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{0,2,2};
        capacity = 8;
        expAns = 0;
        //System.out.printf("\ncapacity=%d\n  values=%s\n  weights=%s\n  quantities=%s\n",
        //        capacity, Arrays.toString(v), Arrays.toString(w), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForCapacity2(v, w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{5, 11, 30};
        q = new int[]{0,3,2};
        capacity = 32;
        expAns = 33; // 10,10,10
        //System.out.printf("\ncapacity=%d\n  values=%s\n  weights=%s\n  quantities=%s\n",
        //        capacity, Arrays.toString(v), Arrays.toString(w), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForCapacity2(v, w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{5, 11, 30};
        q = new int[]{0,2,2};
        capacity = 32;
        expAns = 30;//30
        //System.out.printf("\ncapacity=%d\n  values=%s\n  weights=%s\n  quantities=%s\n",
        //        capacity, Arrays.toString(v), Arrays.toString(w), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForCapacity2(v, w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{2, 2, 2};
        capacity = 10;
        expAns = 22;//5,5
        //System.out.printf("\ncapacity=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForCapacity2(v, w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{1, 2, 2};
        capacity = 10;
        expAns = 11; //5
        //System.out.printf("\ncapacity=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForCapacity2(v, w, q, capacity);
        assertEquals(expAns, ans);
    }

    public void testMaxValueForTarget() {

        int[] w, v, q;
        int target, ans, expAns;

        w = new int[]{5, 10, 30};
        v = new int[]{5, 10, 30};
        q = new int[]{1, 2, 2};
        target = 35;
        expAns = 35;//5,30
        //System.out.printf("\ntarget=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        target, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForTarget(v, w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForTarget2(v, w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{1, 2, 2};
        target = 10;
        expAns = 10; //10
        //System.out.printf("\ntarget=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        target, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForTarget(v, w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForTarget2(v, w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{2, 2, 2};
        target = 10;
        expAns = 22; //2 5's
        //System.out.printf("\ntarget=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        target, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForTarget(v, w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForTarget2(v, w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{1, 2, 2};
        target = 25;
        expAns = 31;// 1 5 avail + 2 10s = 11 + 20 = 31
        //System.out.printf("\ntarget=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        target, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForTarget(v, w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForTarget2(v, w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{3, 2, 2};
        target = 25;
        expAns = 43;//3 5s and 1 10 = 33+10
        //System.out.printf("\ntarget=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        target, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForTarget(v, w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForTarget2(v, w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{3, 2, 2};
        target = 1;
        expAns = 0;
        //System.out.printf("\ntarget=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        target, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForTarget(v, w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForTarget2(v, w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{10, 30};
        v = new int[]{10, 30};
        q = new int[]{1, 1};
        target = 50;
        expAns = 0;
        //System.out.printf("\ntarget=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
        //        target, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForTarget(v, w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.maxValueForTarget2(v, w, q, target);
        assertEquals(expAns, ans);
    }

    public void testNumberOfWaysForTarget() {
        int[] w, q;
        int target, ans, expAns;

        w = new int[]{1, 2, 5, 10};
        q = new int[]{11, 5, 2, 1};
        target = 11;
        expAns = 12;//[10,1], [5,5,1],[5,2,2,2],[5,2,2,1,1],[5,2,1,1,1,1],[5,1,1,1,1,1,1],
        //[2,2,2,2,2,1], [2,2,2,2,1,1,1], [2,2,2,1,1,1,1,1], [2,2,1,1,1,1,1,1,1]
        // [2,1,1,1,1,1,1,1,1,1], [11 1's]
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackBounded.numberOfWaysForTarget(w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{1, 2, 5, 10};
        q = new int[]{3, 5, 1, 1};
        target = 11;
        expAns = 5;//[10,1], [5,2,2,2],[5,2,2,1,1],
                //[2,2,2,2,2,1], [2,2,2,2,1,1,1],
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackBounded.numberOfWaysForTarget(w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{1, 2, 5, 10};
        q = new int[]{3, 2, 1, 1};
        target = 11;
        expAns = 2;//[10,1], [5,2,2,1,1],
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackBounded.numberOfWaysForTarget(w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{2, 5, 10};
        q = new int[]{2, 1, 1};
        target = 11;
        expAns = 0;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackBounded.numberOfWaysForTarget(w, q, target);
        assertEquals(expAns, ans);
    }

    public void testMinNumberOfItemsForTarget() {

        int[] w, q;
        int target, ans, expAns;

        w = new int[]{1, 2, 5, 10};
        q = new int[]{3, 2, 1, 1};
        target = 11;
        expAns = 2;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackBounded.minNumberOfItemsForTarget(w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.minNumberOfItemsForTarget2(w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{1, 2, 5, 10};
        q = new int[]{3, 2, 1, 1};
        target = 9;
        expAns = 3; // 2, 2, 5
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackBounded.minNumberOfItemsForTarget(w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.minNumberOfItemsForTarget2(w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{1, 2, 5, 10};
        q = new int[]{3, 1, 1, 1};
        target = 9;
        expAns = 4;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackBounded.minNumberOfItemsForTarget(w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.minNumberOfItemsForTarget2(w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10};
        q = new int[]{1, 1};
        target = 4;
        expAns = 0;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackBounded.minNumberOfItemsForTarget(w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.minNumberOfItemsForTarget2(w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10};
        q = new int[]{1, 1};
        target = 20;
        expAns = 0;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackBounded.minNumberOfItemsForTarget(w, q, target);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.minNumberOfItemsForTarget2(w, q, target);
        assertEquals(expAns, ans);
    }

    public void testMinNumberOfItemsForCapacity() {
        //System.out.printf("testMinNumberOfItemsForCapacity");

        int[] w, q;
        int capacity, ans, expAns;

        w = new int[]{1, 2, 5, 10};
        q = new int[]{3, 2, 1, 1};
        capacity = 11;
        expAns = 2; // 1, 10
        //System.out.printf("\ncapacity=%d\nweights=%s\nq=%s\n", capacity, Arrays.toString(w), Arrays.toString(q));
        ans = KnapsackBounded.minNumberOfItemsForCapacity(w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.minNumberOfItemsForCapacity2(w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{1, 2, 5, 10};
        q = new int[]{3, 1, 1, 1};
        capacity = 9;
        expAns = 4; // 1,1,2,5
        //System.out.printf("\ncapacity=%d\nweights=%s\nq=%s\n", capacity, Arrays.toString(w), Arrays.toString(q));
        ans = KnapsackBounded.minNumberOfItemsForCapacity(w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.minNumberOfItemsForCapacity2(w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10};
        q = new int[]{1, 1};
        capacity = 4;
        expAns = 0;
        //System.out.printf("\ncapacity=%d\nweights=%s\nq=%s\n", capacity, Arrays.toString(w), Arrays.toString(q));
        ans = KnapsackBounded.minNumberOfItemsForCapacity(w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.minNumberOfItemsForCapacity2(w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10};
        q = new int[]{1, 1};
        capacity = 10;
        expAns = 1;//10
        //System.out.printf("\ncapacity=%d\nweights=%s\nq=%s\n", capacity, Arrays.toString(w), Arrays.toString(q));
        ans = KnapsackBounded.minNumberOfItemsForCapacity(w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.minNumberOfItemsForCapacity2(w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10};
        q = new int[]{1, 2};
        capacity = 30;
        expAns = 3;//5,10,10
        //System.out.printf("\ncapacity=%d\nweights=%s\nq=%s\n", capacity, Arrays.toString(w), Arrays.toString(q));
        ans = KnapsackBounded.minNumberOfItemsForCapacity(w, q, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackBounded.minNumberOfItemsForCapacity2(w, q, capacity);
        assertEquals(expAns, ans);
    }

}
