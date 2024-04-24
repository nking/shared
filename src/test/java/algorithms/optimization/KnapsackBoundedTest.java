package algorithms.optimization;

import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
        //[10,1],   => math.factorial(2)
        // [5,5,1],  => math.factorial(3)/(math.factorial(2))
        // [5,2,2,2], => math.factorial(4)/(math.factorial(3))
        // [5,2,2,1,1], => math.factorial(5)/(math.factorial(2)*math.factorial(2))
        // [5,2,1,1,1,1], => math.factorial(6)/(math.factorial(4))
        // [5,1,1,1,1,1,1], => math.factorial(7)/(math.factorial(6))
        //[2,2,2,2,2,1], => math.factorial(6)/(math.factorial(5))
        // [2,2,2,2,1,1,1], => math.factorial(7)/(math.factorial(3)*math.factorial(4))
        // [2,2,2,1,1,1,1,1], => math.factorial(8)/(math.factorial(3)*math.factorial(5))
        // [2,2,1,1,1,1,1,1,1] => math.factorial(9)/(math.factorial(2)*math.factorial(7))
        // [2,1,1,1,1,1,1,1,1,1], => math.factorial(10)/(math.factorial(9))
        // [11 1's] => 1
        // math.factorial(2) + (math.factorial(3)/(math.factorial(2))) + (math.factorial(4)/(math.factorial(3))) + (math.factorial(5)/(math.factorial(2)*math.factorial(2))) + (math.factorial(6)/(math.factorial(4))) + (math.factorial(7)/(math.factorial(6))) + (math.factorial(6)/(math.factorial(5))) + (math.factorial(7)/(math.factorial(3)*math.factorial(4))) + (math.factorial(8)/(math.factorial(3)*math.factorial(5))) + (math.factorial(9)/(math.factorial(2)*math.factorial(7))) + (math.factorial(10)/(math.factorial(9))) + 1
        //expAns = 220;
        //ans = KnapsackBounded.numberOfSequencesForTarget(w, q, target);
        //assertEquals(expAns, ans);

        w = new int[]{1, 2, 5, 10};
        q = new int[]{3, 5, 1, 1};
        target = 11;
        expAns = 5;//[10,1], [5,2,2,2],[5,2,2,1,1],
                //[2,2,2,2,2,1], [2,2,2,2,1,1,1],
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackBounded.numberOfWaysForTarget(w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{1,2};
        q = new int[]{1,1};
        target = 3;
        expAns = 1;//[1,2]
        ans = KnapsackBounded.numberOfWaysForTarget(w, q, target);
        assertEquals(expAns, ans);
        //System.out.printf("\nw=%s\nq=%s\ntarget=%d\n", Arrays.toString(w), Arrays.toString(q), target);
        // [1,2] =>2
        long expAns2 = 2;
        long ans2 = KnapsackBounded.numberOfSequencesForTarget(w, q, target);
        assertEquals(expAns2, ans2);

        w = new int[]{1,2};
        q = new int[]{3,1};
        target = 3;
        expAns = 2;//[1,2]
        ans = KnapsackBounded.numberOfWaysForTarget(w, q, target);
        assertEquals(expAns, ans);
        // [1,2] => 12,21
        // [1,1,1] => 1
        //System.out.printf("\nw=%s\nq=%s\ntarget=%d\n", Arrays.toString(w), Arrays.toString(q), target);
        expAns2 = 3;
        System.out.printf("q=3,1\n");
        ans2 = KnapsackBounded.numberOfSequencesForTarget(w, q, target);
        assertEquals(expAns2, ans2);

        w = new int[]{1, 2, 5, 10};
        q = new int[]{3, 2, 1, 1};
        target = 11;
        expAns = 2;//[10,1], [5,2,2,1,1],
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackBounded.numberOfWaysForTarget(w, q, target);
        assertEquals(expAns, ans);
        // 10,1 => 2
        // [5,2,2,1,1] => math.factorial(5)/(math.factorial(2)*math.factorial(2)) = 30
        //System.out.printf("\nw=%s\nq=%s\ntarget=%d\n", Arrays.toString(w), Arrays.toString(q), target);
        expAns2 = 32;
        ans2 = KnapsackBounded.numberOfSequencesForTarget(w, q, target);
        assertEquals(expAns2, ans2);

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
