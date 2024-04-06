package algorithms.optimization;

import junit.framework.TestCase;

import java.util.Arrays;

public class KnapsackBoundedTest extends TestCase {

    public void testMaxValueForCapacity() {

        int[] w, v, q;
        int capacity, ans, expAns;

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{1, 2, 2};
        capacity = 8;
        expAns = 11;
        System.out.printf("\ncapacity=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
                capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{0,2,2};
        capacity = 8;
        expAns = 0;
        System.out.printf("\ncapacity=%d\n  values=%s\n  weights=%s\n  quantities=%s\n",
                capacity, Arrays.toString(v), Arrays.toString(w), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{5, 11, 30};
        q = new int[]{0,3,2};
        capacity = 32;
        expAns = 33;
        System.out.printf("\ncapacity=%d\n  values=%s\n  weights=%s\n  quantities=%s\n",
                capacity, Arrays.toString(v), Arrays.toString(w), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{5, 11, 30};
        q = new int[]{0,2,2};
        capacity = 32;
        expAns = 30;
        System.out.printf("\ncapacity=%d\n  values=%s\n  weights=%s\n  quantities=%s\n",
                capacity, Arrays.toString(v), Arrays.toString(w), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{2, 2, 2};
        capacity = 10;
        expAns = 22;
        System.out.printf("\ncapacity=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
                capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{1, 2, 2};
        capacity = 10;
        expAns = 11;
        System.out.printf("\ncapacity=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
                capacity, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForCapacity(v, w, q, capacity);
        assertEquals(expAns, ans);
    }

    public void testMaxValueForTarget() {

        int[] w, v, q;
        int target, ans, expAns;

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{1, 2, 2};
        target = 10;
        expAns = 10; // if using maxValueForCapacity, this would be 11
        System.out.printf("\ntarget=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
                target, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForTarget(v, w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{1, 2, 2};
        target = 25;
        expAns = 31;
        System.out.printf("\ntarget=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
                target, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForTarget(v, w, q, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        q = new int[]{3, 2, 2};
        target = 25;
        expAns = 43;//3 5s and 1 10 = 33+10
        System.out.printf("\ntarget=%d\n  weights=%s\n values=%s\n   quantities=%s\n",
                target, Arrays.toString(w), Arrays.toString(v), Arrays.toString(q));
        ans = KnapsackBounded.maxValueForTarget(v, w, q, target);
        assertEquals(expAns, ans);
    }

    public void testNumberOfWaysForTarget() {

    }

    public void testMinNumberOfItemsForTarget() {

    }

}
