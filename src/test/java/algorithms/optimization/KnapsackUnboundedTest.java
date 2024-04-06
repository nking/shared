package algorithms.optimization;

import junit.framework.TestCase;
import java.util.Arrays;

public class KnapsackUnboundedTest extends TestCase {

    public void testMaxValueForTarget() {
        int[] w, v;
        int target, ans, expAns;

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        target = 8;
        expAns = 0;
        System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForTarget(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{1, 5, 10, 30};
        v = new int[]{2, 5, 10, 30};
        target = 10;
        expAns = 20;
        System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForTarget(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{1, 5, 10, 30};
        v = new int[]{2, 11, 10, 30};
        target = 10;
        expAns = 22;
        System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForTarget(v, w, target);
        assertEquals(expAns, ans);
    }

    public void testMaxValueForCapacity() {

        int[] w, v;
        int target, ans, expAns;

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        target = 8;
        expAns = 11;
        //System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForCapacity(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        target = 4;
        expAns = 0;
        //System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForCapacity(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{7, 10, 30};
        target = 30;
        expAns = 42;
        //System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForCapacity(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{5, 10, 30};
        target = 30;
        expAns = 30;
        //System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForCapacity(v, w, target);
        assertEquals(expAns, ans);
    }

    public void testMinNumberOfItemsForTarget() {

        int[] w;
        int target, ans, expAns;

        w = new int[]{1, 5, 10};
        target = 8;
        expAns = 4;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{1, 5};
        target = 8;
        expAns = 4;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{2, 5, 10};
        target = 8;
        expAns = 4;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{1, 5, 7, 10};
        target = 12;
        expAns = 2;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{2, 5, 10};
        target = 1;
        expAns = 0;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{1, 5, 7, 10};
        target = 7;
        expAns = 1;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForTarget(w, target);
        assertEquals(expAns, ans);
    }

    public void testCountNumberOfWaysForTarget() {

        int[] w;
        int target, ans, expAns;

        // 4 7s
        // 14 2's
        // 2 7's, 7 2's
        w = new int[]{2,7};
        target = 28;
        expAns = 3;
        ans = KnapsackUnbounded.numberOfWaysForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{2,7};
        target = 7;
        expAns = 1;
        ans = KnapsackUnbounded.numberOfWaysForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{2,7};
        target = 4;
        expAns = 1;
        ans = KnapsackUnbounded.numberOfWaysForTarget(w, target);
        assertEquals(expAns, ans);
    }

}
