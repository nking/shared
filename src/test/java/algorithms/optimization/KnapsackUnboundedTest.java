package algorithms.optimization;

import junit.framework.TestCase;
import java.util.Arrays;

public class KnapsackUnboundedTest extends TestCase {

    public void testMaxValueForTarget() {
        int[] w, v;
        int target, ans, expAns;

        w = new int[]{5, 11, 30};
        v = new int[]{11, 10, 30};
        target = 26;
        expAns = 43;//11 + 15 = 3 5's + 1 10 = 3*11 + 1*10 = 33+10=43
        //System.out.printf("\ntarget=%d\n  weights=%s\n values=%s\n",
        //        target, Arrays.toString(w), Arrays.toString(v));
        ans = KnapsackUnbounded.maxValueForTarget(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        target = 8;
        expAns = 0; //
        //System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForTarget(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{1, 5, 10, 30};
        v = new int[]{2, 5, 10, 30};
        target = 10;
        expAns = 20;
        //System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForTarget(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{1, 5, 10, 30};
        v = new int[]{2, 11, 10, 30};
        target = 10;
        expAns = 22;
        //System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForTarget(v, w, target);
        assertEquals(expAns, ans);
    }

    public void testMaxValueForCapacity() {

        int[] w, v;
        int target, ans, expAns;

        w = new int[]{5, 11, 30};
        v = new int[]{11, 23, 30};
        target = 26;
        expAns = 56;// 2 11's = 46. 5 5's =55;  3 5's and 1 11 = 33+23=56
        //System.out.printf("\ncapacity=%d\n  weights=%s\n values=%s\n",
        //        target, Arrays.toString(w), Arrays.toString(v));
        ans = KnapsackUnbounded.maxValueForCapacity(v, w, target);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.maxValueForCapacity2(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        target = 8;
        expAns = 11;
        //System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForCapacity(v, w, target);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.maxValueForCapacity2(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{11, 10, 30};
        target = 4;
        expAns = 0;
        //System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForCapacity(v, w, target);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.maxValueForCapacity2(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{7, 10, 30};
        target = 30;
        expAns = 42;
        //System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForCapacity(v, w, target);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.maxValueForCapacity2(v, w, target);
        assertEquals(expAns, ans);

        w = new int[]{5, 10, 30};
        v = new int[]{5, 10, 30};
        target = 30;
        expAns = 30;
        //System.out.printf("\ntarget=%d\n  values=%s\n  weights=%s\n", target, Arrays.toString(v), Arrays.toString(w));
        ans = KnapsackUnbounded.maxValueForCapacity(v, w, target);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.maxValueForCapacity2(v, w, target);
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

    public void testCountNumberOfSetsAndSequencesForTarget() {

        int[] w;
        int target, ans, expAns;

        // 4 7s
        // 14 2's
        // 2 7's, 7 2's
        w = new int[]{2,7};
        target = 28;
        expAns = 3;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.numberOfWaysForTarget(w, target);
        assertEquals(expAns, ans);
        // count sequences instead of sets:
        // 4 7s         => 1 seq
        // 14 2's       => 1 seq
        // 2 7's, 7 2's =>  9!/(2!*7!) = 36 sequences
        expAns = 38;
        ans = KnapsackUnbounded.numberOfSequencesForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{2,7};
        target = 7;
        expAns = 1;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.numberOfWaysForTarget(w, target);
        assertEquals(expAns, ans);
        expAns = 1;
        ans = KnapsackUnbounded.numberOfSequencesForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{2,7};
        target = 4;
        expAns = 1;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.numberOfWaysForTarget(w, target);
        assertEquals(expAns, ans);
        expAns = 1;
        ans = KnapsackUnbounded.numberOfSequencesForTarget(w, target);
        assertEquals(expAns, ans);

        w = new int[]{2,7};
        target = 5;
        expAns = 0;
        //System.out.printf("\ntarget=%d, weights=%s\n", target, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.numberOfWaysForTarget(w, target);
        assertEquals(expAns, ans);
        expAns = 0;
        ans = KnapsackUnbounded.numberOfSequencesForTarget(w, target);
        assertEquals(expAns, ans);


    }

    public void testMinNumberOfItemsForCapacity() {

        int[] w;
        int capacity, ans, expAns;

        w = new int[]{2, 5, 10};
        capacity = 12;
        expAns = 2;
        //System.out.printf("\ncapacity=%d, weights=%s\n", capacity, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{1, 5, 10};
        capacity = 8;
        expAns = 4;
        //System.out.printf("\ncapacity=%d, weights=%s\n", capacity, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{1, 5};
        capacity = 8;
        expAns = 4;
        //System.out.printf("\ncapacity=%d, weights=%s\n", capacity, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{1, 5, 7, 10};
        capacity = 12;
        expAns = 2;
        //System.out.printf("\ncapacity=%d, weights=%s\n", capacity, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{2, 5, 10};
        capacity = 1;
        expAns = 0;
        //System.out.printf("\ncapacity=%d, weights=%s\n", capacity, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{1, 5, 7, 10};
        capacity = 7;
        expAns = 1;
        //System.out.printf("\ncapacity=%d, weights=%s\n", capacity, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{1, 5, 7, 10};
        capacity = 8;
        expAns = 2;
        //System.out.printf("\ncapacity=%d, weights=%s\n", capacity, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);

        w = new int[]{1, 5, 7, 10};
        capacity = 28;
        expAns = 4;
        //System.out.printf("\ncapacity=%d, weights=%s\n", capacity, java.util.Arrays.toString(w));
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity(w, capacity);
        assertEquals(expAns, ans);
        ans = KnapsackUnbounded.minNumberOfItemsForCapacity2(w, capacity);
        assertEquals(expAns, ans);
    }

}
