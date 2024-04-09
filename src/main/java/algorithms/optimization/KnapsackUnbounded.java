package algorithms.optimization;

import java.util.Arrays;

/**
 * various unbounded knapsack problems
 */
public class KnapsackUnbounded {

    /**
     * find maximum value from summing the values associated with an unbounded quantity
     * of weights that sum to exactly EQ target.
     * @param values non-negative values for items. values and weights describe the same items.
     * @param weights non-negative array of weights for items
     * @param target the exact sum that a combination of an unbounded number of weights should sum to
     * @return the maximum value for items whose weights sum up to exactly equal target.
     */
    public static int maxValueForTarget(int[] values, int[] weights, int target) {
        return maxValue(values, weights, target, true);
    }

    /**
     * find maximum value from summing the values associated with an unbounded quantity of
     * weights that sum to LEQ to capacity.
     * @param values non-negative values for items. values and weights describe the same items.
     * @param weights non-negative array of weights for items
     * @param capacity upper limit to the sum of weights for an unbounded number of items
     * @return
     */
    public static int maxValueForCapacity(int[] values, int[] weights, int capacity) {
        return maxValue(values, weights, capacity, false);
    }

    public static int maxValue(int[] values, int[] weights, int target, boolean solveForExact) {
        int n = weights.length;

        if (values.length != n) {
            throw new IllegalStateException("values and weights must be same length");
        }

        // tab[wc] holds the largest sum of values for the items whose weights sum to wc
        int[] tab = new int[target + 1];
        int sentinel = Integer.MIN_VALUE;;
        Arrays.fill(tab, sentinel);
        tab[0] = 0;

        int i, wc, wc2;
        for (i = 0; i < n; i++) {
            if (weights[i] > target) continue;
            // traverse wc from 1 to target to be able to add current wc2 counts to current count sum
            for (wc = 1; wc <= target; wc++) {
                wc2 = wc - weights[i];

                if (wc2 == 0) {
                    tab[wc] = Math.max(tab[wc], values[i]);
                } else if (wc2 > 0)  {
                    int s = (tab[wc2] == sentinel) ? sentinel : tab[wc2] +  values[i];
                    tab[wc] = Math.max(tab[wc], s);
                }
            }
        }

        if (solveForExact) {
            return tab[target] == sentinel ? 0 : tab[target];
        }

        // search backwards for max
        int max = tab[target];
        i = target;
        while (i > 0) {
            max = Math.max(max, tab[--i]);
        }

        return max == sentinel ? 0 : max;
    }

    /**
     * find maximum value from summing the values associated with an unbounded quantity of
     * weights that sum to LEQ to capacity.
     * @param values non-negative values for items. values and weights describe the same items.
     * @param weights non-negative array of weights for items
     * @param capacity upper limit to the sum of weights for an unbounded number of items
     * @return
     */
    public static int maxValueForCapacity2(int[] values, int[] weights, int capacity) {

        int n = weights.length;

        if (values.length != n) {
            throw new IllegalStateException("values and weights must be same length");
        }

        // tab[wc] holds the largest sum of values for the items whose weights sum to wc
        int[] tab = new int[capacity + 1];

        int i, wc, wc2;
        for (i = 0; i < n; i++) {
            if (weights[i] > capacity) continue;
            // cannot traverse weight sum in reverse order here.
            // traverse wc from 1 to target to be able to add current wc2 counts to current count sum
            for (wc = 1; wc <= capacity; wc++) {
                wc2 = wc - weights[i];

                //System.out.printf("i=%d, wc=%d, wc2=%d, weights[i]=%d  tab=%s\n", i, wc, wc2, weights[i], Arrays.toString(tab));

                if (wc2 >= 0) {
                    tab[wc] = Math.max(tab[wc], tab[wc2] + values[i]);
                }
                //System.out.printf("    tab[%d]=%d\n", wc, tab[wc]);
            }
        }
        return tab[capacity];
    }

    /**
     * find the minimum number of items whose weights sum to exactly EQ target.
     * @param weights non-negative array of item weights
     * @param target the exact sum that a combination of and unbounded quantity of weights should sum to
     * @return
     */
    public static int minNumberOfItemsForTarget(int[] weights, int target) {
        return minNumberOfItems(weights, target, true);
    }

    /**
     * find the minimum number of items whose weights sum to LEQ capacity.
     * @param weights non-negative array of item weights
     * @param capacity the sum that a combination of and unbounded quantity of weights should sum to if possible, else
     *                 the min for the greatest lower sum which was met.
     * @return
     */
    public static int minNumberOfItemsForCapacity(int[] weights, int capacity) {
        return minNumberOfItems(weights, capacity, false);
    }

    public static int minNumberOfItems(int[] weights, int target, boolean solveForExact) {
        int n = weights.length;

        // tab[wc] holds the min number of items for the items whose weights sum to wc
        int[] tab = new int[target+1];
        int sentinel = Integer.MAX_VALUE;//target+1;
        Arrays.fill(tab, sentinel);
        tab[0] = 0;

        int i, wc, wc2;
        for (i = 0; i < n; i++) {
            if (weights[i] > target) continue;
            for (wc = 1; wc <= target; wc++) {
                wc2 = wc - weights[i];
                if (wc2 == 0) {
                    tab[wc] = Math.min(tab[wc], 1);
                } else if (wc2 > 0)  {
                    int s = (tab[wc2] == sentinel) ? sentinel : 1 + tab[wc2];
                    tab[wc] = Math.min(tab[wc], s);
                }

                //System.out.printf("i=%d, wc=%d, wc2=%d, weights[i]=%d  tab=%s\n", i, wc, wc2, weights[i], Arrays.toString(tab));
            }
        }
        //System.out.printf("tab=%s\n", Arrays.toString(tab));
        if (solveForExact) {
            return tab[target] == sentinel ? 0 : tab[target];
        }

        // search backwards for last sum entered
        i = target;
        while (i > 0 && tab[i] == sentinel) {
            --i;
        }

        return tab[i] == sentinel ? 0 : tab[i];
    }

    /**
     * find the minimum number of items whose weights sum to LEQ capacity.
     * @param weights non-negative array of item weights
     * @param capacity the sum that a combination of and unbounded quantity of weights should sum to if possible, else
     *                 the min for the greatest lower sum which was met.
     * @return
     */
    public static int minNumberOfItemsForCapacity2(int[] weights, int capacity) {

        int n = weights.length;

        // tab[wc] holds the min number of items for the items whose weights sum to wc
        int[] tab = new int[capacity+1];
        int sentinel = Integer.MAX_VALUE;
        Arrays.fill(tab, sentinel);
        tab[0] = 0;

        int i, wc, wc2;
        for (i = 0; i < n; i++) {
            if (weights[i] > capacity) continue;
            for (wc = 1; wc <= capacity; wc++) {
                wc2 = wc - weights[i];
                if (wc2 >= 0) {
                    int s = (tab[wc2] == sentinel) ? sentinel : 1 + tab[wc2];
                    tab[wc] = Math.min(tab[wc], s);
                }
            }
        }
        return tab[capacity] == sentinel ? 0 : tab[capacity];
    }

    /**
     * count the number of ways that a combination of an unbounded quantity of weights
     * can sum up to exactly EQ target.
     * The combinations are counted as sets, that is, same numbers in a set are counted only once.
     * @param target the exact sum that a combination of and unbounded quantity of weights should sum to
     * @param weights non-negative array of item weights
     * @return
     */
    public static int numberOfWaysForTarget(int[] weights, int target) {

        int n = weights.length;

        // tab[wc] holds the number of ways that the item weights sum to wc.

        int[] tab = new int[target + 1];
        tab[0] = 1;

        int wc, wc2;
        for (int i = 0; i < n; ++i) {
            if (weights[i] > target) continue;
            // traverse wc from 1 to target to be able to add current wc2 counts to current count sum
            for (wc = 1; wc <= target; ++wc) {
                // the remaining sum after weight subtracted
                wc2 = wc - weights[i];
                if (wc2 >= 0) {
                    // adds counts from remaining sum
                    tab[wc] += tab[wc2];
                }
            }
        }
        //System.out.printf("tab=%s\n", Arrays.toString(tab));
        return tab[target];
    }

    /**
     * count the number of ways that a combination of an unbounded quantity of weights
     * can sum up to exactly EQ target, where the sequences are counted rather than sets,
     * e.g. [1,2] is counted and [2,1] is counted.
     * @param target the exact sum that a combination of and unbounded quantity of weights should sum to
     * @param weights non-negative array of item weights
     * @return
     */
    public static int numberOfSequencesForTarget(int[] weights, int target) {

        int n = weights.length;

        // tab[wc] holds the number of ways that the item weights sum to wc.

        int[] tab = new int[target + 1];
        tab[0] = 1;

        int i, wc, wc2;
        for (wc = 1; wc <= target; ++wc) {
            for (i = 0; i < n; ++i) {
                // the remaining sum after weight subtracted
                wc2 = wc - weights[i];
                if (wc2 >= 0) {
                    // adds counts from remaining sum
                    tab[wc] += tab[wc2];
                }
            }
        }
        //System.out.printf("tab=%s\n", Arrays.toString(tab));
        return tab[target];
    }
}
