package algorithms.optimization;

import java.util.Arrays;
import java.util.List;

public class KnapsackBounded {

    /**
     * find maximum value from summing the values associated with bounded quantities of
     * weights that sum to LEQ to capacity.
     *
     <pre>runtime complexity is O(n*capacity*max(quantities) where n = weights.length
     space complexity is O(capacity).
     </pre>
     *
     * @param values   values for items. values and weights describe the same items.
     * @param weights  non-negative array of weights for items
     * @param the      non-negative quantities of each item available.
     * @param capacity upper limit to the sum of weights for bounded quantities of items
     * @return
     */
    public static int maxValueForCapacity2(int[] values, int[] weights, int[] quantities, int capacity) {

        int n = weights.length;

        if (values.length != n) {
            throw new IllegalStateException("values and weights must be same length");
        }
        if (quantities.length != n) {
            throw new IllegalStateException("quantities and weights must be same length");
        }

        int[] tab = new int[capacity + 1];

        int wc, q, wc2;
        for (int i = 0; i < n; ++i) {
            // since tab holds current and prev i results,
            // need to traverse weights from high to low
            // to avoid including an updated low wc2 in current wc
            for (wc = capacity; wc >= weights[i]; --wc) {
                for (q = 1; q <= quantities[i]; ++q) {
                    wc2 = wc - q * weights[i];
                    if (wc2 >= 0) {
                        tab[wc] = Math.max(tab[wc], tab[wc2] + q * values[i]);
                    }
                }
            }
        }
        return tab[capacity];
    }

    public static int maxValueForCapacity(int[] values, int[] weights, int[] quantities, int capacity) {

        return maxValue(values, weights, quantities, capacity, false);
    }

    /**
     * find maximum value from summing the values associated with bounded quantities
     * of weights that sum to exactly EQ target.
     * For best results, items should be in ascending order by weights.
     <pre>runtime complexity is O(n*capacity*max(quantities) where n = weights.length
     space complexity is O(capacity).
     </pre>
     * @param values values for items. values and weights describe the same items.
     * @param weights non-negative array of weights for items
     * @param quantities
     * @param target the exact sum that a combination of an unbounded number of weights should sum to
     * @return
     */
    public static int maxValueForTarget(int[] values, int[] weights, int[] quantities, int target) {
        return maxValue(values, weights, quantities, target, true);
    }

    public static int maxValue(int[] values, int[] weights, int[] quantities, int target, boolean solveForExact) {
        int n = weights.length;

        if (values.length != n) {
            throw new IllegalStateException("values and weights must be same length");
        }
        if (quantities.length != n) {
            throw new IllegalStateException("quantities and weights must be same length");
        }

        int[] prevTab = null;
        int[] currTab = new int[target + 1];
        int sentinel = Integer.MIN_VALUE;
        Arrays.fill(currTab, sentinel);
        currTab[0] = 0;

        int wc, q, wc2;
        for (int i = 0; i < n; ++i) {
            if (weights[i] > target) continue;

            prevTab = currTab;
            currTab = new int[currTab.length];
            Arrays.fill(currTab, sentinel);
            currTab[0] = 0;

            for (wc = 1; wc <= target; ++wc) {
                for (q = 1; q <= quantities[i]; ++q) {
                    wc2 = wc - q * weights[i];

                    //System.out.printf("i=%d, wc=%d, q=%d, wc2=%d, weights[i]=%d\n",
                    //        i, wc, q,  wc2, weights[i]);

                    if (wc2 == 0) {
                        currTab[wc] = Math.max(currTab[wc], q * values[i]);
                    } else if (wc2 > 0)  {
                        int s = (prevTab[wc2] == sentinel) ? sentinel : prevTab[wc2] + q * values[i];
                        currTab[wc] = Math.max(currTab[wc], s);
                    } else {
                        // carry forward best of prev and curr
                        currTab[wc] = Math.max(currTab[wc], prevTab[wc]);
                    }

                    //System.out.printf("    currTab[wc]=%d\n", currTab[wc]);
                }
            }
            //System.out.printf("%d) prevTab=\n%s\n", i, Arrays.toString(prevTab));
            //System.out.printf("%d) currTab=\n%s\n", i, Arrays.toString(currTab));
        }

        //if (prevTab != null)
        //    System.out.printf("\nprevTab=\n%s\n", Arrays.toString(prevTab));
        //System.out.printf("currTab=\n%s\n", Arrays.toString(currTab));

        if (solveForExact) {
            return currTab[target] == sentinel ? 0 : currTab[target];
        }

        // search backwards for max
        int max = currTab[target];
        int i = target;
        while (i > 0) {
            max = Math.max(max, currTab[--i]);
        }

        return max == sentinel ? 0 : max;
    }

    public static int minNumberOfItemsForTarget(int[] weights, int[] quantities, int target) {
        return minNumberOfItems(weights, quantities, target, true);
    }

    public static int minNumberOfItemsForCapacity(int[] weights, int[] quantities, int capacity) {
        return minNumberOfItems(weights, quantities, capacity, false);
    }

    public static int minNumberOfItems(int[] weights, int[] quantities, int target, boolean solveForExact) {

        int n = weights.length;
        int[] prevTab = null;
        int[] currTab = new int[target + 1];
        int sentinel = Integer.MAX_VALUE;
        Arrays.fill(currTab, sentinel);
        currTab[0] = 0;

        int wc, q, wc2;
        for (int i = 0; i < n; ++i) {
            if (weights[i] > target) continue;

            prevTab = currTab;
            currTab = new int[currTab.length];
            Arrays.fill(currTab, sentinel);
            currTab[0] = 0;

            for (wc = 1; wc <= target; ++wc) {
                for (q = 1; q <= quantities[i]; ++q) {
                    wc2 = wc - q * weights[i];

                    //System.out.printf("i=%d, wc=%d, q=%d, wc2=%d, weights[i]=%d\n",
                    //        i, wc, q,  wc2, weights[i]);

                    /*
                    for each item with wc2==0, replace current tab[wc] if LT
                                       wc2>0, add to previous remaining (wc2) if exists
                                       wc2<0, carry forward prev wc
                     */
                    if (wc2 == 0) {
                        currTab[wc] = Math.min(currTab[wc], q);
                    } else if (wc2 > 0)  {
                        int s = (prevTab[wc2] == sentinel) ? sentinel : q + prevTab[wc2];
                        currTab[wc] = Math.min(currTab[wc], s);
                    } else {
                        // carry forward best of prev and curr
                        currTab[wc] = Math.min(currTab[wc], prevTab[wc]);
                    }

                    //System.out.printf("    prevTab=%s\n", toString(prevTab));
                    //System.out.printf("    currTab=%s\n", toString(currTab));
                }
            }
        }
        //System.out.printf("currTab=%s\n", toString(currTab));

        if (solveForExact) {
            return currTab[target] == sentinel ? 0 : currTab[target];
        }

        // search backwards for last sum entered
        int i = target;
        while (i > 0 && currTab[i] == sentinel) {
            --i;
        }

        return currTab[i] == sentinel ? 0 : currTab[i];
    }

    public static int numberOfWaysForTarget(int[] weights, int[] quantities, int target) {

        int n = weights.length;

        if (weights.length != n) {
            throw new IllegalStateException("values and weights must be same length");
        }
        if (quantities.length != n) {
            throw new IllegalStateException("quantities and weights must be same length");
        }

        /*
        adapted from post by lee215
        https://leetcode.com/problems/number-of-ways-to-earn-points/solutions/3258120/java-c-python-knapsack-dp

        tab[wc] holds the number of ways to create sum wc from coins.

        NOTE that visiting wc (== the "current sum") in reverse order is needed to reuse sub-problems of smaller
        value wcs (== current sums) without having edited them for the current round of i
        It also means that one can break early when wc falls below the coin value types[i][1].

        tab[wc] holds the number of ways to create sum wc from coins.
        tab[wc] is initially [1, 0, 0, 0, ...] of length target+1.
        one complete round of i=0, tab[wc] is populated by coin [0] (types[0][1])
        at intervals for its value. since only tab[0] has a 1, the result for the first coin
        will be at intervals of coin[i] from index 0.
        e.g. for a coin[0] of value 2 and a quantity large enough to cover entire target, say 3 here,
        after one loop of i=0, tab will be:
        tab[wc] = 1 0 1 0 1 0,.... where the last 2 zeros are from coin[0]=2.
        The next coin adds to those non-zero positions at intervals of coin[i+1] (==types[i+1][1]).
        e.g. for coin[i=1] of value 5 and a quantity of at least 1, the result after the loop through i for i=1
        would be tab[wc] = 1 0 1 0 1 1,.... where the last zero is from coin[1]=5.
        so far, for target sum of 2, there is only 1 way to have formed it
        from coins [2,5] and quantities [3,1]
        and for target sum of 5, there is only 1 way to have formed it.
         */

        int[] tab = new int[target + 1];
        tab[0] = 1;
        for (int i = 0; i < n; ++i) {
            // since tab holds current and prev i results,
            // need to traverse weights from high to low
            // to avoid including an updated low wc2 in current wc
            for (int wc = target; wc >= weights[i]; --wc) {
                for (int q = 1; q <= quantities[i]; ++q) {
                    // wc2 is the remaining sum after q coins subtracted
                    int wc2 = wc - weights[i] * q;
                    if (wc2 < 0) break;
                    // adds counts from current sum and remaining sum
                    tab[wc] += tab[wc2];
                }
            }
        }
        //System.out.printf("tab=%s\n", Arrays.toString(tab));
        return tab[target];
    }

    private static String toString(int[][] tab) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < tab.length; ++i) {
            sb.append(String.format("%d) %s\n", i, Arrays.toString(tab[i])));
        }
        sb.delete(sb.length()-1, sb.length());
        return sb.toString();
    }

    private static String toString(int[] tab) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < tab.length; ++i) {
            if (tab[i] == Integer.MAX_VALUE) {
                sb.append("1000");
            } else if (tab[i] == Integer.MIN_VALUE) {
                sb.append("-1");
            } else {
                sb.append(tab[i]);
            }
            sb.append(",");
        }
        sb.delete(sb.length()-1, sb.length());
        return sb.toString();
    }
}
