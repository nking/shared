package algorithms.optimization;

import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 The knapsack problem or rucksack problem is a problem in combinatorial optimization:
 *     Given a set of items, each with a mass and a value,
 *     determine the number of each item to include in a collection
 *     so that the total weight is less than or equal to a given limit
 *     and the total value is as large as possible.
 *
 * The knapsack problem is also called Integer Linear Programming: given a
 * vector of variables and a set of linear constraints on the variables.
 *
 * ****
 * The runtime complexity of Knapsack 0 or 1 is O(n*W) where W is the capacity.
 * But if W is as large as n!, a brute force approach is better.
 * the brute force approach tries every permutation so RT is O(2^n).
 * n is the number of set items.
 *
 * Note that I adapted a version of the dynamic programming solution to use
 * intervals of weight equal the minimum weights instead of intervals of +1
 * when integrating up to the capacity, and that reduces the runtime by a factor
 * of roughly capacity/minimum(weight) as long a minimum weight .gt. 0.
 *
 * http://en.wikipedia.org/wiki/Knapsack_problem

 w_i is weight of each item
 v_i is value of each item
 W is the amount that the knapsack can hold.  (has to be an integer because
     it's used as a counter in code below.)
 n is the number of set items.

 Goal is to create a set S of items maximizing the total v for total w less than capacity W

Variants of Knapsack:

  0-1 knapsack problem:
     restricts the number v_i of copies of each kind of item to zero or one

  bounded knapsack problem (BKP):
     restricts the number of v_i copies to a quantity c_i

  unbounded knapsack problem (UKP):
     no restriction to the number of v_i copies


  There is a  pseudo-polynomial sumValues approach using dynamic programming for the 0-1 knapsack problem.
  There are also branch-and-bound approaches.  There is a naive brute force approach which is unfortunately O(2^n).

  wikipedia page was helpful in finding solutions and also
  http://www.cs.rit.edu/~zjb/courses/800/lec7.pdf

  Below is a dynamic programming solution for the 0-1 variant.

  The solution for one space in the memo matrix is composed of the
  memo from the (last weight bin - current item weight) plus potentially
  the addition of the current item value if it is larger current matrix space result.

* The solution below has runtime complexity O(n*W).
*
* The reason for runtime linear dependence on W is because W intervals of size 1 
* are used in the weight iteration.

* Note that if the capacity is very large, the user may want to renormalize
* the data to a smaller capacity as long as the items are still distinguishable
* from one another.
*
* @author nichole
*/
public class Knapsack01 {

    /**
     * solve the knapsack 0-1 problem with values and weights of items and capacity
     * for the knapsack.
     * Each item is put into the knapsack 0 or 1 times to maximize the total
     * value in the sack while keeping the total weight equal to or less than
     * capacity.
     * 
     * The worse case runtime complexity is O(n*W) where n is the number of
     * the items and W is the capacity.
     *
     * NOTE: consider reducing weights and capacity by the Greatest Common Denominator
     * between capacity and weights.
     * 
     @param values
     @param weights
     @param capacity
     @return 
     */
    public static int maxValueForCapacity(int[] values, int[] weights, int capacity) {

        int n = values.length;

        int[] prevTab = null;
        int[] currTab = new int[capacity + 1];

        int i, wc, c, wc2;
        for (i = 0; i < n; i++) {
            prevTab = currTab;
            currTab = new int[capacity + 1];
            for (wc = 1; wc <= capacity; wc++) {
                wc2 = wc-weights[i];
                // compare to previous wc.  0-1 means we cannot add to current wc2
                if (wc2 >= 0) {
                    currTab[wc] = Math.max(prevTab[wc], prevTab[wc2] + values[i]);
                } else {
                    currTab[wc] = prevTab[wc];
                }
            }
        }

        /*System.out.printf("weights=\n%s\n", FormatArray.toString(weights, "%d"));
        System.out.printf("values=\n%s\n", FormatArray.toString(values, "%d"));
        System.out.printf("memo=\n%s\n", FormatArray.toString(memo, "%d"));
        System.out.flush();*/

        return currTab[capacity];
    }

    /**
     * solve the knapsack 0-1 problem with values and weights of items and capacity
     * for the knapsack.
     * Each item is put into the knapsack 0 or 1 times to maximize the total
     * value in the sack while keeping the total weight equal to or less than
     * capacity.
     *
     * The worse case runtime complexity is O(n*W) where n is the number of
     * the items and W is the capacity.
     *
     * NOTE: consider reducing weights and capacity by the Greatest Common Denominator
     * between capacity and weights.
     *
     @param values
     @param weights
     @param capacity
     @return
     */
    public static int maxValueForCapacity2(int[] values, int[] weights, int capacity) {

        int n = values.length;

        int[] tab = new int[capacity + 1];

        int i, wc, t, c, wc2;
        for (i = 0; i < n; i++) {
            // since tab holds current and prev i results,
            // need to traverse weights from high to low
            // to avoid including an updated low wc2 in current wc
            for (wc = capacity; wc >= weights[i]; --wc) {
                wc2 = wc - weights[i];
                tab[wc] = Math.max(tab[wc], tab[wc2] + values[i]);
            }
        }

        return tab[capacity];
    }

    /**
     * find maximum value or sum of items whose weights sum to exactly equal the target, but with the restricted quantity
     * of 0 or 1 for each item.
     * @param values non-negative values for items. values and weights describe the same items.
     * @param weights non-negative array of weights for items
     * @param target the exact sum of weights for the minimum number of items which sum to that
     *               weight.  each item can only be added up to a quantity of 1.
     * @return the maximum value for items whose weights sum up to exactly equal target.
     */
    public static int maxValueForTarget(int[] values, int[] weights, int target) {

        int n = values.length;

        int[] prevTab = null;
        int[] currTab = new int[target + 1];
        int sentinel = Integer.MIN_VALUE;
        Arrays.fill(currTab, sentinel);
        currTab[0] = 0;

        int i, wc,c, wc2;
        for (i = 1; i <= n; i++) {
            prevTab = currTab;
            currTab = new int[target + 1];
            Arrays.fill(currTab, sentinel);
            currTab[0] = 0;

            for (wc = 1; wc <= target; wc++) {
                wc2 = wc-weights[i-1];

                if (wc2 == 0) {
                    currTab[wc] = Math.max(prevTab[wc],  values[i - 1]);
                } else if (wc2 > 0)  {
                    int s = (prevTab[wc2] == sentinel) ? sentinel : prevTab[wc2] + values[i-1];
                    currTab[wc] = Math.max(prevTab[wc], s);
                } else {
                    currTab[wc] = prevTab[wc];
                }

            }
        }

        return currTab[target] == sentinel ? 0 : currTab[target];
    }

    /**
     * find maximum value or sum of items whose weights sum to exactly equal the target, but with the restricted quantity
     * of 0 or 1 for each item.
     * @param values non-negative values for items. values and weights describe the same items.
     * @param weights non-negative array of weights for items
     * @param target the exact sum of weights for the minimum number of items which sum to that
     *               weight.  each item can only be added up to a quantity of 1.
     * @return the maximum value for items whose weights sum up to exactly equal target.
     */
    public static int maxValueForTarget2(int[] values, int[] weights, int target) {

        int n = values.length;

        int[] tab = new int[target + 1];
        int sentinel = Integer.MIN_VALUE;
        Arrays.fill(tab, sentinel);
        tab[0] = 0;

        int i, wc, c, wc2;
        for (i = 0; i < n; i++) {
            // since tab holds current and prev i results,
            // need to traverse weights from high to low
            // to avoid including an updated low wc2 in current wc
            for (wc = target; wc >= weights[i]; --wc) {
                wc2 = wc-weights[i];
                if (wc2 == 0) {
                    tab[wc] = Math.max(tab[wc],  values[i]);
                } else if (wc2 > 0)  {
                    int s = (tab[wc2] == sentinel) ? sentinel : tab[wc2] + values[i];
                    tab[wc] = Math.max(tab[wc], s);
                }
            }
        }

        return tab[target] == sentinel ? 0 : tab[target];
    }

    /**
     * find the minimum number of items whose weights sum to LEQ capacity.
     * @param weights non-negative array of item weights
     * @param capacity the largest sum of weights for the minimum number of items which sum to that
     * weight.  each item can only be added up to a quantity of 1.  if capacity isn't achieved, the solution
     *                 for the next lowest summed amount is returned.
     * @return
     */
    public static int minNumberOfItemsForCapacity(int[] weights, int capacity) {
        return minNumberOfItems(weights, capacity, false);
    }

    /**
     * find the minimum number of items whose weights sum to exact target.
     * @param weights non-negative array of item weights
     * @param target the sum of weights for the minimum number of items which sum to that
     * target.
     * @return
     */
    public static int minNumberOfItemsForTarget(int[] weights, int target) {
        return minNumberOfItems(weights, target, true);
    }

    public static int minNumberOfItems(int[] weights, int target, boolean solveForExact) {
        int n = weights.length;

        // tab[wc] holds the min number of items for the items whose weights sum to wc
        int[] prevTab = null;
        int[] currTab = new int[target + 1];
        int sentinel = Integer.MAX_VALUE;
        Arrays.fill(currTab, sentinel);
        currTab[0] = 0;

        int i, wc, wc2;
        for (i = 0; i < n; i++) {
            if (weights[i] > target) continue;
            prevTab = currTab;
            currTab = new int[target + 1];
            Arrays.fill(currTab, sentinel);
            currTab[0] = 0;

            for (wc = 1; wc <= target; wc++) {
                wc2 = wc - weights[i];

                if (wc2 == 0) {
                    currTab[wc] = Math.min(currTab[wc], 1);
                } else if (wc2 > 0)  {
                    int s = (prevTab[wc2] == sentinel) ? sentinel : 1 + prevTab[wc2];
                    currTab[wc] = Math.min(currTab[wc], s);
                } else {
                    currTab[wc] = prevTab[wc];
                }

                //System.out.printf("i=%d, wc=%d, wc2=%d, weights[i]=%d  tab=%s\n", i, wc, wc2, weights[i], Arrays.toString(tab));
            }
        }
        //System.out.printf("tab=%s\n", Arrays.toString(tab));
        if (solveForExact) {
            return currTab[target] == sentinel ? 0 : currTab[target];
        }

        // search backwards for last sum entered
        i = target;
        while (i > 0 && currTab[i] == sentinel) {
            --i;
        }

        return currTab[i] == sentinel ? 0 : currTab[i];
    }

    public static int minNumberOfItemsForCapacity2(int[] weights, int capacity) {
        return minNumberOfItems2(weights, capacity, false);
    }

    public static int minNumberOfItemsForTarget2(int[] weights, int target) {
        return minNumberOfItems2(weights, target, true);
    }

    public static int minNumberOfItems2(int[] weights, int capacity, boolean solveForTarget) {
        int n = weights.length;

        // tab[wc] holds the min number of items for the items whose weights sum to wc
        int[] tab = new int[capacity + 1];
        int sentinel = Integer.MAX_VALUE;
        Arrays.fill(tab, sentinel);
        tab[0] = 0;

        int i, wc, wc2;
        for (i = 0; i < n; i++) {
            // since tab holds current and prev i results,
            // need to traverse weights from high to low
            // to avoid including an updated low wc2 in current wc
            for (wc = capacity; wc >= weights[i]; --wc) {
                wc2 = wc - weights[i];
                if (wc2 == 0) {
                    tab[wc] = Math.min(tab[wc], 1);
                } else if (wc2 > 0)  {
                    int s = (tab[wc2] == sentinel) ? sentinel : 1 + tab[wc2];
                    tab[wc] = Math.min(tab[wc], s);
                }
            }
        }

        i = capacity;
        if (solveForTarget) {
            return tab[i] == sentinel ? 0 : tab[i];
        }

        while (i > 0 && tab[i] == sentinel) {
            --i;
        }
        //System.out.printf("tab=%s\n", Arrays.toString(tab));
        return tab[i] == sentinel ? 0 : tab[i];
    }

    public static int numberOfWaysForTarget(int[] weights, int target) {

        int n = weights.length;

        if (weights.length != n) {
            throw new IllegalStateException("values and weights must be same length");
        }

        int[] tab = new int[target + 1];
        tab[0] = 1;
        int wc, wc2;
        for (int i = 0; i < n; ++i) {
            // since tab holds current and prev i results,
            // need to traverse weights from high to low
            // to avoid including an updated low wc2 in current wc
            for (wc = target; wc > 0; --wc) {
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
     * a quick look at integrating over weight intervals of the size of the smallest
     * weight up to capacity (in contrast to +1 intervals in standard impl).
     * It's scaling the data to reduce the runtime complexity, but at the
     * expense of accuracy.
     @param values
     @param weights
     @param capacity
     @return 
     */
    static int approxDynamically(int[] values, int[] weights, int capacity) {

        int factor;
        int q1 = findQ1Diff(weights);
        if (q1 < 2) {
            factor = 1;
        } else {
            factor = q1;
        }
        //int minW = MiscMath0.findMin(weights);
        //factor = minW;
        
        int n = values.length;
        int i, j, t, c, t2, wc, w;
        
        int[] ws = getIntervals(capacity, factor);
        
        //System.out.printf("capacity=%d, factor=%d, q1Diff=%d, intervals=%s\n", capacity,
        //    factor, q1, Arrays.toString(ws));
                
        int[][] memo = new int[n + 1][ws.length + 1];
        for (i = 0; i < memo.length; i++) {
            memo[i] = new int[ws.length  + 1];
        }
        
        for (i = 1; i <= n; i++) {
            for (j = 0; j < ws.length; ++j) {//10,8,6,4,2
                w = ws[j];//capacities
                wc = j+1;
                t = memo[i - 1][wc];
                t2 = (int) Math.ceil((float) weights[i - 1] / (float) factor);
                if (w >= weights[i - 1]) {// to avoid exceeding capacity
                    c = memo[i - 1][wc - t2] + values[i - 1];
                    memo[i][wc] = Math.max(c, t);
                } else {
                    memo[i][wc] = t;
                }
            }
        }

        //System.out.printf("memo=\n%s\n", FormatArray.toString(memo, "%d"));
        
        return memo[memo.length - 1][memo[0].length - 1];
    }

    private static int[] getIntervals(int capacity, int minW) {
        
        TIntList wcs = new TIntArrayList();
        for (int wc = capacity; wc >=1; wc-=minW) {
            wcs.add(wc);
        }
        wcs.reverse();  
        
        return wcs.toArray();
    }

    private static int findQ1Diff(int[] a) {
        a = Arrays.copyOf(a, a.length);
        Arrays.sort(a);
        int[] diffs = new int[a.length-1];
        int i, d;
        int minDiff = Integer.MAX_VALUE;
        for (i = 1; i < a.length; ++i) {
            d = a[i] - a[i - 1];
            if (d < minDiff) {
                minDiff = d;
            }
            diffs[i-1] = d;
        }
        int q1 = diffs[(int)(diffs.length*.25)];
        //System.out.printf("median diff=%d  Q1 diff = %d\n", diffs[diffs.length/2],
        //    q1);
        return q1;
    }

}
