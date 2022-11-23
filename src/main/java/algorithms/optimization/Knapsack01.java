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

 Goal is to create a set S of items maximizing the total v for total w less than capacity W

Variants of Knapsack:

  0-1 knapsack problem:
     restricts the number v_i of copies of each kind of item to zero or one

  bounded knapsack problem (BKP):
     restricts the number of v_i copies to a quantity c_i

  unbounded knapsack problem (UKP):
     no restriction to the number of v_i copies


  There is a  pseudo-polynomial sumValues using dynamic programming for 0-1 knapsack problem.
  There are also branch-and-bound approaches.  There is a naive brute force approach which is unfortunately O(2^n).

  wikipedia page was helpful in finding solutions and also
  http://www.cs.rit.edu/~zjb/courses/800/lec7.pdf

  Below is a dynamic programming solution for the 0-1 variant.

  the solution for one space in the memo matrix is composed of the
  memo from the (last weight bin - current item weight) plus potentially
  the addition of the current item value if it is larger current matrix space result.

* The solution below has runtime complexity O(n*W).
*
* The W is a factor because W intervals of size 1 are used in the weight iteration.

* Note that if the capacity is very large, the user may want to renormalize
* the data to a smaller capacity as long as the items are still
* differentiable from one another.
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
     * If the minimum value in the weights array is larger than 1, it is used
     * as an interval in summing up to the capacity, hence reduces the runtime
     * complexity by roughly that as a factor.
     * 
     * @param values
     * @param weights
     * @param capacity
     */
    public static int solveDynamically(int[] values, int[] weights, int capacity) {

        return solveDynamically0(values, weights, capacity);
    }
    
    static int solveDynamically0(int[] values, int[] weights, int capacity) {

        int n = values.length;
        
        int[][] memo = new int[n + 1][capacity + 1];
        for (int i = 0; i < memo.length; i++) {
            memo[i] = new int[capacity + 1];
        }
        int i, wc, t, c;
        for (i = 1; i <= n; i++) {
            for (wc = 1; wc <= capacity; wc++) {
                t = memo[i-1][wc];
                if (wc >= weights[i-1]) { // to avoid exceeding capacity
                    c = memo[i-1][wc-weights[i-1]] + values[i-1];
                    memo[i][wc] = Math.max( c, t);
                } else {
                    memo[i][wc] = t;
                }
            }
        }
        
        System.out.printf("memo=%n%s%n", FormatArray.toString(memo, "%d"));

        return memo[memo.length - 1][memo[0].length - 1];
    }
    
    /**
     * a quick look at integrating over weight intervals of the size of the smallest
     * weight up to capacity (in contrast to +1 intervals in standard impl).
     * It's scaling the data to reduce the runtime complexity, but at the
     * expense of accuracy.
     * @param values
     * @param weights
     * @param capacity
     * @return 
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
        
        System.out.printf("capacity=%d, factor=%d, q1Diff=%d, intervals=%s%n", capacity, 
            factor, q1, Arrays.toString(ws));
                
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

        System.out.printf("memo=%n%s%n", FormatArray.toString(memo, "%d"));
        
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
        System.out.printf("median diff=%d  Q1 diff = %d%n", diffs[diffs.length/2],
            q1);
        return q1;
    }

}
