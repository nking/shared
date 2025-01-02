package algorithms.puzzles;

import junit.framework.TestCase;

public class TriplestepDPTest extends TestCase {

    /*
    problem 8.1 from Cracking the Coding Interview book by Gayle McDowell

    Triple Step: A child is running up a staircase with n steps and can hop either 1 step, 2 steps, or 3 steps at a time.
    Implement a method to count how many possible ways the child can run up the stairs.

    Note that if one could take any number of steps between 1 and n,
    the answer would be by composition: 2^(n-1).

    The author estimates the recursive solution as having
    r.t.c. as O(3^n), but a sketch of the recursion tree,
    can be summarized as O(((1<<n) + (1<<(n-1)) + (1<<(n-2)))
    which closely approximates the nIter tracked below.
     */

    int nIter = 0;
    // my recursive solution
    public long nWaysR(int n) {
        nIter = 0;
        // [0,n-1]
        return r(0, n);
    }//
    public long r(int i, int n) {
        ++nIter;
        if (i == n) {
            return 1;
        }
        long c = 0;
        for (int j = i+1; j <= i+3; ++j) {
            if (j > n) continue;
            c += r(j, n);
        }
        return c;
    }

    // my tabular solution
    // r.t.c. O(n)
    public long nWaysTabular(int n) {
        // keep adding the previous 3 numbers
        long[] tab = new long[]{1,2,4};// n=1, n=2, n=3
        if (n < 4) {
            return tab[n-1];
        }
        // 5 memory loads, 3 memory stores, 3 flops
        for (int i = 4; i <= n; ++i) {
            long s = tab[2] + tab[1] + tab[0];
            tab[0] = tab[1];
            tab[1] = tab[2];
            tab[2] = s;
        }
        return tab[2];
        /*
        3 memory loads, 3 memory stores, + 3*2 flops (includes index math)
        so this method is not better.
        int j0 = 0;
        for (int i = 4; i <= n; ++i) {
            // add the other 2 to tab[j0];
            for (int j = 1; j < 3; ++j) {
                tab[j0] += tab[(j0 + j)%3];
            }
            // increment j0
            j0 = (j0 + 1)%3;
        }
        return tab[(j0+2)%3];
         */
    }

    // the book solution
    public long countWays(int n) {
        ++nIter;
        if (n < 0) {
            return 0;
        } else if (n == 0) {
            return 1;
        } else {
            return countWays(n - 1) + countWays(n - 2) + countWays(n - 3);
        }
    }
    public void test0()  {
        int[] costs = new int[3];
        int[][] c = new int[][] {
            {1, costs[0]}, {7, costs[1]}, {30, costs[2]}
        };
        int rtc;
        for (int n = 3; n < 20; ++n) {
            rtc = (1<<n) + (1<<(n-1)) + (1<<(n-2));
            long c1 = nWaysR(n);
            long c2 = nWaysTabular(n);
            nIter = 0;
            long c3 = countWays(n);
            System.out.printf("n=%2d  c1=%4d (nIter=%4d)  c2=%4d  c3=%4d (nIter=%4d); recursive r.t.c. est =%d\n",
                    n, c1, nIter, c2, c3, nIter, rtc);
        }
    }

}
