package algorithms.exponential;

import algorithms.util.FormatArray;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * a look at solving TSP dynamically. The algorithms is NP-Hard, no known
 * polynomial time algorithms for TSP, so an approximation such as TSP-Prim's
 * MST should be used or one of the 1.5*optimal algorithms (see wikipedia).
 *
 * This dynamic version is just to be instructive, not to be used as the
 * runtime complexity is O(n^2 * 2^n) where n is the number of vertices.
 *
 * The implementation is from
 * https://www.interviewbit.com/blog/travelling-salesman-problem/
 * There is a copyright on the page, so one should assume the material is
 * not open source.
 * For that reason, this class is packaged in the test code, not part of
 * the project api jar.
 * 
 *
 * @author nichole
 */
public class TSPDynamic {

    private final int N, start;
    private final double[][] distance;
    private List<Integer> tour = new ArrayList<>();
    private double minTourCost = Double.POSITIVE_INFINITY;
    private boolean solverFinished = false;

    public TSPDynamic(double[][] distance) {
        this(0, distance);
    }

    public TSPDynamic(int start, double[][] distance) {
        
        N = distance.length;

        if (N <= 2) {
            throw new IllegalStateException("N <= 2 not yet supported.");
        }
        if (N != distance[0].length) {
            throw new IllegalStateException("Matrix must be square (n x n)");
        }
        if (start < 0 || start >= N) {
            throw new IllegalArgumentException("Invalid start node.");
        }
        this.start = start;
        this.distance = distance;
    }

    // Returns the optimal tour for the traveling salesman problem.
    public List<Integer> getTour() {
        if (!solverFinished) {
            solveIteratively();
        }
        return tour;
    }
    // Returns the minimal tour cost.

    public double getTourCost() {
        if (!solverFinished) {
            solveIteratively();
        }
        return minTourCost;
    }

    // Solves the traveling salesman problem and caches solution.
    public void solveIteratively() {

        if (solverFinished) {
            return;
        }

        final int END_STATE = (1 << N) - 1;
        double[][] memo = new double[N][END_STATE + 1];

        // Add all outgoing edges from the starting node to memo table.
        for (int end = 0; end < N; end++) {
            if (end == start) {
                continue;
            }
            memo[end][(1 << start) | (1 << end)] = distance[start][end];
            /*
            start  end    1 << start)  (1 << end)   (1 << start) | (1 << end)
            0      0      1            1             1
            0      1      1            2             3
            0      2      1            4             5
            0      3      1            8             9
            */
        }
        
        System.out.printf("memo=\n%s\n", FormatArray.toString(memo, "%.1f"));

        // r is the number of vertexes within n vertexes, in which the subset bits are set to 1.
        for (int r = 3; r <= N; r++) {
            for (int subset : combinations(r, N)) {
                
                System.out.printf("r=%d subset=%d(%s)\n", r, subset, Integer.toBinaryString(subset));
                
                if (notIn(start, subset)) {
                    continue;
                }
                for (int next = 0; next < N; next++) {
                    System.out.printf("   next=%d(%s)", next, Integer.toBinaryString(next));
                    if (next == start || notIn(next, subset)) {
                        System.out.printf("\n");
                        continue;
                    }
                    int subsetWithoutNext = subset ^ (1 << next);
                    
                    System.out.printf("   subsetWithoutNext=%d(%s)\n", subsetWithoutNext, Integer.toBinaryString(subsetWithoutNext));
                    
                    double minDist = Double.POSITIVE_INFINITY;
                    for (int end = 0; end < N; end++) {
                        if (end == start || end == next || notIn(end, subset)) {
                            continue;
                        }
                        double newDistance = memo[end][subsetWithoutNext] + distance[end][next];
                        if (newDistance < minDist) {
                            minDist = newDistance;
                            System.out.printf("      end=%d(%s) minDist=%.1f\n", end, Integer.toBinaryString(end), minDist);
                        }
                    }
                    memo[next][subset] = minDist;
                }
            }
        }

        // Connect tour back to starting node and minimize cost.
        for (int i = 0; i < N; i++) {
            if (i == start) {
                continue;
            }
            double tourCost = memo[i][END_STATE] + distance[i][start];
            if (tourCost < minTourCost) {
                minTourCost = tourCost;
                System.out.printf("      i=%d minTourCost=%.1f\n", i, minTourCost);
            }
        }
        int lastIndex = start;
        int state = END_STATE;
        tour.add(start);

        // Reconstruct TSP path from memo table.
        for (int i = 1; i < N; i++) {

            int index = -1;
            for (int j = 0; j < N; j++) {
                if (j == start || notIn(j, state)) {
                    continue;
                }
                if (index == -1) {
                    index = j;
                }
                double prevDist = memo[index][state] + distance[index][lastIndex];
                double newDist = memo[j][state] + distance[j][lastIndex];
                if (newDist < prevDist) {
                    index = j;
                }
            }

            tour.add(index);
            state = state ^ (1 << index);
            lastIndex = index;
        }

        tour.add(start);
        Collections.reverse(tour);

        solverFinished = true;
    }

    // test bit operation to see if elem is in sbset
    private static boolean notIn(int elem, int subset) {
        return ((1 << elem) & subset) == 0;
    }

    // This method generates all bit sets of size n where r bits 
    // are set to one. The result is returned as a list of integer masks.
    //TODO: improve this
    public static List< Integer> combinations(int r, int n) {
        List<Integer> subsets = new ArrayList<>();
        combinations(0, 0, r, n, subsets);
        return subsets;
    }

    private static void combinations(int set, int at, int r, int n, List< Integer> subsets) {

        // Return early if there are more elements left to select than what is available.
        int elementsLeftToPick = n - at;
        if (elementsLeftToPick < r) {
            return;
        }

        // We selected 'r' elements so we found a valid subset!
        if (r == 0) {
            subsets.add(set);
        } else {
            for (int i = at; i < n; i++) {
                // Try including this element
                set |= 1 << i; // set bit operation sets bit i in 'set'

                combinations(set, i + 1, r - 1, n, subsets);

                set &= ~(1 << i); // clear bit operation, unsets bit i in 'set'
            }
        }
    }
    
    //===================================
    // https://www.interviewbit.com/blog/travelling-salesman-problem/
      private int[] completed = null;
    private int sentinel = Integer.MAX_VALUE;

    public void solveRecursively() {
        this.completed = new int[N];
        tour.clear();
        minTourCost = 0;
        solverFinished = false;

        mincost(start);

        solverFinished = true;

        System.out.printf("\ntour cost=%.0f\n", minTourCost);
    }

    private void mincost(int city) {
        int ncity;

        completed[city] = 1;

        System.out.printf("%d--->", city);
        tour.add(city);
        ncity = least(city);

        if (ncity == sentinel) {
            ncity = 0;
            System.out.printf("%d", ncity);
            minTourCost += distance[city][ncity];
            tour.add(ncity);

            return;
        }

        mincost(ncity);
    }

    private int least(int c) {
        int i, nc = sentinel;
        double min = sentinel;
        double kmin = sentinel;

        for (i = 0; i < N; i++) {
            if ((distance[c][i] != 0) && (completed[i] == 0)) {
                if (distance[c][i] + distance[i][c] < min) {
                    min = distance[i][0] + distance[c][i];
                    kmin = distance[c][i];
                    nc = i;
                }
            }
        }

        if (min != sentinel) {
            minTourCost += kmin;
        }

        return nc;
    }
}
