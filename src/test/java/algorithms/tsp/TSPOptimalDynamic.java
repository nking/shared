package algorithms.tsp;

import algorithms.util.FormatArray;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Dynamic algorithms for the traveling Salesman problem.. 
 * The algorithms is NP-Hard, no known
 * polynomial time algorithms exist for TSP, so an approximation such as TSP-Prim's
 * MST should be used or one of the 1.5*optimal algorithms (see wikipedia).
 *
 * This dynamic version is instructive, and not to be used for most datasets as the
 * runtime complexity is O(n^2 * 2^n) where n is the number of vertices.
 * 
 * The iterative method is limited to data holding 31 cities or less.
 * 
 *  The approximate TSP algorithms are in the curvature scale space project.

 * <pre>
 * references:
 * 
 * The implementation is from
 * https://www.interviewbit.com/blog/travelling-salesman-problem/
 * There is a copyright on the page, so one should assume the material is
 * not open source.
 * For that reason, this class is packaged in the test code, not part of
 * the project api jar.
 * 
 * A variant with a cost limit rather than minimizaing the cost is
 * in the thesis of Sokkappa, 1990.
 * "The cost-constrained traveling salesman problem"
 * https://doi.org/10.2172/6223080
 * 
 * 
 * The thesis refers to the
 * Held–Karp algorithm a.k.a. Bellman–Held–Karp algorithm
 * https://en.wikipedia.org/wiki/Held%E2%80%93Karp_algorithm
 * pp 58 - 
 * 
 * for start city = 1,
 * calculate for each set of cities  S={2,3,4,...n} and every city
 * t not contained in S and not being 1
 * the shortest one-way path from 1 to t that passes through every city in 
 * S in some order (but not through any other cities). 
 * Denote this distance g(S, t), 
 * and write d(u,v) for the length of the direct edge from u to v.
 * Compute values of g(S, t) starting with the smallest sets S and finishing
 * with the largest.
 * 
 * When {S has two or fewer elements, then calculating g(S, t)}
 * requires looking at one or two possible shortest paths. 
 * For example, g(emptyset, t) is simply d(1, t), and g(2, 3) 
 * is just the length of 1-->2-->3.
 * Likewise, g({2, 3}, 4) is the length of either 1-->2-->3-->4 OR
 * 1-->3-->2-->4, whichever is shorter.
 * 
 * Once S contains 3 or more cities, the number of paths thru S rises quickly,
 * but only a few such paths need to be examined to find the shortest.
 * For example, if 1-->2-->3-->4 is shorter than 1-->3-->2-->4, 
 * then 1-->2-->3-->4-->5 mst be shorter than 1-->3-->2-->4-->5
 * and the length of the later is not a possible value of g({2, 3, 4}, 5).
 * 
 * Similarly, if the shortest path from 1 to {2,3,4} to 5 is 1-->4-->3-->2-->5
 * and the shortest path from 1 to {2,3,4,5} to 6 ends with the edge 5-->6
 * then the whole path 1 to 6 must be 1-->4-->3-->2-->5-->6 and not any
 * other of the 5 paths created by visiting {2,3,4} in different order.
 * 
 * More generally, let the set of k cities be S = {s1, s2, s3, s4, ...sk).
 * For every integer 1 .leq. i .leq. k, write 
 * S_i is S with si removed: {s1, s2, ...s(i-1), s(i+1), ..sk).
 * Then if the shortest path from 1 to S through t
 * has si as its 2nd to last city, then removing the final edge from this path
 * must give the shortest path from 1 to si thru S_i.
 * This means that there are only k possible shortest paths from 1 to t thru S,
 * one for each possible 2nd to last city si with length g(S_i, si) + d(si, t)
 * and g(S, t) = min_{1 .leq. i .leq. k}( g(S_i, si) + d(si, t) ).
 * 
 * This stage of the algorithm finishes when g({2, ... ,i-1,i+1,n}, i)
 * is known for every integer 2 .leq. i .leq. n
 * giving the shortest distance from city 1 to city i that passes through every
 * other city. 
 * 
 * The much shorter second stage adds these distances to the edge lengths 
 * d(i,1) to give n-1 possible shortest cycles, and then finds the shortest.

 * 
 * runtime complexity:
 * 
 * Computing one value of g(S, t) for a k-element subset S of {2, ...n}
 * requires finding the shortest of k possible paths, each found by adding a 
 * known value of g and an edge length from the original graph; that is, 
 * it requires time proportional to k.
 * 
 * There are C(n-1, k) k-element subsets of {2,...n} and each subset gives
 * n-k-1 possible values of t.
 * 
 * Computing all values of G(S, t) where |S| = k, 
 * requires time k *(n-k-1)*C(n-1, k) = (n-1)*(n-2) * C(n-3, k-1).
 * Then the total time across all subset sizes is 
 * (n-1)*(n-2) * summation_over_k_from_1_to_{n-2}( C(n-3, k-1) )
 * =  (n-1)*(n-2) * 2^(n-3).
 * The runtime complexity for Held–Karp algorithm is then O(n^2 * 2^n).
 * The 2nd stage is O(n) and does affect the runtime complexity.
 * 
 * 
 * NOTE that the method below has a larger runtime complexity for now as
 * some of the inner loops need to be edited to visit only the nodes not
 * already in the current solution path.
 * 
 * <pre>
 * NK:
 * Keeping notes here on looking at bitstring patterns to find ways to 
 * condense the sums of 3-paths used in the subproblems.
 * 
    let start node = 0.  n=5
    The total number of permutations of a path with fixed start node = 0 and n=5
    is n!/n = 24.

    for n=127, nPerm = 2.4e+211, so if wanted to store each permutation, indexed,
    would need to use java's BigInteger as a bitstring or this project's VeryLongBitString.java

    to re-use sub-problems, can determine how many sets of paths of 3 nodes are in
    the set of n-1 which is (n-1)/3.
    (sets of size 2 are given by the distance matrix.)

    Store all 3 node combinations, excluding the number 0, for reuse.
    The number of k-permutations without excluding a node is: C(n,k) = n!/(k!*(n-k)!).
    For example, 4 nodes, k-permutations size=3:  C(n,k) = n!/(k!*(n-k)!) = 4
      7 (    111)
     11 (   1011)
     13 (   1101)
     14 (   1110)
    Excluding the permutations with bit=0 set: use n2=n-1, C(n-1,k) = 1
     14 (   1110)

    So, one can generate the k=3 permutations of subset of size n, and exclude the 1st bit
    by generating for n2=n-1, and the left shift by 1 the result.
    (can use SubsetChooser.java for the permutations)

    To further look at bit patterns to see if memoization can be used, but with
    condensed intermediate steps to use less storage, one more example:

    For n=10, nPerm=84:
    as the 3 sequence bitstrings are generated, the inverse bitstring w/ 0 cleared
    has the set of nodes this bitstring can be added to
         bitstring    inverse, excl 0
    14 (0000001110)  (1111110000)
    22 (0000010110)  (1111101000) perm for n3=125 can be generated, exclude bit 1, then left shift by 3 to get paths to add this to.
    26 (0000011010)  (1111100100) perm for n3=249 can be generated, exclude bits 1,2, then left shift by 2 to get paths to add this to.
    28 (0000011100)  (1111100000)
    38 (0000100110)  (1111011000)
    42 (0000101010)  (1111010100)
    44 (0000101100)  (1111010000)
    50 (0000110010)  (1111001100)
    52 (0000110100)  (1111001000)
    56 (0000111000)  (1111000100)
    70 (0001000110)  (1110111000)
    74 (0001001010)  (1110110100)
    ...
    
    Once all of the 3-node path subsets are generated and stored,
    one can consider again that each of the complete permutations of paths 
    is composed of disjoint combinations of (n-1)/3 of the 3-node path subsets,
    (plus up to 2 nodes if n = 12, we have (n-1)/3 3-node path subsets 
     + 1 fixed node + 2 free nodes that belong in the permuation).

    For now, will consider only the cases for n = 1 + a multiple of 3 to look at
    ways to condense the problem.
    
    have 84 k=3 sequences and have the unset bits, excluding bit 0
    
    note: as each complete path of the(n-1)/3 sections is totaled, compare it to minTotal.

paused here
n=10, k=3
ni = n-1
for (p = 0; p &lt np-1; ++p) { 0,1 (the subsetchooser needs n>k)
    subsetchooser(ni, k) 
    nSetBits = ni-3 
    store to reuse... see TSPDynamicOutlineTest
    ni -= 3;    
}                   
//nSetBits = n-1-3-3
//X--> stop before this as n-1-6 is not &gt 6 subsetchooser(n-1-6, k)

    iterate over the 84 k=3 sequences:
        have si     and siInv
        (0000010110)    (1111101000)
        only add si to permutations of the set bits: 1111101000
           *each of those is composed of 3-paths already calculated,
           but as each permutation is generated, and si added to it,
           need to store that for reuse.

           consider the datastructure needed for storage:               
               set(bitstring, sum), get(bitstring), contains(bitstring).
           the bitstring keys might be very large (BigInteger or VeryLongBitString).
           the number of items might be very large.
           java language restricts array sizes to (1 &lt; &lt; 31)-1 and that's true for the
           structures internal to the collections.
           the total number of permutations being n!/n means can only solve for 13
           cities if need to store each permutation.
           will revisit this after outline...
           
    pausing here...

 * </pre>
 * 
 * </pre>
 * @author nichole
 */
public class TSPOptimalDynamic {

    private final int N, start;
    private final double[][] distance;
    private List<Integer> tour = new ArrayList<>();
    private double minTourCost = Double.POSITIVE_INFINITY;
    private boolean solverFinished = false;
    
    private static Logger logger;
    private static Level LEVEL = Level.FINE;

    /**
     * construct the TSP instance using the distance matrix.  the starting
     * vertex will always be 0, as the start node does not matter for a
     * solution including a cycle with all vertexes.
     * 
     * @param distance 
     */
    public TSPOptimalDynamic(double[][] distance) {
        this(0, distance);
    }
    
    public TSPOptimalDynamic(int start, double[][] distance) {
        N = distance.length;

        if (N <= 2) {
            throw new IllegalStateException("N <= 2 not yet supported.");
        }
        if (N != distance[0].length) {
            throw new IllegalStateException("Matrix must be square (n x n)");
        }
        this.start = start;
        this.distance = distance;
        
        logger = Logger.getLogger(this.getClass().getSimpleName());
        
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

    /** Solves the traveling salesman problem and caches solution.
     * NOTE that algorithm can only handle N less than 32.
     * 
     * TODO: should be able to re-order the k-permutations and visits
        so that the embedded 0 to N loops are over the unvisited vertexes only,
        hence decreasing the loop iteratively.
     */
    public void solveIteratively() {
        
        if (N > 31) {
            throw new IllegalArgumentException("the dynamic iterative algorithm "
                + "has a limit of 31 for the number of cities.");
        }

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
        logger.log(LEVEL, String.format(
        //System.out.printf(
            "memo[%d][%d]=\n%s\n", memo.length, memo[0].length, 
            FormatArray.toString(memo, "%.1f")));
        
        //TODO: should be able to re-order the k-permutations and visits
        // so that the embedded 0 to N loops are over the unvisited vertexes only,
        //    hence decreasing the loop iteratively.

        // r is the number of vertexes within n vertexes, in which the subset bits are set to 1.
        for (int r = 3; r <= N; r++) {
            
            for (int subset : combinations(r, N)) {
                
                logger.log(LEVEL, String.format(
                //System.out.printf(
                    "r=%d subset=%d(%s)\n", r, subset, Integer.toBinaryString(subset)));
                
                if (notIn(start, subset)) {
                    continue;
                }
                                
                for (int next = 0; next < N; next++) {
                    logger.log(LEVEL, String.format(
                        //System.out.printf(
                        "   next=%d(%s)", next, Integer.toBinaryString(next)));
                    if (next == start || notIn(next, subset)) {
                        logger.log(LEVEL, String.format(
                            //System.out.printf(
                            "\n"));
                        continue;
                    }
                    
                    int subsetWithoutNext = subset ^ (1 << next);
                    
                    logger.log(LEVEL, String.format(
                    //System.out.printf(
                        "   subsetWithoutNext=%d(%s)\n", subsetWithoutNext, 
                        Integer.toBinaryString(subsetWithoutNext)));
                    
                    double minDist = Double.POSITIVE_INFINITY;
                    for (int end = 0; end < N; end++) {
                        if (end == start || end == next || notIn(end, subset)) {
                            continue;
                        }
                        double newDistance = memo[end][subsetWithoutNext] + distance[end][next];
                        if (newDistance < minDist) {
                            minDist = newDistance;
                            logger.log(LEVEL, String.format(
                            //System.out.printf(
                               "      end=%d(%s) "
                               + "minDist = memo[end %d][subWONxt %d] + distance[end %d][next %d] = %.1f\n", 
                               end, Integer.toBinaryString(end), end, subsetWithoutNext, end, next, minDist));
                        }
                    }
                    logger.log(LEVEL, String.format(
                    //System.out.printf(
                        "      stored as memo[next %d][subset %d]=%.1f\n",
                        next, subset, minDist));
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
                logger.log(LEVEL, String.format(
                //System.out.printf(
                    "      i=%d minTourCost=%.1f\n", i, minTourCost));
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
    
}
