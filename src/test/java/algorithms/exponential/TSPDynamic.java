package algorithms.exponential;

import algorithms.util.FormatArray;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
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
 * The iterative method is limited to data holding 31 cities or less,
 * and the recursive method is limited to 2^31 -1 cities.
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
 * The runtime complexity is then O(n^2 * 2^n).
 * 
 * The 2nd stage is O(n) and does affect the runtime complexity.
 * 
 * 
 * NOTE that summation_over_k_from_1_to_{n-2}( C(n-3, k-1) )
 * </pre>
 * @author nichole
 */
public class TSPDynamic {

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
     * NOTE that algorithm can only handle N < 32.
     */
    public void solveIteratively() {
        
        if (N > 31) {
            throw new IllegalArgumentException("the dynamic iterative algorithm "
                + "has a limit of 31 for the number of cities.   Consider the "
                + "recursive algorithm which has a limit of 2^31 - 1");
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

        // r is the number of vertexes within n vertexes, in which the subset bits are set to 1.
        for (int r = 3; r <= N; r++) {
            // The number of ways to select an ordered sequence
            // of k items out of n distinct items for a fixed length of k
            // (a.k.a.  k-permutations of n)
            // = n! / (n − k)!
            //https://en.m.wikipedia.org/wiki/Permutation#k-permutations_of_n
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
    
    //===================================
    // https://www.interviewbit.com/blog/travelling-salesman-problem/
    //   with bug fixes here to the method least, and corrections to use
    // the start node.  sent the changes to interviewbit too.
    //NOTE: if implementing in C code and do not have a library with a hash set,
    //   it would be easier to replace uncompleted set with a doubly linked list
    //   which is easy to make, and then one must use the node back as an argument
    //   so the removal from the uncompleted linked list is O(1).
    // SEE notes below the "least" method for C programming implementation of DoubleLinkedList
    Set<Integer> uncompleted = null;
    private final int sentinel = Integer.MAX_VALUE;

    public void solveRecursively() {
        // init
        tour.clear();
        minTourCost = 0;
        solverFinished = false;
        uncompleted = new HashSet<Integer>();
        for (int i = 0; i < N; ++i) {
            uncompleted.add(i);
        }
        
        mincost(start);

        solverFinished = true;

        logger.log(LEVEL, String.format(
        //System.out.printf(
            "\ntour cost=%.0f\n", minTourCost));
    }
    /*
    mincost(start)
        mincost(  s1 = min_Arg_i (dist[start][i] + dist[c][i]) )
        mincost( s2 = ...)
    */
    
    // 2T(N-1)+1  is Θ(2^n) 
    private void mincost(int city) {
        
        //O(1)
        uncompleted.remove(city);

        logger.log(LEVEL, String.format(
        //System.out.printf(
            "%d--->", city));
        tour.add(city);
        
        // O(N-nCompleted)...   T(N-1)
        //minTourCost += dist[city][ncity] where ncity is the minimum i in dist[start][i] + dist[city][i]
        int ncity = least(city);

        if (ncity == sentinel) {
            logger.log(LEVEL, String.format(
            //System.out.printf(
                "%d", start));
            minTourCost += distance[city][start];
            tour.add(start);

            return;
        }

        mincost(ncity);
    }

    private int least(int c) {
        int nc = sentinel;
        double min = sentinel;
        double kmin = sentinel;
        
        // using i=0:N is O(N), so changed to loop only over "uncompleted" items
        //for (i = 0; i < N; i++) {
        for (int i : uncompleted) {
            if (distance[c][i] != 0) {
                if (distance[start][i] + distance[c][i] < min) {
                    min = distance[start][i] + distance[c][i];
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
    
    /*
    For C programming:
    ------------------

    the UML for the Node and DoubleLinkedList classes:

        DoubleLinkedList
    =========================
       head:Node
       tail:Node
       n:int
      ------------------------
       +insert(Node node):void
       +remove(Node node):void
       +search(int key)
       +first() : Node
       +last() : Node
       +isEmpty() : boolean
       +size() : int
      ------------------------

    Node
    =========================
       +key:int
       +next:Node
       +prev:Node
      ------------------------
      +Node(key:int)
      ------------------------

   Java code:
     public class DoubleLinkedList {
        protected Node head = null;
        protected Node tail = null;
        protected int n = 0;
        public DoubleLinkedList() {
        }
        // add node to end of double linked list
        public void insert(Node node) {
            if (head == null) {
                assert(tail == null);
                head = node;
                node.prev = null;
                tail = node;
                node.next = null;
                n++;
                return;
            }
            assert(tail != null);
            assert(tail.next == null);
            tail.next = node;
            node.prev = tail;
            node.next = null;
            n++;
        }

        public void remove(Node node) {
            // connect node.prev with node.next
            if (node.equals(head)) {
                if (n == 1) {
                    head = null;
                    assert(head.equals(tail));
                    tail = null;
                } else {
                    head = head.next;
                }
            } else if (node.equals(tail)) {
                assert(n > 1);
                tail.prev.next = null;
                tail = tail.prev;
          } else {
                assert(node.prev != null); // node is not head
                assert(node.next != null); // node is not tail
                Node p = node.prev;
                Node n2 = node.next;
                p.next = n2;
                n2.prev = p;
            }
            n--;
        }
        public Node first() {
            return head;
        }
        public Node last() {
            return tail;
        }
        public boolean isEmpty() {
            return (n == 0);
        }
        public int size() {
            return n;
        }
    }

    C Programming:
    // in C programming, to create an object that has itself as members, use
    // "forward declaration"

    // typedef'd struct are instantiated on the heap
    typedef struct Node Node;

    // consider declaring the node pointer (it's held in the function dynamic stack frame):
    typedef Node * NodePtr;

    struct Node {
      int key;
      Node *next;
      Node *prev;
    };
    
    OR
    
    struct Node {
      int key;
      NodePtr next;
      NodePtr prev;
    };

    // instantiate on the heap:
    NodePtr node = (NodePtr) malloc (sizeof(Node));

    // access members of node through ->
    node->next = node2;
}
*/
}
