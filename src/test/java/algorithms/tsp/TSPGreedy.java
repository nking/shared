package algorithms.tsp;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * This is a greedy solution for TSP and is not guaranteed to give the 
 * optimal solution.
 * 
 * This is a corrected version of the InterviewBit code referenced below.
 * 
 * The approximate TSP algorithms are in the curvature scale space project.
 * 
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
 * Edits have been made to their code below.
 * 
 * also see
 * Held–Karp algorithm a.k.a. Bellman–Held–Karp algorithm
 * https://en.wikipedia.org/wiki/Held%E2%80%93Karp_algorithm
 * pp 58 - 
 * 
 * The runtime complexity is O(N^2).
 * 
 * </pre>
 * @author nichole
 */
public class TSPGreedy {

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
    public TSPGreedy(double[][] distance) {
        this(0, distance);
    }
    
    public TSPGreedy(int start, double[][] distance) {
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
            solveRecursively();
        }
        return tour;
    }
    // Returns the minimal tour cost.

    public double getTourCost() {
        if (!solverFinished) {
            solveRecursively();
        }
        return minTourCost;
    }

    /**
     * use a greedy iterative approach to find a TSP solution that is not
     * guaranteed to be optimal.
     * The runtime complexity is O(N^2).
     */
    public void solveIteratively() {
        // init
        tour.clear();
        minTourCost = 0;
        solverFinished = false;
        uncompleted = new HashSet<Integer>();
        for (int i = 0; i < N; ++i) {
            uncompleted.add(i);
        }
        
        prev = start;
        int city = start;
        
        int ncity;
        
        //N*(N-1)/2.
        while (true) {
        
            //O(1)
            uncompleted.remove(city);

            logger.log(LEVEL, String.format(
            //System.out.printf(
                "%d--->", city));
            tour.add(city);

            // O(N-nCompleted)... 
            ncity = least(city);

            if (ncity == sentinel) {
                logger.log(LEVEL, String.format(
                //System.out.printf(
                    "%d", start));
                minTourCost += distance[city][start];
                tour.add(start);

                break;
            }

            city = ncity;
        }
        
        solverFinished = true;

        logger.log(LEVEL, String.format(
        //System.out.printf(
            "%ntour cost=%.0f%n", minTourCost));
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

    /**
     * use a greedy recursive approach to find a TSP solution that is not
     * guaranteed to be optimal.
     * The runtime complexity is Θ(n^2) from a recurrence relation that is T(n–1) + n;
     */
    public void solveRecursively() {
        // init
        tour.clear();
        minTourCost = 0;
        solverFinished = false;
        prev = start;
        
        uncompleted = new HashSet<Integer>();
        for (int i = 0; i < N; ++i) {
            uncompleted.add(i);
        }
        //T(n–1) + n  = Θ(n^2)
        mincost(start);

        solverFinished = true;

        logger.log(LEVEL, String.format(
        //System.out.printf(
            "%ntour cost=%.0f%n", minTourCost));
    }
    /*
    mincost(start)
        mincost(  s1 = min_Arg_i (dist[start][i] + dist[c][i]) )
        mincost( s2 = ...)
    */
    
    private int prev;
        
    private void mincost(int city) {
        
        uncompleted.remove(city);

        logger.log(LEVEL, String.format(
        //System.out.printf(
            "%d--->", city));
        tour.add(city);
                
        // O(N-nCompleted)
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
        
        for (int i : uncompleted) {
            if (distance[c][i] != 0) {
                if (distance[prev][c] + distance[c][i] < min) {
                    min = distance[prev][c] + distance[c][i];
                    kmin = distance[c][i];
                    nc = i;
                    //System.out.printf("    nc=%d%n", nc);
                }
            }
        }

        if (min != sentinel) {
            minTourCost += kmin;
        }
        
        prev = c;

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

*/
}
