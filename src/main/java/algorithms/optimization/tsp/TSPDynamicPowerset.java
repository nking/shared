package algorithms.optimization.tsp;

import java.util.Arrays;

/**
 * Given a list of distances between pairs of cities, the travelling salesperson problem finds the shortest
 * path that visits each city exactly once.
 * The path is a Hamiltonian cycle of minimum weight.
 * This class uses a dynamic programming approach
 * of a bitstring to hold the 2^n possible states of visiting cities while re-using
 * smaller solutions for bigger ones.
 */
public class TSPDynamicPowerset {

    /**
     * find the minimum sum Hamiltonian cycle for the nodes given by the dist matrix.
     * Note that node 0 has been arbitrarily adopted as the start vertex, but the
     * code could easily be adapted for another src node.
     * r.t.c. O((n^2)*(2^n)).  s.c. is O(n*(2^n)).
     * @param dist
     * @return an array of the minimum cycle nodes in path order.
     */
    public static int[] findMinCycle(double[][] dist) {
        int n = dist.length;
        if (dist[0].length != n) {
            throw new IllegalArgumentException("expected a square matrix");
        }
        if (n > 12) {
            throw new IllegalArgumentException("for dynamic programming approach and a COT computer, " +
                    "n should be <= 12");
        }

        /*
        first dimension holds the bits set, that is the set bits of the path so far.
        second dimension holds the last node set in the path.  need for book-keeping to know that next
        node added is an edge starting with this last node.
         */
        double[][] tab = new double[1<<n][n];

        //NOTE: can consider more efficient book keeping because n is small, but will use string for readability
        String[][] tabPath = new String[1<<n][n];

        //NOTE: have edited the below so that we always start at vertex 0 since it is a cycle,
        // and the relative order of an answer will be the same going around the cycle.

        /*
        considering subproblems.

        simplest approach is start w/ s and walk thru what other vars and dimensions are needed to hold state:
           init
              fill all of tab with sentinel = Integer.MAX_VALUE.
              tab[1<<i] = 0 for i=0 to n

              for s = 0 to 1<<n
                 we have a state s
                 for i=
                    subtract a node i to get prev state
                    tab[s] = max ( tab[s], tab[prev] + edge from last node in s to i

                    so now we see we need to store last node in tab too
                    and iterate over each j as last j where it makes sense

         */
        double sentinel = Double.POSITIVE_INFINITY;
        for (int i = 0; i < (1<<n); ++i) {
            Arrays.fill(tab[i], sentinel);
        }
        /* // this solves for any node as the start vertex.  need to adapth i loop below for it if enable.
        for (int i = 0; i < n; ++i) {
            // init with 0 for all single node states.  the end node is itself
            tab[1<<i][i] = 0;
            tabPath[1<<i][i] = Integer.toString(i);
        }*/
        // init for starting at i=0
        tab[1<<0][0] = 0;
        tabPath[1<<0][0] = Integer.toString(0);

        for (int s = 1; s < (1<<n); ++s) {

            // i is the potential node to add to the state prev = S - i, that is S \ i
            // we skip i=0 since we already set that as the starting node
            for (int i = 1; i < n; ++i) {
                // make sure s and i have at least 1 bit in common
                if ((s & (1<<i)) == 0) continue;

                int sPrev = (s ^ (1<<i));

                // loop over all set bits in sPrev as j
                int tmp = sPrev;

                while (tmp > 0) {
                    int j = (int)(Math.log(tmp & -tmp)/Math.log(2)); //lsb = tmp & -tmp

                    if (tab[sPrev][j] != sentinel) {
                        double d = tab[sPrev][j] + dist[j][i];

                        if (tab[s][i] > d) {
                            tab[s][i] = d;
                            tabPath[s][i] = String.format("%s,%d", tabPath[sPrev][j], i);
                        }

                    }
                    tmp &= (tmp - 1);
                }
            }
        }

        // find min of all states covering all nodes.
        // need to add the dist cost from last to first point to make it a cycle
        int s = (1<<n)-1;
        double min = Integer.MAX_VALUE;
        String[] minPath = null;
        int minIdx = 0;
        for (int i = 0; i < n; ++i) {
            if (tab[s][i] == sentinel) continue;
            String[] path = tabPath[s][i].split(",");
            double add = dist[Integer.parseInt(path[path.length-1])][Integer.parseInt(path[0])];
            if (min > (tab[s][i] + add)) {
                min = tab[s][i] + add;
                minIdx = i;
                minPath = path;
            }
        }
        int[] path = new int[n+1];
        for (int i = 0; i < n; ++i) {
            path[i] = Integer.parseInt(minPath[i]);
        }
        path[n] = path[0];
        return path;
    }
}
