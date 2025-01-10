package algorithms.range;

import algorithms.trees.LeastCommonAncestor;

public class RMQ {

    protected final LeastCommonAncestor lca;

    /**
     * given array of values in a, build internal datastructure to
     * handle minimum within index range queries.
     * The r.t.c. of constructor is O(n) where n = a.length.
     * @param a
     */
    public RMQ(int[] a) {
        this.lca = new LeastCommonAncestor(a);
    }

    /**
     * given an array of queries, where a query is
     * an index range, return an array of the indices of the
     * minimum values in a for the queries.
     * The  r.t.c. is O(queries.length).
     * @param queries an array of queries in format queries[0] = [index0_0, index0_1],
     *                queries[1] = [index1_0, index1_1], etc.
     * @return
     */
    public int[] min(int[][] queries) {
        int[] out = new int[queries.length];
        for (int i = 0; i < queries.length; ++i) {
            out[i] = lca.find(queries[i][0], queries[i][1]);
        }
        return out;
    }

    /**
     * given a query array holding 2 indices as
     * an index range of array 'a', return the index of the
     * minimum value in array 'a' for the query index range.
     * The  r.t.c. is O(1).
     * @param query query in format  [index0, index1]
     * @return
     */
    public int min(int[] query) {
        return lca.find(query[0], query[1]);
    }
}
