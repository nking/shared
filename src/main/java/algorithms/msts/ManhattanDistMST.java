package algorithms.msts;

import algorithms.disjointSets.UnionFind;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Given a set of points, find the subset of "edges" that connect all points
 * resulting in a total sum of weights that is minimum.
 * The weights are the calculated Manhattan distance between the points.
 * The r.t.c. of this algorithm is O(n*log(n)).
 *
 * A detailed explanation of the lagorithm can be found in
 * https://cp-algorithms.com/geometry/manhattan-distance.html
 */
public class ManhattanDistMST {

    // implementing a Node instead of using a long[] so that the Set hash function finds
    // same data to be equivalent.
    protected static class Node {
        final long[] weightIJ;
        public Node(long weight, int i, int j) {
            weightIJ = new long[]{weight, i, j};
        }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof Node)) return false;
            return Arrays.equals(((Node)obj).weightIJ, weightIJ);
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(weightIJ);
        }
    }

    /**
     * calculate the manhattan based MST for the set of data points ps which are
     * rows of x, y coordinates as integers.
     * r.t.c. is O(n*log(n)).
     * @param ps
     * @return the MST as rows of {edge start index, edge end index, manhattan dist of edge}
     * where index is index of ps (hence represent a point).
     */
    public static List<long[]> manhattanMST(int[][] ps) {
        int n = ps.length;
        UnionFind uf = new UnionFind(n);
        // Node has weight, and indexes of coordinates
        List<Node> edges = new ArrayList<>(manhattanMSTEdges(ps));
        Collections.sort(edges, (o1, o2) -> {return Long.compare(o1.weightIJ[0], o2.weightIJ[0]);});
        List<long[]> out = new ArrayList<>();
        for (Node edge : edges) {
            int i = (int) edge.weightIJ[1];
            int j = (int) edge.weightIJ[2];
            if (uf.find(i) != uf.find(j)) {
                uf.union(i, j);
                out.add(new long[]{i, j, edge.weightIJ[0]});
            }
        }
        assert(out.size() == n-1);
        return out;
    }

    /**
     * given an array of points where each row is a 2-dimensional point, find
     * the disjoint edges.
     * r.t.c. is O(n*log(n)).
     <pre>
     reference:
     code was adapted from
     https://ebi-fly13.github.io/Library/graph/manhattan_mst.hpp.html
     The repository uses creative commons zero v.10 universal license https://github.com/ebi-fly13/Library/blob/main/LICENSE

     </pre>
     * @param ps
     * @return
     */
    protected static Set<Node> manhattanMSTEdges(int[][] ps) {

        int n = ps.length;

        List<Integer> ids = IntStream.range(0, n).boxed()
                .collect(Collectors.toCollection(ArrayList::new));

        // each element length is 3.  1st=, 2nd=, 3rd=
        Set<Node> edges = new HashSet<>();

        for (int s = 0; s < 2; s++) {
            for (int t = 0; t < 2; t++) {

                Collections.sort(ids, (i, j) -> {
                    return Integer.compare(
                            (ps[i][0] + ps[i][1]), (ps[j][0] + ps[j][1])
                    );
                });

                // using negative keys to make descending order
                TreeMap<Integer, Integer> sweep = new TreeMap<>();
                for (int i : ids) {

                    Set<Map.Entry<Integer, Integer>> rm = new HashSet<>();
                    for (Map.Entry<Integer, Integer> entry
                            : sweep.tailMap(-ps[i][1], true).entrySet()) {

                        rm.add(entry);

                        int j = entry.getValue();
                        if (ps[i][0] - ps[j][0] < ps[i][1] - ps[j][1]) break;

                        edges.add(new Node(
                                Math.abs(ps[i][1] - ps[j][1]) + Math.abs(ps[i][0] - ps[j][0]),
                                i, j));
                    }

                    for (Map.Entry<Integer, Integer> entry : rm) {
                        sweep.remove(entry.getKey(), entry.getValue());
                    }

                    sweep.put(-ps[i][1], i);
                }// end i from ids

                // swap all of x with all of y
                for (int i = 0; i < n; ++i) {
                    ps[i][0] ^= ps[i][1];
                    ps[i][1] ^= ps[i][0];
                    ps[i][0] ^= ps[i][1];
                }
            }// end t
            // mult x by -1
            for (int i = 0; i < n; ++i) {
                ps[i][0] *= -1;
            }
        }

        return edges;
    }
}
