package algorithms.sort;

import algorithms.misc.MiscMath0;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

public class BucketSort {

    /**
     * sort a, and for ties, prefer the element with a smaller index on a.
     * runtime complexity is on average O(n) where n is the length of a, but
     * if the distribution of values in a has a very large range, and large
     * clusters distibuted over small number of integers, the runtime complexity
     * will be larger than O(n), but still less than O(n*log_2(n)).
     * @param a
     * @return the sorted indexes of a.  these can be used to sort other arrays of data that are associated with a.
     */
    public static int[] sortAndReturnIndexes(final double[] a) {

        int n = a.length;
        double[] minMax = MiscMath0.getMinMax(a);
        double binWidth = (minMax[1] - minMax[0])/(n - 1.);

        int[] indexes = new int[n];

        int i;

        List<LinkedList<Node>> b = new ArrayList<>();
        for (i = 0; i < n; ++i) {
            b.add(new LinkedList<Node>());
        }
        int iB;
        for (i = 0; i < n; ++i) {
            iB = (int) ((a[i] - minMax[0])/binWidth);
            b.get(iB).add(new Node(i, a[i]));
        }
        int iOut = 0;
        double[] out = new double[n];
        LinkedList<Node> d;
        Node node;
        for (i = 0; i < n; ++i) {
            d = b.get(i);
            if (d.isEmpty()){
                continue;
            }
            Collections.sort(d);
            while (!d.isEmpty()) {
                node = d.removeFirst();
                out[iOut] = node.val;
                indexes[iOut] = node.idx;
                ++iOut;
            }
        }
        System.arraycopy(out, 0, a, 0, n);
        return indexes;
    }

    /**
     * sort a, and for ties, prefer the element with a smaller index on a.
     * runtime complexity is on average O(n) where n is the length of a, but
     * if the distribution of values in a has a very large range, and large
     * clusters distibuted over small number of integers, the runtime complexity
     * will be larger than O(n), but still less than O(n*log_2(n)).
     * @param a
     * @return the sorted indexes of a.  these can be used to sort other arrays of data that are associated with a.
     */
    public static int[] sortAndReturnIndexes(final int[] a) {

        int n = a.length;
        int[] minMax = new int[]{MiscMath0.findMin(a), MiscMath0.findMax(a)};
        double binWidth = ((double)(minMax[1] - minMax[0]))/n;

        int[] indexes = new int[n];

        int i;

        List<LinkedList<Node2>> b = new ArrayList<>();
        for (i = 0; i < n+1; ++i) {
            b.add(new LinkedList<Node2>());
        }
        int iB;
        for (i = 0; i < n; ++i) {
            iB = (int) (((double)(a[i] - minMax[0]))/binWidth);
            b.get(iB).add(new Node2(i, a[i]));
        }
        int iOut = 0;
        int[] out = new int[n];
        LinkedList<Node2> d;
        Node2 node;
        for (i = 0; i < b.size(); ++i) {
            d = b.get(i);
            if (d.isEmpty()){
                continue;
            }
            Collections.sort(d);
            while (!d.isEmpty()) {
                node = d.removeFirst();
                out[iOut] = node.val;
                indexes[iOut] = node.idx;
                ++iOut;
            }
        }
        System.arraycopy(out, 0, a, 0, n);
        return indexes;
    }

    private static class Node implements Comparable<Node> {
        public final int idx;
        public final double val;
        public Node(int i, double v) {
            this.idx = i;
            this.val = v;
        }

        @Override
        public int compareTo(Node o) {
            if (this.val < o.val) {
                return -1;
            } else if (this.val == o.val) {
                return Integer.compare(this.idx, o.idx);
            }
            return +1;
        }
    }

    private static class Node2 implements Comparable<Node2> {
        public final int idx;
        public final int val;
        public Node2(int i, int v) {
            this.idx = i;
            this.val = v;
        }

        @Override
        public int compareTo(Node2 o) {
            if (this.val < o.val) {
                return -1;
            } else if (this.val == o.val) {
                return Integer.compare(this.idx, o.idx);
            }
            return +1;
        }
    }
}
