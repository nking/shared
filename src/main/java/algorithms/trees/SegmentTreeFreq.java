package algorithms.trees;

import java.util.*;

public class SegmentTreeFreq {

    protected static class Node {
        Map<Integer, Integer> countMap = new HashMap<>();
        public Node(int v) {
            countMap.put(v, 1);
        }
        public Node(Node nodeL, Node nodeR) {
            countMap.putAll(nodeL.countMap);
            for (Map.Entry<Integer, Integer> entry : nodeR.countMap.entrySet()) {
                countMap.put(entry.getKey(), countMap.getOrDefault(entry.getKey(), 0) +
                        entry.getValue());
            }
        }
    }

    protected final Node[] tree;
    protected final int n;

    /**
     * constuctor for segment tree.  r.t.c. O(n*log(n))
     * @param a array of numbers that frequency queries will be performed on
     */
    public SegmentTreeFreq(int[] a) {
        this.n = a.length;
        this.tree = new Node[2*n];
        build(a, 1, 0, n-1);
    }

    protected void build(int[] a, int tIdx, int treeL, int treeR) {
        if (treeL == treeR) {
            tree[tIdx] = new Node(a[treeL]);
            return;
        }
        int mid = (treeL+treeR)/2;
        //int tLeftChild = 2*idx;
        //int tRightChild = iLeftChild + 1;
        int tLeftChild = tIdx + 1;
        int tRightChild = tIdx + (2 * (mid - treeL + 1));
        // this is tRightChild = tIdx + treeR - treeL + 2

        build(a,tLeftChild, treeL, mid);

        build(a,tRightChild, mid + 1, treeR);

        tree[tIdx] = new Node(tree[tLeftChild], tree[tRightChild]);
    }

    /**
     * query for frequencies of integer counts for the given ranges.
     * r.t.c. is O(Q*log(N)) where Q = q.length and N = a.length from array given to constructor.
     * (note that construction r.t.c. is O(N*log(N))) and so total is
     * O((Q+N)log(N)).
     * In contrast, MosAlgorithm.queryFrequencies which is O((Q+N)sqrt(N)).
     * @param queries an array of query ranges.  each row is range [i1, i2] inclusive
     *                where i1 and i2 are indexes of array a.  the indexes should be
     *                0-based.
     * @return array of frequency maps, 1 for each row of queries
     */
    public List<Map<Integer, Integer>> query(int[][] queries) {
        List<Map<Integer, Integer>> out = new ArrayList<>();
        for (int[] q : queries) {
            out.add(query(q[0], q[1]));
        }
        return out;
    }

    public Map<Integer, Integer> query(int i0, int i1) {
        return query(1, 0, n-1, i0, i1);
    }

    protected Map<Integer, Integer> query(int tIdx, int treeL, int treeR, int qL, int qR) {
        if (tIdx < 1) {
            throw new IllegalArgumentException("tIdx must be >= 1");
        }
        if (qL > qR) {
            return new HashMap<>();
        }

        if (treeL == qL && treeR == qR) {
            return tree[tIdx].countMap;
        }

        int mid = (treeL + treeR)/2;
        //int tLeftChild = 2*tIdx;
        //int tRightChild = tLeftChild + 1;
        int tLeftChild = tIdx + 1;
        int tRightChild = tIdx + (2 * (mid - treeL + 1));
        // this is tRightChild = tIdx + treeR - treeL + 2

        Map<Integer, Integer> map = new HashMap<>();
        Map<Integer, Integer> leftMap = query(tLeftChild, treeL, mid, qL, Math.min(qR, mid));
        Map<Integer, Integer> rightMap = query(tRightChild, mid + 1, treeR, Math.max(qL, mid+1), qR);
        map.putAll(leftMap);
        for (Map.Entry<Integer, Integer> entry : rightMap.entrySet()) {
            map.put(entry.getKey(), map.getOrDefault(entry.getKey(), 0) +
                    entry.getValue());
        }

        return map;
    }
}