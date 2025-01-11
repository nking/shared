package algorithms.trees;

import java.util.Arrays;

/**
 * a heap-like data structure for intervals.
 * It's update and query operations have r.t.c. O(log(n)) where n is the largest
 * index that the tree can hold.
 * pre-processing is done during construction and has r.t.c. O(n).
 * s.c. is O(n*log(n))
 *
 <pre>
 references:
 http://cp-algorithms.com/data_structures/segment_tree.html
 http://www.topcoder.com/thrive/articles/Range%20Minimum%20Query%20and%20Lowest%20Common%20Ancestor
 </pre>
 */
public class SegmentTree {

    // M holds indices
    protected final long[] tree;
    protected final int n;

    public SegmentTree(int[] a) {

        n = a.length;
        this.tree = new long[2*n];

        build(a, 1, 0, n-1);

        System.out.printf("tree=%s\n", Arrays.toString(tree));
    }

    /*
    with 1-based indexing:
        i = node
        i.left = 2*i
        i.right = 2*i + 1
        i.parent = i/2

    with 0-based indexing:
        i = node
        i.left = 2*i + 1
        i.right = 2*i + 2
        i.parent = (i-1)/2
     */

    protected void build(int[] a, int tIdx, int i, int j) {
        if (i == j) {
            tree[tIdx] = a[i];
            return;
        }
        int mid = (i+j)/2;
        int tLeftChild = 2*tIdx;
        int tRightChild = tLeftChild + 1;

        build(a,tLeftChild, i, mid);

        build(a,tRightChild, mid + 1, j);

        tree[tIdx] = tree[tLeftChild] + tree[tRightChild];
    }

    /**
     * return the sum from a[i0] to a[i1], inclusive.
     * @param i0
     * @param i1
     * @return
     */
    public long sum(int i0, int i1) {
        long s = sum(1, 0, n-1, i0, i1);
        return s;
    }

    protected long sum(int tIdx, int treeL, int treeR, int qL, int qR) {
        if (qL > qR) {
            return 0;
        }

        if (treeL == qL && treeR == qR) {
            return tree[tIdx];
        }

        int mid = (treeL + treeR)/2;
        int tLeftChild = 2*tIdx;
        int tRightChild = tLeftChild + 1;

        return sum(tLeftChild, treeL, mid, qL, Math.min(qR, mid))
                + sum(tRightChild, mid + 1, treeR, Math.max(qL, mid+1), qR);
    }
/*
    public void updateSet(int i, int newVal) {
        updateSet(1, 0, n-1, i, newVal);
    }

    protected void updateSet(int idx, int treeL, int treeR, int pos, int newVal) {
        if (treeL == treeR) {
            tree[idx] = newVal;
            return;
        }

        int mid = (treeL + treeR)/2;
        int iLeftChild = 2*idx;
        int iRightChild = iLeftChild + 1;

        if (pos <= mid) {
            updateSet(iLeftChild, treeL, mid, pos, newVal);
        } else {
            updateSet(iRightChild, mid+1, treeR, pos, newVal);
        }

        tree[idx] = tree[iLeftChild] + tree[iRightChild];
    }*/
}
