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
    /*
    given an array of numbers a,
    stores the sum within index ranges in a tree of non-overlapping ranges that
    proceeed from root being index range [i=0, j=n-1].
    left being [i, (i+j)/2] and right being [((i+j)/2)+1, j].
    The leaf nodes are single indices.

    The tree can be built for sums or other operations,
    though BigInteger is not used internally, so the input should be such that operations
     won't overflow for a long (values < (1<<63)-1).
     */
    public static enum TYPE {
        SUM
        //can add more types such as mult, max, freq, etc
    }

    // M holds indices
    protected final long[] tree;
    protected final int n;

    protected final TYPE type;

    protected final boolean use0BasedIndexes;

    public SegmentTree(int[] a, TYPE type, boolean use0BasedIndexes) {
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }
        n = a.length;
        this.tree = new long[4*n];
        this.type = type;
        this.use0BasedIndexes = use0BasedIndexes;

        int idx = use0BasedIndexes ? 0 : 1;

        build(a, idx, 0, n-1);
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

    protected void build(int[] a, int idx, int i, int j) {
        if (i == j) {
            tree[idx] = a[i];
            return;
        }
        int mid = (i+j)/2;
        int iLeftChild = 2*idx;
        int iRightChild = iLeftChild + 1;
        if (use0BasedIndexes) {
            ++iLeftChild;
            ++iRightChild;
        }
        build(a,iLeftChild, i, mid);
        build(a,iRightChild, mid + 1, j);
        tree[idx] = tree[iLeftChild] + tree[iRightChild];
    }

    public long sum(int i0, int i1) {
        int idx = use0BasedIndexes ? 0 : 1;
        long s = sum(idx, 0, n-1, i0, i1);
        return s;
    }

    protected long sum(int idx, int treeL, int treeR, int qL, int qR) {
        if (qL > qR)
            return 0;

        if (treeL == qL && treeR == qR) {
            return tree[idx];
        }

        int mid = (treeL + treeR)/2;
        int iLeftChild = 2*idx;
        int iRightChild = iLeftChild + 1;
        if (use0BasedIndexes) {
            ++iLeftChild;
            ++iRightChild;
        }

        //long s1 = sum(iLeftChild, treeL, mid, qL, Math.min(qR, mid));
        //long s2 = sum(iRightChild, mid + 1, treeR, Math.max(qL, mid+1), qR);
        //return s1 + s2;
        return sum(iLeftChild, treeL, mid, qL, Math.min(qR, mid))
                + sum(iRightChild, mid + 1, treeR, Math.max(qL, mid+1), qR);
    }

    public void update(int idx, int treeL, int treeR, int pos, int newVal) {
        if (treeL == treeR) {
            tree[idx] = newVal;
            return;
        }

        int mid = (treeL + treeR)/2;
        int iLeftChild = 2*idx;
        int iRightChild = iLeftChild + 1;
        if (use0BasedIndexes) {
            ++iLeftChild;
            ++iRightChild;
        }

        if (pos <= mid) {
            update(iLeftChild, treeL, mid, pos, newVal);
        } else {
            update(iRightChild, mid+1, treeR, pos, newVal);
        }

        tree[idx] = tree[iLeftChild] + tree[iRightChild];
    }
}
