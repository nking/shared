package algorithms.trees;

import algorithms.graphs.HierholzersEulerCircuit;

import javax.print.attribute.standard.MediaSize;
import java.util.*;

/**
 * The Least Common Ancestor (LCA) problem (a.k.a. lowest common ancestor):
 * given a tree, find the Least Common Ancestor (LCA) of a pair of nodes.
 * The LCA of nodes and in a tree is the shared ancestor of and that is located farthest
 * from the root.
 * <p>
 * The goal is to pre-process the data so that subsequent queries are fast.
 * <p>
 * This has been written for 0-based indexing, that is, the query indexes start at 0 and
 * are w.r.t. the original array 'a' given to the constructor.
 *
 * <pre>
 * reference:
 * "The LCA Problem Revisited"
 * Michael A. Bender, Martın Farach-Colton
 * SUNY Stony Brook, Rutgers University
 * May 16, 2000
 * </pre>
 */
public class LeastCommonAncestor {

    /*
    LCA:
        given:  a rooted tree T having n nodes.

        query:
           For nodes u and v and of tree T, query LCA returns the least common ancestor
           of u and v in T, that is, it returns the node furthest from the root that
           is an ancestor of both and u and v.

        a dFS traversal of the tree for paths that encounter u and v finds the
        LCA(u,v) as the deepest node on both paths before reaching u and v.

    related:
    Range Minimum Query (RMQ) Problem,
        given:  A length array of numbers

        query: for indices i and j between 1 and n, query RMQ_A(x,y) returns the index
               of the smallest element in the subarray A[i...j]

        pre-processing time: f(n)
        query time: q(n)
        r.t.c. stated as <f(n), q(n)>

     */

    protected final int[] tree;
    protected final int treeSrc;
    protected final int[][] eulerCircuitAndDepth;
    protected int[] firstVisit;
    protected int n;

    protected final int[] log2;
    protected final int blockSize;
    protected final int blockCnt;
    protected final int[][] st;
    protected final int[][][] blocks;
    protected final int[] blockMask;

    /**
     * constructor r.t.c. is O(n) where n = a.length.
     *
     * @param a
     */
    public LeastCommonAncestor(int[] a) {

        n = a.length;

        // indices of parent nodes.  e.g. tree[v]=u says a[v] is a child of a[u]
        this.tree = makeCartesianTree(a);
        assert (tree.length == n);

        int src = -1;
        for (int i = 0; i < n; ++i) {
            if (tree[i] == -1) {
                src = i;
                break;
            }
        }
        if (src == -1) {
            throw new IllegalArgumentException("Error in algorithm. see tree=" + Arrays.toString(tree));
        }
        this.treeSrc = src;


        // create a bi-directional adjacency map out of the indices parent child relationships:
        Map<Integer, LinkedList<Integer>> adjMap = createBiDirectionalMap(tree);

        // form a euler circuit (traversing edges once to connect all nodes and return to start).
        HierholzersEulerCircuit euler = new HierholzersEulerCircuit();
        // row 0 : the indices of array 'a' in order of the euler circuit.
        // row 1 : the depth of the node in tree.  depth of src = 0.
        // the difference between adjacent elements in eulerCircuitAndDepth[1] is exactly 1
        this.eulerCircuitAndDepth = euler.createCircuitAndDepth(adjMap, src);

        assert (eulerCircuitAndDepth[0].length == (2 * n - 1));

        //Store occurence of first occurrence each vertex in E another array firstVisit.
        // Size will be n as there are n vertices.
        this.firstVisit = new int[n];
        Arrays.fill(firstVisit, -1);
        int aIdx;
        for (int i = 0; i < eulerCircuitAndDepth[0].length; ++i) {
            aIdx = eulerCircuitAndDepth[0][i];
            if (firstVisit[aIdx] == -1) {
                firstVisit[aIdx] = i;
            }
        }

        //https://cp-algorithms.com/graph/lca_farachcoltonbender.html

        // precompute all log values
        int m = eulerCircuitAndDepth[0].length;
        this.log2 = new int[m + 1];
        log2[0] = -1;
        for (int i = 1; i <= m; i++) {
            log2[i] = (log2[i/2] + 1);
        }

        this.blockSize = Math.max(1, log2[m]/2);
        // round-up m/blocksize:
        this.blockCnt = (m + blockSize - 1) / blockSize;

        // precompute minimum of each block and build sparse table
        //st.assign(block_cnt, vector<int>(log_2[block_cnt] + 1));
        this.st = new int[blockCnt][log2[blockCnt] + 1];
        for (int i = 0, j = 0, b = 0; i < m; i++, j++) {
            if (j == blockSize) {
                j = 0;
                b++;
            }
            if (j == 0 || minByH(i, st[b][0]) == i) {
                st[b][0] = i;
            }
        }
        for (int l = 1; l <= log2[blockCnt]; l++) {
            for (int i = 0; i < blockCnt; i++) {
                int ni = i + (1 << (l - 1));
                if (ni >= blockCnt) {
                    st[i][l] = st[i][l - 1];
                } else {
                    st[i][l] = minByH(st[i][l - 1], st[ni][l - 1]);
                }
            }
        }

        // precompute mask for each block
        //block_mask.assign(blockCnt, 0);
        this.blockMask = new int[blockCnt];
        for (int i = 0, j = 0, b = 0; i < m; i++, j++) {
            if (j == blockSize) {
                j = 0;
                b++;
            }
            if (j > 0 && (i >= m || minByH(i - 1, i) == i - 1)) {
                blockMask[b] += 1 << (j - 1);
            }
        }

        // precompute RMQ for each unique block
        int possibilities = 1 << (blockSize - 1);
        //blocks.resize(possibilities);
        this.blocks = new int[possibilities][][];

        for (int b = 0; b < blockCnt; b++) {
            int mask = blockMask[b];
            if (blocks[mask] != null) {
                continue;
            }
            //blocks[mask].assign(block_size, vector<int>(block_size));
            blocks[mask] = new int[blockSize][blockSize];
            for (int l = 0; l < blockSize; l++) {
                blocks[mask][l][l] = l;
                for (int r = l + 1; r < blockSize; r++) {
                    blocks[mask][l][r] = blocks[mask][l][r - 1];
                    if (b * blockSize + r < m) {
                        blocks[mask][l][r] = minByH(b * blockSize + blocks[mask][l][r],
                                b * blockSize + r) - b * blockSize;
                    }
                }
            }
        }
    }

    protected int minByH(int i, int j) {
        //return height[euler_tour[i]] < height[euler_tour[j]] ? i : j;
        return eulerCircuitAndDepth[1][i] < eulerCircuitAndDepth[1][j] ? i : j;
    }

    //O(1) method, not yet corrected for possible 1-based indexing
    public int find(int i0, int i1) {
        //https://cp-algorithms.com/graph/lca_farachcoltonbender.html

        int l = firstVisit[i0];
        int r = firstVisit[i1];
        if (l > r) { // swap
            l ^= r;
            r ^= l;
            l ^= r;
        }
        int bl = l / blockSize;
        int br = r / blockSize;
        if (bl == br) {
            int idx = lcaInBlock(bl, l % blockSize, r % blockSize);
            return eulerCircuitAndDepth[0][idx];
        }
        int ans1 = lcaInBlock(bl, l % blockSize, blockSize - 1);
        int ans2 = lcaInBlock(br, 0, r % blockSize);
        int ans = minByH(ans1, ans2);
        if (bl + 1 < br) {
            int l2 = log2[br - bl - 1];
            int ans3 = st[bl+1][l2];
            int ans4 = st[br - (1 << l2)][l2];
            ans = minByH(ans, minByH(ans3, ans4));
        }
        return eulerCircuitAndDepth[0][ans];
    }

    int lcaInBlock(int b, int l, int r) {
        return blocks[blockMask[b]][l][r] + b * blockSize;
    }

    /**
     * find the least common ancestor for indices i0 and i1 where the indices are w/ respect to the
     * original array give to constructor.
     * It finds the index in array 'a' for the minimum value in the query index range [i0, i1].
     * The r.t.c. is O(log_2(n)) where n is a.length.
     *
     * @param i0
     * @param i1
     * @return index into array 'a' of the least common ancestor of
     */
    public int findWithLogN(int i0, int i1) {

        Set<Integer> common = new HashSet<>();

        int iLeft = i0;
        int iRight = i1;
        while (iLeft != -1 || iRight != -1) {
            if (iLeft == iRight) {
                return iLeft;
            } else if (common.contains(iLeft)) {
                return iLeft;
            } else if (common.contains(iRight)) {
                return iRight;
            }
            if (iLeft != -1) {
                common.add(iLeft);
                iLeft = tree[iLeft];
            }
            if (iRight != -1) {
                common.add(iRight);
                iRight = tree[iRight];
            }
        }
        return -1;
    }

    /**
     * makes the parent array of a as a cartesion tree.
     *
     * @param a
     * @return
     */
    protected int[] makeCartesianTree(int[] a) {
        int n = a.length;
        int[] parent = new int[n];
        Arrays.fill(parent, -1);

        Stack<Integer> s = new Stack<>();

        // build cartesian tree using montonic decr queue
        for (int i = 0; i < n; ++i) {
            int last = -1;
            while (!s.isEmpty() && a[s.peek()] >= a[i]) {
                last = s.pop();
            }
            if (!s.isEmpty()) {
                parent[i] = s.peek();
            }
            if (last > -1) {
                parent[last] = i;
            }
            s.push(i);
        }
        return parent;
    }

    private Map<Integer, LinkedList<Integer>> createBiDirectionalMap(int[] tree) {
        Map<Integer, LinkedList<Integer>> adjMap = new HashMap<>();
        int u;
        for (int v = 0; v < tree.length; ++v) {
            u = tree[v];
            if (u == -1) {
                // root node
                continue;
            }
            adjMap.putIfAbsent(u, new LinkedList<>());
            adjMap.putIfAbsent(v, new LinkedList<>());
            adjMap.get(u).add(v);
            adjMap.get(v).add(u);
        }
        return adjMap;
    }

    /**
     * find the distance between 2 nodes in a tree where root is the tree, and node1Val and
     * node2Val are the nodes in the tree to find.
     * r.t.c. is O(n) because the tree is not necessarily a binary search tree.
     <pre>
     method is adapted from Chapter 18 of Competitive Programming Handbook by Antti Laaksonen, Chap 18.
     </pre>
     * @param root
     * @param node1Val
     * @param node2Val
     * @return the distance between the first nodes found in root tree with values node1Val and
     * node2Val.  Note that if either node1Val or node2VAl are not found, the return is Long.MAX_VALUE.
     */
    protected static long distBetweenNodes(NAryTreeNode root, int node1Val, int node2Val) {
        /*
        find the nodes in the tree and store the depths of nodes along the way.
        when both are found, ascend the tree for the node deeper than the other until
        both nodes are same depth.
        then while nodes are not the same node, ascend to find common parent.
        the total dist = depth(node1) + depth(node2) - 2*depth(lca)
         */
        Map<NAryTreeNode, Integer> depthMap = new HashMap<>();
        NAryTreeNode node1 = null;
        NAryTreeNode node2 = null;

        // use level order traversal until find both
        Queue<NAryTreeNode> q = new ArrayDeque<>();
        q.offer(root);
        depthMap.put(root, 1);
        if (root.getData() == node1Val) {
            node1 = root;
        }
        if (root.getData() == node2Val) {
            node2 = root;
        }
        while (!q.isEmpty() && (node1 == null || node2 == null)) {
            root = q.poll();
            int level = depthMap.get(root);
            for (NAryTreeNode ch : root.getChildren()) {
                depthMap.put(ch, level + 1);
                if (ch.getData() == node1Val) {
                    node1 = ch;
                }
                if (ch.getData() == node2Val) {
                    node2 = ch;
                }
                if (node1 != null && node2 != null) break;
                q.add(ch);
            }
        }
        if (node1 == null || node2 == null) {
            return Long.MAX_VALUE;
        }
        int d1 = depthMap.get(node1);
        int d2 = depthMap.get(node2);

        while (node1 != null && depthMap.get(node1) > depthMap.get(node2)) {
            node1 = node1.getParent();
        }
        while (node1 != null && node2 != null && depthMap.get(node2) > depthMap.get(node1)) {
            node2 = node2.getParent();
        }
        while (node1 != null && node2 != null && !node1.equals(node2)) {
            node1 = node1.getParent();
            node2 = node2.getParent();
        }
        int dLCA = depthMap.get(node1);

        return d1 + d2 - 2*dLCA;
    }

    /**
     * given tree 'root', find the distance between node1 and node2 in the tree.
     * @param root
     * @param node1
     * @param node2
     * @return the distance (number of edges) between node1 and node2.  If either node
     * is not found in the tree, Long.MAX_VALUE is returned.
     */
    protected static long distBetweenNodes(NAryTreeNode root, NAryTreeNode node1, NAryTreeNode node2) {
         /*
        find the nodes in the tree and store the depths of nodes along the way.
        when both are found, ascend the tree for the node deeper than the other until
        both nodes are same depth.
        then while nodes are not the same node, ascend to find common parent.
        the total dist = depth(node1) + depth(node2) - 2*depth(lca)
         */
        Map<NAryTreeNode, Integer> depthMap = new HashMap<>();

        // use level order traversal until find both
        Queue<NAryTreeNode> q = new ArrayDeque<>();
        q.offer(root);
        depthMap.put(root, 1);
        if (root.equals(node1)) {
            node1 = root;
        }
        if (root.equals(node2)) {
            node2 = root;
        }
        while (!q.isEmpty() && !depthMap.containsKey(node1) || !depthMap.containsKey(node2)) {
            root = q.poll();
            int level = depthMap.get(root);
            for (NAryTreeNode ch : root.getChildren()) {
                depthMap.put(ch, level + 1);
                if (depthMap.containsKey(node1) && depthMap.containsKey(node2)) break;
                q.add(ch);
            }
        }
        if (!depthMap.containsKey(node1) || !depthMap.containsKey(node2)) {
            return Long.MAX_VALUE;
        }

        int d1 = depthMap.get(node1);
        int d2 = depthMap.get(node2);
        //
        while (node1 != null && depthMap.get(node1) > depthMap.get(node2)) {
            node1 = node1.getParent();
        }
        while (node1 != null && node2 != null && depthMap.get(node2) > depthMap.get(node1)) {
            node2 = node2.getParent();
        }
        while (node1 != null && node2 != null && !node1.equals(node2)) {
            node1 = node1.getParent();
            node2 = node2.getParent();
        }
        int dLCA = depthMap.get(node1);

        return d1 + d2 - 2*dLCA;
    }

}
