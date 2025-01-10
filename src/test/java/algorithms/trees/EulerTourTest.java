package algorithms.trees;

import algorithms.graphs.DFS;
import algorithms.graphs.GraphUtil;
import algorithms.graphs.HierholzersEulerCircuit;
import junit.framework.TestCase;

import java.util.*;

public class EulerTourTest extends TestCase {

    /*
                           7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12

    pre-order visits: 7,3,1,0,2,5,4,6,11,9,8,10,13,12
       (root, left subtree, right subtree)
    in-order visits: 0,1,2,3,4,5,6,7,8,9,10,11,12,13
      (left subtree, root, right subtree)
    post-order visits: 0,2,1,4,6,5,3, 8,10,9,12,13,11,7
      (left subtree, right subtree, root)
     */
    public void test0() {

        int[] expected = new int[]{
                7, 3, 1, 0, 2, 5, 4, 6, 11, 9, 8, 10, 13, 12
        };

        Map<Integer, Set<Integer>> adjMap = buildTreeAsGraph();

        // euler circuit starting from root=7, is same as pre-order tree traversal.
        int src = 7;
        HierholzersEulerCircuit euler = new HierholzersEulerCircuit();
        int[] eulerCircuit = euler.createCircuit(GraphUtil.convertGraph2(adjMap), src);

        DFS dfs = new DFS(GraphUtil.convertGraph(adjMap));
        dfs.walk(src);

        int[] eulerTourDFS = dfs.getEulerTourFromEndTimes();

        BinaryTreeNode<Integer> tree = buildTree();
        TreeTraversal traversals = new TreeTraversal();
        int[] eulerTourPreOrder = traversals.getPreorder(tree);

        assertTrue(Arrays.equals(expected, eulerCircuit));
        assertTrue(Arrays.equals(expected, eulerTourDFS));
        assertTrue(Arrays.equals(expected, eulerTourPreOrder));
    }
    /*
                           7
                3                        11
           1         5              9         13
        0   2     4   6           8   10    12
     */
    protected BinaryTreeNode<Integer> buildTree() {
        BinaryTreeNode<Integer> tree = new BinaryTreeNode<>(7);
        tree.setRight(new BinaryTreeNode<>(11));
        tree.getRight().setRight(new BinaryTreeNode<>(13));
        tree.getRight().getRight().setLeft(new BinaryTreeNode<>(12));

        tree.getRight().setLeft(new BinaryTreeNode<>(9));
        tree.getRight().getLeft().setLeft(new BinaryTreeNode<>(8));
        tree.getRight().getLeft().setRight(new BinaryTreeNode<>(10));

        tree.setLeft(new BinaryTreeNode<>(3));

        tree.getLeft().setLeft(new BinaryTreeNode<>(1));
        tree.getLeft().getLeft().setLeft(new BinaryTreeNode<>(0));
        tree.getLeft().getLeft().setRight(new BinaryTreeNode<>(2));

        tree.getLeft().setRight(new BinaryTreeNode<>(5));
        tree.getLeft().getRight().setLeft(new BinaryTreeNode<>(4));
        tree.getLeft().getRight().setRight(new BinaryTreeNode<>(6));

        return tree;
    }
    protected Map<Integer, Set<Integer>> buildTreeAsGraph() {
        Map<Integer, Set<Integer>> adjMap = new HashMap<>();

        adjMap.put(7, new HashSet<Integer>());
        adjMap.get(7).add(3);  adjMap.get(7).add(11);

        adjMap.put(11, new HashSet<Integer>());
        adjMap.get(11).add(9);  adjMap.get(11).add(13);

        adjMap.put(9, new HashSet<Integer>());
        adjMap.get(9).add(8);  adjMap.get(9).add(10);

        adjMap.put(13, new HashSet<Integer>());
        adjMap.get(13).add(12);

        adjMap.put(3, new HashSet<Integer>());
        adjMap.get(3).add(1);  adjMap.get(3).add(5);

        adjMap.put(1, new HashSet<Integer>());
        adjMap.get(1).add(0);  adjMap.get(1).add(2);

        adjMap.put(5, new HashSet<Integer>());
        adjMap.get(5).add(4);  adjMap.get(5).add(6);

        return adjMap;
    }

    protected Map<Integer, Set<Integer>> createTree2() {

        Map<Integer, Set<Integer>> adjMap = new HashMap<>();

        int u=1, v=2;
        adjMap.putIfAbsent(u, new HashSet<>());
        adjMap.putIfAbsent(v, new HashSet<>());
        adjMap.get(u).add(v);
        adjMap.get(v).add(u);

        u=2;
        v=5;
        adjMap.putIfAbsent(u, new HashSet<>());
        adjMap.putIfAbsent(v, new HashSet<>());
        adjMap.get(u).add(v);
        adjMap.get(v).add(u);

        u=2;
        v=6;
        adjMap.putIfAbsent(u, new HashSet<>());
        adjMap.putIfAbsent(v, new HashSet<>());
        adjMap.get(u).add(v);
        adjMap.get(v).add(u);

        u=1;
        v=3;
        adjMap.putIfAbsent(u, new HashSet<>());
        adjMap.putIfAbsent(v, new HashSet<>());
        adjMap.get(u).add(v);
        adjMap.get(v).add(u);

        u=1;
        v=4;
        adjMap.putIfAbsent(u, new HashSet<>());
        adjMap.putIfAbsent(v, new HashSet<>());
        adjMap.get(u).add(v);
        adjMap.get(v).add(u);

        u=4;
        v=7;
        adjMap.putIfAbsent(u, new HashSet<>());
        adjMap.putIfAbsent(v, new HashSet<>());
        adjMap.get(u).add(v);
        adjMap.get(v).add(u);

        return adjMap;

    }
}
