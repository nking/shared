package algorithms.disjointSets;

import java.util.*;

/**
 * a version of Tarjan's Disjoint Forest, union find for a fixed number of vertices
 */
public class UnionFind {

    protected final int[] parent;
    protected final int[] rank;
    // nComponents will be wrong once an i or j in union have been merged more than once
    protected int nComponents = 0;
    protected final int n;

    public UnionFind(int n) {
        this.n = n;
        rank = new int[n];
        parent = new int[n];
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
        nComponents = n;
    }

    public int find(int i) {
        // path compression while searching up until parent[ii]==ii
        if (parent[i] != i) {
            parent[i] = find(parent[i]);
        }
        return parent[i];
    }

    public boolean union(int i, int j) {
        int pI = find(i);
        int pJ = find(j);
        if (pI == pJ) return false;

        if (rank[pI] > rank[pJ]) {
            parent[pJ] = pI;
        } else if (rank[pJ] > rank[pI]) {
            parent[pI] = pJ;
        } else {
            // choose pI
            parent[pJ] = pI;
            ++rank[pI];
        }
        --nComponents;
        return true;
    }

    public int[] getParent() {
        return parent;
    }

    public Map<Integer, Set<Integer>> getComponents() {
        Map<Integer, Set<Integer>> map = new HashMap();
        for (int i = 0; i < parent.length; ++i) {
            map.putIfAbsent(parent[i], new HashSet<>());
            map.get(parent[i]).add(i);
        }
        return map;
    }

}
