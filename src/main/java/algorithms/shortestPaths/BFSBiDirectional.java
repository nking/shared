package algorithms.shortestPaths;

import java.util.*;

/**
 * perform Breadth First search bi-directionally from src to destination, a.k.a. 2-end BFS.
 *
 * Note that the search by one direction, then the other, was adapted from the book
 * "Cracking the Code Interview, 6th Ed. by Gayle Laakmann McDowell.
 *
 * Bi-directional BFS runtime complexity is O(b^(d/2)) where b is the branching factor and d is the distance from
 * start to destination.
 * BFS without the bi-directional search is O(b^d).
 *
 * Bi-directional search works best when the branching factors from both ends are similar.
 *
 * @author Nichole
 */
public class BFSBiDirectional {
    Map<Integer, Integer> visitedS;
    Map<Integer, Integer> visitedD;
    Map<Integer, Integer> dS;
    Map<Integer, Integer> dD;
    Map<Integer, Integer> prevS;
    Map<Integer, Integer> prevD;
    Map<Integer, Set<Integer>> adjMap;

    public List<Integer> search(Map<Integer, Set<Integer>> adjMap, int src, int dest) {
        if (adjMap == null || adjMap.isEmpty()) {
            throw new IllegalArgumentException("adjMap cannot be null or empty");
        }
        this.adjMap = adjMap;
        // init
        int n = adjMap.size();
        visitedS = new HashMap<>();
        visitedD = new HashMap<>();
        dS = new HashMap<>();
        dD = new HashMap<>();
        prevS = new HashMap<>();
        prevD = new HashMap<>();

        java.util.Queue<Integer> qS = new ArrayDeque<Integer>();
        java.util.Queue<Integer> qD = new ArrayDeque<Integer>();
        visitedS.put(src, 1);
        dS.put(src, 0);
        visitedD.put(dest, 1);
        dD.put(dest, 0);
        qS.offer(src);
        qD.offer(dest);

        while (!qS.isEmpty() && !qD.isEmpty()) {
            Integer collisionId = searchLevel(qS, visitedS, visitedD, dS, prevS);
            if (collisionId != null) {
                return mergePaths(collisionId);
            }
            collisionId = searchLevel(qD, visitedD, visitedS, dD, prevD);
            if (collisionId != null) {
                return mergePaths(collisionId);
            }
        }
        return null;
    }

    protected Integer searchLevel(Queue<Integer> q0,
                                  Map<Integer, Integer> visited0, Map<Integer, Integer> visited1,
                                  Map<Integer, Integer> d0,
                                  Map<Integer, Integer> prev0) {

        int nS = q0.size();
        for (int i = 0; i < nS; i++) {
            int u = q0.poll();
            if (visited1.containsKey(u)) {
                return u;
            }
            Set<Integer> vS = adjMap.get(u);
            if (vS != null) {
                for (int v : vS) {
                    if (!d0.containsKey(v) || (d0.get(v) > d0.get(u) + 1)) {
                        d0.put(v, d0.get(u) + 1);
                        prev0.put(v, u);
                        if (!visited0.containsKey(v)) {
                            visited0.put(v, 1);
                            q0.offer(v);
                        }
                    }
                }
            }
        }
        return null;
    }

    public List<Integer> mergePaths(int connection) {
        List<Integer> path = new ArrayList<>();
        Integer prev = connection;
        while (prev != null) {
            path.add(prev);
            prev = prevS.get(prev);
        }
        // reverse order to put src at beginning
        int n = path.size();
        int swap;
        for (int i = 0; i < n/2; ++i) {
            swap = path.get(n-i-1);
            path.set(n-i-1, path.get(i));
            path.set(i, swap);
        }
        prev = prevD.get(connection);
        while (prev != null) {
            path.add(prev);
            prev = prevD.get(prev);
        }
        return path;
    }
}
