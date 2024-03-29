package algorithms.graphs;

import algorithms.misc.MiscMath0;
import algorithms.sort.MiscSorter;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.HashSet;
import java.util.Set;

/**
The Reverse Cuthill-McKee (RCM) is an algorithm to reduce the bandwidth of a
symmetric matrix (remember that an adjacency matrix for an undirected graph
is symmetric).   The reduction of the bandwidth of a matrix
reduces storage and computational costs.
 
From "The Reverse Cuthill-McKee Algorithm in Distributed-Memory"
2016 Azad, Jacquelin, Buluc and Ng
 
Since obtaining a reordering to minimize bandwidth is an NP-complete problem, 
various heuristics are used in practice such as Cuthill-McKee, 
Reverse Cuthill-McKee (RCM), and Sloan’s algorithms [4], [5], [6]. 
This paper solely focuses on the RCM algorithm [5] because, 
with careful algorithm design, it is amenable to 
massive distributed-memory parallelism – the primary topic of interest of this paper.

 [4] E. Cuthill and J. McKee, “Reducing the bandwidth of sparse symmetric matrices,” in Proc. of 24th national conference. ACM, 1969, pp. 157– 172.
 [5] A.GeorgeandJ.W.-H.Liu,Computer Solution of Large Sparse Positive Definite Systems. Englewood Cliffs, New Jersey: Prentice-Hall Inc., 1981.
 [6] S. Sloan, “An algorithm for profile and wavefront reduction of sparse matrices,” International Journal for Numerical Methods in Engineering, vol. 23, no. 2, pp. 239–251, 1986.

TODO: consider implementing the distributed-memory version of Azadi et al. 2016.

 see also skyline matrix storage
 https://en.wikipedia.org/wiki/Skyline_matrix

@author nichole
 */
public class ReverseCuthillMcKee {
    
    /**
     * given the adjacency map as pairs of edges, calculate the
     * Cuthill-McKee ordering to reduce the bandwidth of a
     * symmetric matrix (remember that an adjacency matrix for an undirected graph
     * is symmetric).
     * runtime complexity is O(|V| + |E|*log_2(|E|) where |V| is the number
     * of vertices and |E| is the number of edges where an edge is counted
     * once.
     * <pre>
     * references:
     *https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm
     * 
     * The Reverse Cuthill-McKee Algorithm in Distributed-Memory 
       2016 Azad, Jacquelin, Buluc and Ng

     * </pre>
     @param gE undirected adjacency graph
     @return Cuthill-McKee ordering
     */
    public static int[] cuthillMcKeeOrdering(Set<PairInt> gE) {
        
        TIntObjectMap<TIntSet> adjMap = createSymmetricAdjMap(gE);
        
        return cuthillMcKeeOrdering(adjMap);
    }
    
    /**
     * given the adjacency map as pairs of edges, calculate the reverse
     * Cuthill-McKee ordering.
     * runtime complexity is O(|V| + |E|*log_2(|E|) where |V| is the number
     * of vertices and |E| is the number of edges where an edge is counted
     * once.
     *
     * Reducing the bandwidth of a graph consists of ﬁnding a special labeling of vertices which minimizes the
     * maximum absolute difference between the labels of adjacent vertices.
     * <pre>
     * references:
     * 
     * https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm
     * 
     * The Reverse Cuthill-McKee Algorithm in Distributed-Memory 
       2016 Azad, Jacquelin, Buluc and Ng

     The Bandwidths of a Matrix. A Survey of Algorithms, 2014 Maftelu-Scai
     * </pre>
     *
     @param gE undirected adjacency graph
     @return reverse Cuthill-McKee ordering
     */
    public static int[] rcm(Set<PairInt> gE) {
        
        TIntObjectMap<TIntSet> adjMap = createSymmetricAdjMap(gE);
        
        int[] c = cuthillMcKeeOrdering(adjMap);
       
        MiscMath0.reverse(c);
        int[] r = new int[c.length];
        for (int i = 0; i < c.length; ++i) {
            r[c[i]] = i;
        }

        return r;
    }
    
    /**
     * given the adjacency map of an undirected graph, calculate the
     * Cuthill-McKee ordering.
     * runtime complexity is O(|V| + |E|*log_2(|E|) where |V| is the number
     * of vertices and |E| is the number of edges where an edge is counted
     * once.
     * <pre>
     * references:
     * 
     * https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm
     * 
     * The Reverse Cuthill-McKee Algorithm in Distributed-Memory 
       2016 Azad, Jacquelin, Buluc and Ng

     * </pre>
     @param adjMap undirected adjacency graph
     @return reverse Cuthill-McKee ordering
     */
    public static int[] rcm(TIntObjectMap<TIntSet> adjMap) {
        
        int[] c = cuthillMcKeeOrdering(adjMap);
        
        MiscMath0.reverse(c);
        int[] r = new int[c.length];
        for (int i = 0; i < c.length; ++i) {
            r[c[i]] = i;
        }
        
        return r;
    }
    
    /**
     * given the adjacency map of an undirected graph, calculate the
     * Cuthill-McKee ordering.
     * runtime complexity is O(|V| + |E|*log_2(|E|) where |V| is the number
     * of vertices and |E| is the number of edges where an edge is counted
     * once.
     * <pre>
     * references:
     * 
     * https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm
     * 
     * The Reverse Cuthill-McKee Algorithm in Distributed-Memory 
       2016 Azad, Jacquelin, Buluc and Ng

     * </pre>
     @param adjMap undirected adjacency graph
     @return Cuthill-McKee ordering
     */
    public static int[] cuthillMcKeeOrdering(TIntObjectMap<TIntSet> adjMap) {
                
        // peripheral vertex = vertex with lowest degree
        Set<PairInt> unique = uniqueEdges(adjMap);
        TIntIntMap vertexDegreeMap = createVertexDegreeMap(unique);
        int n = vertexDegreeMap.size();
        //int m = unique.size();
        
        int pIdx = findMinDegreeVertex(vertexDegreeMap);
        
        int[] out = new int[n];
        TIntSet outSet = new TIntHashSet();
        int oIdx = 0;
        out[oIdx] = pIdx;        
        outSet.add(out[oIdx]);
        
        // reversed out for use in sort by earliest predecessor 
        int[] revOut = new int[n];
        revOut[pIdx] = -1;
        
        oIdx++;
        
        ArrayDeque<Integer> q = new ArrayDeque<>();
        q.add(pIdx);
        
        TIntSet nhbrs;
        TIntList adj = new TIntArrayList();
        int[] adjSorted;
        TIntIterator iter;
        int nIdx, i;
        // runtime complexity is O(|V| + |E|*log_2(|E|)
        while (outSet.size() < n) {
            if (q.isEmpty()) {
                // find minimum degree vertex which is not in outSet
                pIdx = findMinDegreeVertex(vertexDegreeMap, outSet);
                assert(pIdx > -1);
                out[oIdx] = pIdx;
                revOut[pIdx] = oIdx;
                outSet.add(pIdx);
                oIdx++;
            } else {
                pIdx = q.poll();
            }
            
            adj.clear();
                        
            // gather neighbors that aren't in outSet and sort them by
            //  their revOut values w/ ties broken by smallest vertex degree
            nhbrs = adjMap.get(pIdx);
            iter = nhbrs.iterator();
            while (iter.hasNext()) {
                nIdx = iter.next();
                if (!outSet.contains(nIdx)) {
                    adj.add(nIdx);
                }
            }
            if (adj.isEmpty()) {
                continue;
            }
            adjSorted = sort(adj, revOut, vertexDegreeMap);
            
            for (i = 0; i < adjSorted.length; ++i) {
                out[oIdx] = adjSorted[i];
                revOut[adjSorted[i]] = oIdx;
                outSet.add(adjSorted[i]);
                q.add(adjSorted[i]);
                oIdx++;
            }
        }
        
        return out;
    }
        
    /**
     given the adjacency matrix of an undirected graph (a is symmetric), calculate the reverse
     * Cuthill-McKee ordering.
     * runtime complexity is O(|V| + |E|*log_2(|E|) where |V| is the number
     * of vertices and |E| is the number of edges where an edge is counted
     * once.
     * <pre>
     * references:
     * 
     * https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm
     * 
     * The Reverse Cuthill-McKee Algorithm in Distributed-Memory 
       2016 Azad, Jacquelin, Buluc and Ng

     * </pre>
     @param a symmetric adjacency matrix where entry a[i][j] .gt. 0 indicates an edge
     * between vertexes i and j and the graph is undirected.
     @return reverse Cuthill-McKee ordering
     */
    public static int[] rcm(int[][] a) {
        
        TIntObjectMap<TIntSet> adjMap = createSymmetricAdjMap(a);
        
        int[] c = cuthillMcKeeOrdering(adjMap);
        
        MiscMath0.reverse(c);
        int[] r = new int[c.length];
        for (int i = 0; i < c.length; ++i) {
            r[c[i]] = i;
        }
        
        return r;
    }
    
    /**
     given the adjacency matrix of an undirected graph (a is symmetric), calculate the reverse
     * Cuthill-McKee ordering.
     * runtime complexity is O(|V| + |E|*log_2(|E|) where |V| is the number
     * of vertices and |E| is the number of edges where an edge is counted
     * once.
     * <pre>
     * references:
     * 
     * https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm
     * 
     * The Reverse Cuthill-McKee Algorithm in Distributed-Memory 
       2016 Azad, Jacquelin, Buluc and Ng

     * </pre>
     @param a symmetric adjacency matrix where entry a[i][j] .gt. 0 indicates an edge
     * between vertexes i and j and the graph is undirected.
     @return Cuthill-McKee ordering
     */
    public static int[] cuthillMcKeeOrdering(int[][] a) {
        TIntObjectMap<TIntSet> adjMap = createSymmetricAdjMap(a);
        return cuthillMcKeeOrdering(adjMap);
    }
    
    /**
     * rewrite the adjacency map edges where for all (u, v) to pairs such that u is less than v
     * (u and v are vertex indexes).  The rewriting is to avoid double counting
     * in other methods.
     @param adjMap
     @return the adjacency map rewritten so that u is less than v
     */
    static Set<PairInt> uniqueEdges(TIntObjectMap<TIntSet> adjMap) {
    
        TIntObjectIterator<TIntSet> iter = adjMap.iterator();
        
        Set<PairInt> out = new HashSet<PairInt>();
        int u, v;
        TIntSet set;
        TIntIterator iter2;
        while (iter.hasNext()) {
            iter.advance();
            u = iter.key();
            set = iter.value();
            iter2 = set.iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                if (u < v) {
                    out.add(new PairInt(u, v));
                } else {
                    out.add(new PairInt(v, u));
                }
            }
        }
        return out;
    }
    
    static TIntIntMap createVertexDegreeMap(Set<PairInt> gE) {
        int u, v;
        TIntIntMap vertexDegreeMap = new TIntIntHashMap();
        for (PairInt uv : gE) {
            u = uv.getX();
            v = uv.getY();
            if (vertexDegreeMap.containsKey(u)) {
                vertexDegreeMap.put(u, vertexDegreeMap.get(u) + 1);
            } else {
                vertexDegreeMap.put(u, 1);
            }
            if (vertexDegreeMap.containsKey(v)) {
                vertexDegreeMap.put(v, vertexDegreeMap.get(v) + 1);
            } else {
                vertexDegreeMap.put(v, 1);
            }
        }
        return vertexDegreeMap;
    }
    
    static int findMinDegreeVertex(TIntIntMap vertexDegreeMap) {
        TIntIntIterator iter = vertexDegreeMap.iterator();
        int minD = Integer.MAX_VALUE;
        int minV = -1;
        while (iter.hasNext()) {
            iter.advance();
            if (iter.value() < minD) {
                minV = iter.key();
                minD = iter.value();
            }
        }
        return minV;
    }
    
    static int findMinDegreeVertex(TIntIntMap vertexDegreeMap, TIntSet outSet) {
        TIntIntIterator iter = vertexDegreeMap.iterator();
        int minD = Integer.MAX_VALUE;
        int minV = -1;
        int v;
        while (iter.hasNext()) {
            iter.advance();
            v = iter.key();
            if (outSet.contains(v)) {
                continue;
            }
            if (iter.value() < minD) {
                minV = iter.key();
                minD = iter.value();
            }
        }
        return minV;
    }

    static TIntObjectMap<TIntSet> createSymmetricAdjMap(Set<PairInt> gE) {
        
        TIntObjectMap<TIntSet> out = new TIntObjectHashMap<TIntSet>();
        
        int u, v;
        TIntSet set;
        for (PairInt uv : gE) {
            u = uv.getX();
            v = uv.getY();
            
            set = out.get(u);
            if (set == null) {
                set = new TIntHashSet();
                out.put(u, set);
            }
            set.add(v);
            
            set = out.get(v);
            if (set == null) {
                set = new TIntHashSet();
                out.put(v, set);
            }
            set.add(u);
        }
        
        return out;
    }

    static int[] sort(TIntList adj, int[] revOut, TIntIntMap vertexDegreeMap) {
        //sort by smallest revOut values w/ ties broken by smallest vertex degree
        int[] a = new int[adj.size()];
        int[] b = new int[adj.size()];
        int[] idxs = new int[adj.size()];
        int i, v;
        for (i = 0; i < adj.size(); ++i) {
            v = adj.get(i);
            a[i] = revOut[v];
            b[i] = vertexDegreeMap.get(v);
            idxs[i] = i;
        }
        
        MiscSorter.sortBy1stArgThen2nd(a, b, idxs);
        
        int[] out = new int[idxs.length];
        for (i = 0; i < idxs.length; ++i) {
            out[i] = adj.get(idxs[i]);
        }
        
        return out;
    }

    static TIntObjectMap<TIntSet> createSymmetricAdjMap(int[][] a) {
        
        TIntObjectMap<TIntSet> out = new TIntObjectHashMap<TIntSet>();
        
        TIntSet set;
        int u, v;
        for (u = 0; u < a.length; ++u) {
            for (v = 0; v < a[u].length; ++v) {
                if (a[u][v] < 1) {
                    continue;
                }
                set = out.get(u);
                if (set == null) {
                    set = new TIntHashSet();
                    out.put(u, set);
                }
                set.add(v);

                set = out.get(v);
                if (set == null) {
                    set = new TIntHashSet();
                    out.put(v, set);
                }
                set.add(u);
            }
        }
        return out;
    }

}
