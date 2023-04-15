package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.util.*;

/**
 * miscellaneous graph methods.
 * 
 * TODO: consider a method to remap vertex keys for input graphs that start 
 * counting at 1 instead of 0 (they need vertex 0 removed and all vertex
 * numbers and vertex numbers in edges decreased by 1.
 * 
 * @author nichole
 */
public class GraphUtil {
    
    /**
     * create an adjacency list from the given graph knowing that the vertexes 
     * present are numbered 0 through number of vertexes-1.
     * Note that even if more than one edge is present in the direction from 
     * vertex u to vertex v, only one is present in the adjacency list for that
     * directional edge in the returned adjacency list.
     * runtime complexity is O(|V| + |E|).
     @param g
     @return 
     */
    public static SimpleLinkedListNode[] createAdjacencyList(
        NewmanGMLParser.GMLGraph g) {
        
        int nV = g.nodeIdLabelMap.size();
        
        SimpleLinkedListNode[] out = createAdjacencyList(g, nV);
        
        return out;
    }
    
    /**
     * create an adjacency list from a graph which may be missing vertex information
     * and might not be numbered from 0 to |V|-1.
     * runtime complexity is O((|V|*log_2(|V|( * |E|) which is longer 
     * due to a needed sort.
     * Note that even if more than one edge is present in the direction from 
     * vertex u to vertex v, only one is present in the adjacency list for that
     * directional edge in the returned adjacency list.
     @param g
     @return 
     */
    public static SimpleLinkedListNode[] createAdjacencyList2(
        NewmanGMLParser.GMLGraph g) {
        
        int[] vs = g.nodeIdLabelMap.keys();
        Arrays.sort(vs, 0, vs.length);
        int nV = vs[vs.length - 1] + 1;
        
        SimpleLinkedListNode[] out = createAdjacencyList(g, nV);
        
        return out;
    }
    
    private static SimpleLinkedListNode[] createAdjacencyList(
        NewmanGMLParser.GMLGraph g, final int nVertexes) {
        
        SimpleLinkedListNode[] out = new SimpleLinkedListNode[nVertexes];
        for (int v = 0; v < nVertexes; ++v) {
            out[v] = new SimpleLinkedListNode();
        }
        
        //Entry<PairInt, TFloatList>entry;
        PairInt uv;
        //TFloatList wList;
        //Iterator<Entry<PairInt, TFloatList>> iter = g.edgeWeightMap.entrySet().iterator();
        Iterator<PairInt> iter = g.edgeWeightMap.keySet().iterator();
        while (iter.hasNext()) {
            //entry = iter.next();
            //uv = entry.getKey();
            //wList = entry.getValue();
            uv = iter.next();
            out[uv.getX()].insert(uv.getY());
        }
        
        return out;
    }


    /**
     * convert the adjacency graph g in TIntObjectMap TIntSet  into a graph built with
     * SimpleLinkedListNode[].  note that this method assumes that the vertexes are ordered such
     * that the final range of indexes returned is [0, max Vertex number].
     @param g
     @return
     */
    public static SimpleLinkedListNode[] convertGraph(TIntObjectMap<TIntSet> g) {
        int[] minMax = minAndMaxVertexNumbers(g);
        int n = minMax[1] + 1;

        SimpleLinkedListNode[] g2 = new SimpleLinkedListNode[n];
        for (int i = 0; i < n; ++i) {
            g2[i] = new SimpleLinkedListNode();
        }
        TIntObjectIterator<TIntSet> iter = g.iterator();
        TIntIterator iter2;
        int u, v;
        while (iter.hasNext()) {
            iter.advance();
            u = iter.key();
            iter2 = iter.value().iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                g2[u].insert(v);
            }
        }
        return g2;
    }

    /**
     *
     @param g
     @return
     */
    public static int[] minAndMaxVertexNumbers(TIntObjectMap<TIntSet> g) {
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        TIntObjectIterator<TIntSet> iter = g.iterator();
        TIntIterator iter2;
        int u, v;
        while (iter.hasNext()) {
            iter.advance();
            u = iter.key();
            if (u < min) {
                min = u;
            }
            if (u > max) {
                max = u;
            }
            iter2 = iter.value().iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                if (v < min) {
                    min = v;
                }
                if (v > max) {
                    max = v;
                }
            }
        }
        return new int[]{min, max};
    }

    /**
     * convert the adjacency graph g in SimpleLinkedListNode[] into a
     * graph built with TIntObjectMap TIntSet .
     @param g
     @return
     */
    public static TIntObjectMap<TIntSet> convertGraph(SimpleLinkedListNode[] g) {
        TIntObjectMap<TIntSet> g2 = new TIntObjectHashMap<TIntSet>();
        int u, v;
        SimpleLinkedListNode vNode;
        TIntSet uSet;
        for (u = 0; u < g.length; ++u) {
            uSet = g2.get(u);
            if (uSet == null) {
                uSet = new TIntHashSet();
                g2.put(u, uSet);
            }
            vNode = g[u];
            while (vNode != null) {
                uSet.add(vNode.getKey());
                vNode = vNode.getNext();
            }
        }
        return g2;
    }

    /**
     * Find the bandwidth of the given graph as the maximum distance between 2 adjacent
     * vertices where the distance is the absolute difference between the index
     * numbers.
     * @param adjMap
     * @return the maximum difference between vertex numbers of adjacent vertexes.
     * returns -1 for a map without edges.
     */
    public static int measureBandwidth(Map<Integer, Set<Integer>> adjMap) {

        int max = -1;
        Iterator<Map.Entry<Integer, Set<Integer>>> iter = adjMap.entrySet().iterator();
        int u;
        int b;
        while (iter.hasNext()) {
            Map.Entry<Integer, Set<Integer>> entry = iter.next();
            u = entry.getKey();
            for (int v  : entry.getValue()) {
                b = Math.abs(u - v);
                if (b > max) {
                    max = b;
                }
            }
        }
        return max;
    }

    /**
     * Find the bandwidth of the given graph as the maximum distance between 2 adjacent
     * vertices where the distance is the absolute difference between the index
     * numbers.
     * @param adjList
     *  @return the maximum difference between vertex numbers of adjacent vertexes.
     *  returns -1 for a map without edges.
     */
    public static int measureBandwidth(SimpleLinkedListNode[] adjList) {
        int max = -1;
        int v;
        int b;
        SimpleLinkedListNode uNode;
        SimpleLinkedListNode vNode;
        for (int u = 0; u < adjList.length; ++u) {
            uNode = adjList[u];
            vNode = uNode.getNext();
            if (vNode != null && vNode.getKey() != -1) {
                v = vNode.getKey();
                b = Math.abs(u - v);
                if (b > max) {
                    max = b;
                }
            }
        }
        return max;
    }

    /**
     * Find the bandwidth of the given graph as the maximum distance between 2 adjacent
     * vertices where the distance is the absolute difference between the index
     * numbers.
     * @param adjMap
     * @return the maximum difference between vertex numbers of adjacent vertexes.
     * returns -1 for a map without edges.
     */
    public static int measureBandwidth(TIntObjectMap<TIntSet> adjMap) {

        int max = -1;
        TIntObjectIterator<TIntSet> iter = adjMap.iterator();
        int u, v;
        int b;
        TIntIterator iter2;
        while (iter.hasNext()) {
            iter.advance();
            u = iter.key();
            iter2 = iter.value().iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                b = Math.abs(u - v);
                if (b > max) {
                    max = b;
                }
            }
        }
        return max;
    }

    /**
     * Find the bandwidth of the given graph as the maximum distance between 2 adjacent
     * vertices where the distance is the absolute difference between the index
     * numbers.
     * @param edges
     * @return the maximum difference between vertex numbers of adjacent vertexes.
     * returns -1 for a map without edges.
     */
    public static int measureBandwidth(Set<PairInt> edges) {

        int max = -1;
        int u, v;
        int b;
        for (PairInt edge : edges) {
            b = Math.abs(edge.getX() - edge.getY());
            if (b > max) {
                max = b;
            }
        }
        return max;
    }

    /**
     * relabel the graph vertexes to use the numbers in rIdxs
     * @param edges
     * @param rIdxs
     * @return relabeled graph
     */
    public static Set<PairInt> relabel(Set<PairInt> edges, int[] rIdxs) {
        Set<PairInt> r = new HashSet<>();
        for (PairInt edge : edges) {
            r.add(new PairInt(rIdxs[edge.getX()], rIdxs[edge.getY()]));
        }
        return r;
    }

    /**
     * relabel the graph vertexes to use the numbers in rIdxs as reverse indexes.
     *
     * @param edges
     * @param rIdxs
     * @return relabeled graph
     */
    public static Set<PairInt> relabelWithReverse(Set<PairInt> edges, int[] rIdxs) {
        TIntIntMap iSet = new TIntIntHashMap();
        for (int i = 0; i < rIdxs.length; ++i) {
            iSet.put(rIdxs[i], i);
        }
        Set<PairInt> r = new HashSet<>();
        for (PairInt edge : edges) {
            r.add(new PairInt(iSet.get(edge.getX()), iSet.get(edge.getY())));
        }
        return r;
    }

    public static Map<Integer, Set<Integer>> copy(Map<Integer, Set<Integer>> adjMap) {
        Map<Integer, Set<Integer>> c = new HashMap<Integer, Set<Integer>>();

        Iterator<Map.Entry<Integer, Set<Integer>>> iter = adjMap.entrySet().iterator();
        while (iter.hasNext()) {
            Map.Entry<Integer, Set<Integer>> entry = iter.next();
            Set<Integer> set = new HashSet<Integer>(entry.getValue());
            c.put(entry.getKey(), set);
        }
        return c;
    }
}
