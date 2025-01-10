package algorithms.graphs;

import algorithms.trees.BinaryTreeNode;
import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
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
     * convert the adjacency graph g in TIntObjectMap TIntSet  into a graph built with
     * SimpleLinkedListNode[].  note that this method assumes that the vertexes are ordered such
     * that the final range of indexes returned is [0, max Vertex number].
     @param g
     @return
     */
    public static SimpleLinkedListNode[] convertGraph(Map<Integer, Set<Integer>> g) {
        int[] minMax = minAndMaxVertexNumbers(g);
        int n = minMax[1] + 1;

        SimpleLinkedListNode[] g2 = new SimpleLinkedListNode[n];
        for (int i = 0; i < n; ++i) {
            g2[i] = new SimpleLinkedListNode();
        }

        int u;
        for (Map.Entry<Integer, Set<Integer>> entry : g.entrySet()) {
            u = entry.getKey();
            for (int v : entry.getValue()) {
                g2[u].insert(v);
            }
        }
        return g2;
    }

    public static TIntObjectMap<TIntSet> convertGraph2(Map<Integer, Set<Integer>> g) {

        TIntObjectMap<TIntSet> g2 = new TIntObjectHashMap<>();
        int u;
        for (Map.Entry<Integer, Set<Integer>> entry : g.entrySet()) {
            u = entry.getKey();
            for (int v : entry.getValue()) {
                g2.putIfAbsent(u, new TIntHashSet());
                g2.get(u).add(v);
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

    public static int[] minAndMaxVertexNumbers(Map<Integer, Set<Integer>> g) {
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        int u;
        for (Map.Entry<Integer, Set<Integer>> entry : g.entrySet()) {
            u = entry.getKey();
            min = Math.min(u, min);
            max = Math.max(u, max);
            for (int v : entry.getValue()) {
                min = Math.min(v, min);
                max = Math.max(v, max);
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
        for (Map.Entry<Integer, Set<Integer>> entry : adjMap.entrySet()) {
            Set<Integer> set = new HashSet<Integer>(entry.getValue());
            c.put(entry.getKey(), set);
        }
        return c;
    }

    public static Map<Integer, LinkedList<Integer>> convertGraph3(Map<Integer, Set<Integer>> adjMap) {
        Map<Integer, LinkedList<Integer>> c = new HashMap<>();
        for (Map.Entry<Integer, Set<Integer>> entry : adjMap.entrySet()) {
            LinkedList<Integer> list = new LinkedList<Integer>(entry.getValue());
            c.put(entry.getKey(), list);
        }
        return c;
    }

    public static Map<Integer, LinkedList<Integer>> copy2(Map<Integer, LinkedList<Integer>> adjMap) {
        Map<Integer, LinkedList<Integer>> c = new HashMap<>();
        for (Map.Entry<Integer, LinkedList<Integer>> entry : adjMap.entrySet()) {
            LinkedList<Integer> list = new LinkedList<Integer>(entry.getValue());
            c.put(entry.getKey(), list);
        }
        return c;
    }



    /**
     * find the vertex with the largest number of neighbors.
     * @param adjMap an adjacency map with key = vertex index and values = adjacent vertexes
     * @return the vertex index with the largest number of neighbors and the degree
     */
    public static int[] findMaxDegreeVertex(TIntObjectMap<TIntSet> adjMap) {
        int max = Integer.MIN_VALUE;
        int maxIdx = -1;

        TIntObjectIterator<TIntSet> iter = adjMap.iterator();
        while (iter.hasNext()) {
            iter.advance();
            if (iter.value().size() > max) {
                max = iter.value().size();
                maxIdx = iter.key();
            }
        }
        return new int[]{maxIdx, max};
    }

    public static Map<Integer, Integer> createDegreeMapForVertices(Set<Integer> vertices,
                                                                   Map<Integer, Set<Integer>> adjMap) {
        Map<Integer, Integer> degreeMap = new HashMap<Integer, Integer>();

        int nA;
        for (int v : vertices) {
            if (!adjMap.containsKey(v)) {
                nA = 0;
            } else {
                nA = adjMap.get(v).size();
            }
            degreeMap.put(v, nA);
        }
        return degreeMap;
    }

    /**
     * given graph G=(v,E) as the adjacency map adjMap, subtract vertex v from the graph.
     * Note that the vertexes adjacent to v are expected to still be present in the map as
     * keys themselves.
     * More specifically, this method removes v as key in adjMap and v in any values of adjMap.
     * revAdjMap is also updated.
     * @param adjMap
     */
    public static void subtractVertex(final int v, Map<Integer, Set<Integer>> adjMap, Map<Integer, Set<Integer>> revAdjMap) {

        if (!adjMap.containsKey(v)) {
            return;
        }
        Set<Integer> adjSet = adjMap.remove(v);
        Set<Integer> revSet = revAdjMap.remove(v);
        for (int x : adjSet) {
            if (revAdjMap.containsKey(x)) {
                revAdjMap.get(x).remove(v);
                /*if (revAdjMap.get(x).isEmpty()) {// can remove if is in color map
                    revAdjMap.remove(x);
                }*/
            }
        }
        for (int x : revSet) {
            if (adjMap.containsKey(x)) {
                adjMap.get(x).remove(v);
                /*if (adjMap.get(x).isEmpty()) {// can remove if is in color map
                    adjMap.remove(x);
                }*/
            }
        }

        /*
        v = 2
        adjMap: {1:{2, 3}, 2:{3}}
        revMap: {2:{1}, 3:{1,2}}

        adjSet = {3}
        revSet = {1}

        update revMap values using adjSet={3}: revMap: {2:{1}, 3:{1}}
        remove key '2' in revMap: {3:{1}}
        update adjMap values using revSet={1} : adjMap: {1:{3}, 2:{3}}
        remove key '2' in adjMap: {1:{3}}
        */

    }

    public static Map<Integer, Set<Integer>> createReverseMapping(Map<Integer, Set<Integer>> adjMap) {
        //adj mapping: {1:{2,3}, 2:{3}}
        //  1--> 2, 3
        //  2--> 3
        //rev mapping: {2:{1}, 3:{1,2}}
        // 2 is a neighbor of 1
        // 3 is a neighbor of 1 and 2

        Map<Integer, Set<Integer>> revMapping = new HashMap<Integer, Set<Integer>>();

        Iterator<Map.Entry<Integer, Set<Integer>>> iter = adjMap.entrySet().iterator();
        int u;
        Set<Integer> set;
        Set<Integer> vSet;
        Map.Entry<Integer, Set<Integer>> entry;
        while (iter.hasNext()) {
            entry = iter.next();
            u = entry.getKey();
            set = entry.getValue();
            for (int v : set) {
                vSet = revMapping.get(v);
                if (vSet == null) {
                    vSet = new HashSet<Integer>();
                    revMapping.put(v, vSet);
                }
                vSet.add(u);
            }
        }
        return revMapping;
    }

    /**
     * given a degree map holding key=vertex, value=degree and a set of vertexes p,
     * find the vertex in p that has the largest value in degreeMap
     * @param p
     * @param degreeMap
     * @return the vertex in p that has the largest value in degreeMap. if no members of p are in degreeMap, a -1 is returned.
     */
    public static int findMaxDegreeVertex(Set<Integer> p, Map<Integer, Integer> degreeMap) {
        int maxIdx = -1;
        int max = -1;
        int n;
        for (int v : p) {
            if (!degreeMap.containsKey(v)) {
                continue;
            }
            n = degreeMap.get(v);
            if (max < n) {
                max = n;
                maxIdx = v;
            }
        }
        return maxIdx;
    }

    /**
     * given degreeMap holding key=vertex, value=degree, return the vertex with largest value
     * @param degreeMap
     * @return the vertex with largest value
     */
    public static int findMaxDegreeVertex(Map<Integer, Integer> degreeMap) {
        int maxIdx = -1;
        int max = -1;
        Iterator<Map.Entry<Integer, Integer>> iter = degreeMap.entrySet().iterator();
        Map.Entry<Integer, Integer> entry;
        int n;
        while (iter.hasNext()) {
            entry = iter.next();
            n = entry.getValue();
            if (max < n) {
                max = n;
                maxIdx = entry.getKey();
            }
        }
        return maxIdx;
    }

    public static boolean addEdge(int v, int w, Map<Integer, Set<Integer>> adjMap) {
        Set<Integer> set = adjMap.get(v);
        if (set == null) {
            set = new HashSet<Integer>();
            adjMap.put(v, set);
        }
        return set.add(w);
    }
    public static boolean removeEdge(int v, int w, Map<Integer, Set<Integer>> adjMap) {
        Set<Integer> set = adjMap.get(v);
        if (set == null) {
            return false;
        }
        return set.remove(w);
    }

    /**
     * runtime complexity is O(N*log_2(N)) because it sorts the adjacency list if sort=true
     * @param adjMap
     * @return
     */
    public static TIntObjectMap<TIntList> copyToOrderedAdjMap(Map<Integer, Set<Integer>> adjMap, boolean sort) {
        TIntObjectMap<TIntList> map = new TIntObjectHashMap<TIntList>();
        Iterator<Map.Entry<Integer, Set<Integer>>> iter = adjMap.entrySet().iterator();
        Set<Integer> set;
        TIntList list;
        int u;
        Map.Entry<Integer, Set<Integer>> entry;
        while (iter.hasNext()) {
            entry = iter.next();
            u = entry.getKey();
            set = entry.getValue();
            list = new TIntArrayList();
            for (int v : set) {
                list.add(v);
            }
            if (sort) {
                list.sort();
            }
            map.put(u, list);
        }
        return map;
    }

    public static Set<PairInt> extractEdges(Map<Integer, Set<Integer>> adjMap) {
        Set<PairInt> edges = new HashSet<PairInt>();
        Iterator<Map.Entry<Integer, Set<Integer>>> iter = adjMap.entrySet().iterator();
        int u;
        Map.Entry<Integer, Set<Integer>> entry;
        while (iter.hasNext()) {
            entry = iter.next();
            u = entry.getKey();
            for (int v : entry.getValue()) {
                edges.add(new PairInt(u, v));
            }
        }
        return edges;
    }
    public static Set<PairInt> extractEdgesUsingLexicographicOrder(Map<Integer, Set<Integer>> adjMap) {
        Set<PairInt> edges = new HashSet<PairInt>();
        Iterator<Map.Entry<Integer, Set<Integer>>> iter = adjMap.entrySet().iterator();
        int u;
        Map.Entry<Integer, Set<Integer>> entry;
        while (iter.hasNext()) {
            entry = iter.next();
            u = entry.getKey();
            for (int v : entry.getValue()) {
                if (u <= v) {
                    edges.add(new PairInt(u, v));
                } else {
                    edges.add(new PairInt(v, u));
                }
            }
        }
        return edges;
    }

    public <T> int countNodes(Map<Integer, Map<Integer, T>> adjMap) {
        Set<Integer> nodes = new HashSet<>();
        for (int u : adjMap.keySet()) {
            nodes.add(u);
            for (int v : adjMap.get(u).keySet()) {
                nodes.add(v);
            }
        }
        return nodes.size();
    }

    protected static Map<Integer, List<double[]>> createSortedAdjList(int nNodes, Map<Integer, Map<Integer, Double>> adjMap0) {

        Map<Integer, List<double[]>> adjMap = new HashMap<>();

        for (int u : adjMap0.keySet()) {
            for (Map.Entry<Integer, Double> entry : adjMap0.get(u).entrySet()) {
                int v = entry.getKey();
                double w = entry.getValue();
                adjMap.putIfAbsent(u, new ArrayList<double[]>());
                adjMap.putIfAbsent(v, new ArrayList<double[]>());
                adjMap.get(u).add(new double[]{v, w});
                adjMap.get(v).add(new double[]{u, w});
            }
        }

        for (int key : adjMap.keySet()) {
            Collections.sort(adjMap.get(key), (o1, o2) -> Double.compare(o1[1], o2[1]));
        }

        return adjMap;
    }

    public static int[] findMinWeightEdge(Map<Integer, Map<Integer, Double>> adjMap) {
        double min = Double.POSITIVE_INFINITY;
        int[] minEdge = new int[2];
        for (int u : adjMap.keySet()) {
            for (Map.Entry<Integer, Double> entry : adjMap.get(u).entrySet()) {
                if (entry.getValue() < min) {
                    min = entry.getValue();
                    minEdge[0] = u;
                    minEdge[1] = entry.getKey();
                }
            }
        }
        return minEdge;
    }

    public <T> int countNodes2(Map<Integer, Set<Integer>> adjMap) {
        Set<Integer> nodes = new HashSet<>();
        for (int u : adjMap.keySet()) {
            nodes.add(u);
            for (int v : adjMap.get(u)) {
                nodes.add(v);
            }
        }
        return nodes.size();
    }
}
