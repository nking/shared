package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.util.Arrays;
import java.util.Iterator;

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
     * @param g
     * @return 
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
     * @param g
     * @return 
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
     * convert the adjacency graph g in TIntObjectMap<TIntSet> into a graph built with
     * SimpleLinkedListNode[].  note that this method assumes that the vertexes are ordered such
     * that the final range of indexes returned is [0, max Vertex number].
     * @param g
     * @return
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
     * graph built with TIntObjectMap<TIntSet>.
     * @param g
     * @return
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

}
