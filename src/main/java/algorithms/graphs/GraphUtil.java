package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
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
}
