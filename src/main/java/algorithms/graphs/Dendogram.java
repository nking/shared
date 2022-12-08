package algorithms.graphs;

import algorithms.disjointSets.DisjointForest;
import algorithms.graphs.Betweenness.Results;
import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TObjectFloatIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * encapsulates dendogram layers for each number of communities as
 * the graph is increasingly partitioned into larger numbers of 
 * communities.
 * A dendogram layer is the number of components (a.k.a. communities) 
 * in graph layer.
   
 * @author nichole
 */
public class Dendogram {

    private final List<DendogramLayer> layers;
    
    private final SimpleLinkedListNode[] originalGraph;
    
    private Dendogram() {originalGraph=null;layers=null;}
    
    /**
     *
     @param adjList
     */
    public Dendogram(final SimpleLinkedListNode[] adjList) {
        originalGraph = copyGraph(adjList);
        layers = new ArrayList<DendogramLayer>(adjList.length);
    }
    
    private SimpleLinkedListNode[] copyGraph(final SimpleLinkedListNode[] g) {
        SimpleLinkedListNode[] c = g.clone();
        for (int i = 0; i < g.length; ++i) {
            c[i] = new SimpleLinkedListNode(g[i]);
        }
        return c;
    }
    
    /**
     *
     @param src
     @return
     */
    public int createUsingGirvanNewman(int src) {
        
        if (originalGraph == null || layers == null) {
            throw new IllegalArgumentException("instance must be constructed with"
                    + " adjacency matrix");
        }
        
        int nV = originalGraph.length;
        
        Betweenness b = new Betweenness();
        SimpleLinkedListNode[] g = copyGraph(originalGraph);

        List<TIntSet> cc = DisjointForest.connectedComponents(g);
        int nC = cc.size();
        DendogramLayer layer = createDendogramLayer(cc, nV);
        getLayers().add(layer);
        
        while (nC < nV) {
            
            System.out.printf("nComponents=%d (nV=%d)\n", nC, nV);
            System.out.flush();

            Results r = b.girvanNewmanDAG(g, src);

            TObjectFloatMap<PairInt> edgeWeights = r.getEdges();
            PairInt max = getMaxEdge(edgeWeights);
            if (max == null) {
                break;
            }
            System.out.printf("max edge = (%d, %d)\n", max.getX(), max.getY());
            
            g[max.getX()].delete(max.getY());

            cc = DisjointForest.connectedComponents(g);
            nC = cc.size();
            layer = createDendogramLayer(cc, nV);
            getLayers().add(layer);
        }
        
        return nC;
    }
    
    /**
     *
     */
    public static class DendogramLayer {
        /**
         * the number of components (a.k.a. communities) in graph layer
         */
        int nComponents;
        
        /**
         * index is vertex number, value is component number.  if value is -1,
         * there is no assigned component.
         */
        int[] vertexComponents;
        
        /**
         * map with key = component number, value = set of vertex numbers in the component
         */
        TIntObjectMap<TIntSet> componentVertexes;
    }
    
    private DendogramLayer createDendogramLayer(List<TIntSet> cc,int nVertexes) {
        DendogramLayer layer = new DendogramLayer();
        layer.nComponents = cc.size();
        layer.vertexComponents = new int[nVertexes];
        Arrays.fill(layer.vertexComponents, -1);
        layer.componentVertexes = new TIntObjectHashMap<TIntSet>();
        
        TIntIterator iter;
        int v;
        for (int c = 0; c < cc.size(); ++c) {
            
            TIntSet set = cc.get(c);
            
            layer.componentVertexes.put(c, set);
            
            iter = set.iterator();
            while (iter.hasNext()){
                v = iter.next();
                layer.vertexComponents[v] = c;
            }    
        }
        
        return layer;
    }

    private PairInt getMaxEdge(TObjectFloatMap<PairInt> edgeWeights) {
        
        TObjectFloatIterator<PairInt> iter = edgeWeights.iterator();
        
        float max = Float.NEGATIVE_INFINITY;
        PairInt maxP = null;
        float w;
        
        for (int i = 0; i < edgeWeights.size(); ++i) {
            iter.advance();
            w = iter.value();
            if (w > max) {
                max = w;
                maxP = iter.key();
            }
        }
        
        return maxP;
    }

    /**
     @return the layers
     */
    public List<DendogramLayer> getLayers() {
        return layers;
    }
    
}
