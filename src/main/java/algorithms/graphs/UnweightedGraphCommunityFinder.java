package algorithms.graphs;

import algorithms.graphs.Dendogram.DendogramLayer;
import algorithms.misc.MinMaxPeakFinder;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.util.List;

/**
 *
 * @author nichole
 */
public class UnweightedGraphCommunityFinder {

    public static class BestDivision {
        DendogramLayer layer;
        float modularity;
    }
    
    /**
     * given the original graph adjacency list, compute the modularity of the
     * graph using the connected components from an iteration of the girvan-newman
     * algorithm for edge betweenness scores.
     Reference: Newman and Girvan 2004, PHYSICAL REVIEW E 69, 026113,
     "Finding and evaluating community structure in networks"

     NOTE: If the number of within-community edges is no better than random,
     * we will get Q = 0. value approaching the maximum, Q = 1, indicate strong
     * community structure [50]. values tend to be 0.3 to 0.7 and higher values
     * are rare.
     * 
     * NOTE: see notes for Newman 2006 spectral modularity algorithms in
     * this project's directory docs/miscNotes.
     
     * @return the largest modularity in the dendogram made from hierarchical
     * subgraphs made from iterative removal of the maximum edge 
     * from the graph where the score is calculated
     * by the Girvan-Newman 2002 algorithm.
     * 
     * @param adjList unweighted graph adjacency list.
     * @param src the source vertex to start calculations with if known.  has
     * to be a valid vertex in the graph, but a random choice does affect the 
     * results adversely.
     */
    public static BestDivision girvanNewman2002(SimpleLinkedListNode[] adjList,
        int src) {
                
        Dendogram d = new Dendogram(adjList);
        
        int kFinal = d.createUsingGirvanNewman(src);
        
        List<Dendogram.DendogramLayer> layers = d.getLayers();
        
        double nEdges = countEdges(adjList);
        
        // store q and k to look for peak in q
        TIntList ks = new TIntArrayList();
        TFloatList qs = new TFloatArrayList();
        
        Modularity m = new Modularity();
        
        double[][] e;
        int i, j, k;
        double q;
        for (Dendogram.DendogramLayer layer : layers) {
            
            q = m.girvanNewman2002(layer, adjList, nEdges);
            k = layer.nComponents;
            
            ks.add(k);
            qs.add((float)q);
        }
        
        for (i = 0; i < ks.size(); ++i) {
            System.out.printf("k=%d q=%.3f%n", ks.get(i), qs.get(i));
        }
        
        // find the max among the peaks
        MinMaxPeakFinder mpf = new MinMaxPeakFinder();
        int[] maxIdxs = mpf.findPeaks(qs.toArray(), 0.f, 2.5f);
        float maxQ = Float.NEGATIVE_INFINITY;
        int maxQIdx = -1;
        for (i = 0; i < maxIdxs.length; ++i) {
            int idx = maxIdxs[i];
            if (ks.get(idx) > maxQ) {
                maxQ = ks.get(idx);
                maxQIdx = idx;
            }
            System.out.printf("found a peak for q=%.3f k=%d%n", qs.get(idx), ks.get(idx));
        }
        
        if (maxQIdx == -1) {
            return null;
        }
        
        BestDivision best = new BestDivision();
        best.layer = layers.get(maxQIdx);
        best.modularity = maxQ;
        return best;
    }
    
    private static double countEdges(SimpleLinkedListNode[] adjList) {
        
        int n = 0;
        
        for (int u = 0; u < adjList.length; ++u) {
            SimpleLinkedListNode vNode = adjList[u];
            while (vNode != null && vNode.getNumberOfKeys() > 0) {
                int v = vNode.getKey();
                n++;
                vNode = vNode.getNext();
            }
        }
        
        return n;
    }
    
}
