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

    /*
    TODO: consider adding the simplest spectral graph partitioning for comparison.
    caveat is that it hasn't been found to be a good partitoning of networks in general.
    reference:  Newman 2006, "Finding community structure in networks using the eigenvectors of matrices"

   (1) create Laplacian from Degree matrix and adjacency matrix:
      L = D - A
   (2) find the 2nd smallest eifenvector (the Fieldler vector) using
   CUR decomposition or other efficient matrix methods.

   (3) v_i are the normalized eigenvectors v of the Laplacian, but we only need the 2nd smallest eigenvector.

   (4) then assign membership:
       index vector s of n elements.
       s_i = +1 if vertex i is in group 1.  assign if v_i >= 0 for Fielder vector
           = -1 if vertex i is in group 2.  assign if v_i < 0 for Fielder vector
   (5) can repeat that until some criteria are met or exceeded

   see also graph cuts algorithm and girvan newman algorithm.
     */
    /**
     *
     */
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
     
     @return the largest modularity in the dendogram made from hierarchical
     * subgraphs made from iterative removal of the maximum edge 
     * from the graph where the score is calculated
     * by the Girvan-Newman 2002 algorithm.
     * 
     @param adjList unweighted graph adjacency list.
     @param src the source vertex to start calculations with if known.  has
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
            System.out.printf("k=%d q=%.3f\n", ks.get(i), qs.get(i));
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
            System.out.printf("found a peak for q=%.3f k=%d\n", qs.get(idx), ks.get(idx));
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
