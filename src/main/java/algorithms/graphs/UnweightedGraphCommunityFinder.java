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
     * @return the largest modularity in the dendogram made from hierarchical
     * subgraphs made from iterative removal of the maximum edge 
     * from the graph where the score is calculated
     * by the Girvan-Newman 2002 algorithm.
     
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
    
    /**
     * a matrix based eigenvalue algorithm for dividing the graph into 2 subgraphs 
     * and iteratively their subgraphs until they are indivisible by a rule
     * based on the eigenvalues of the modularity matrix..
     * 
     * following:
     * Newman 2006, PNAS June 6, 2006 103 (23) 8577-8582
        “Modularity and community structure in networks”
        https://arxiv.org/pdf/physics/0602124.pdf
     * @param adjList
     * @return 
     */
    public static BestDivision newman2006(SimpleLinkedListNode[] adjList) {
        /*
        B_i_j(g) = A_i_j - ( (k_i * k_j)/ (2*m) ) - delta_i_j*[k_i(g) - k_i*(d_g/(2*m))
        
        Q_g = s^T * B(g) * s
        
        if Q_g does not increase, the division should not be made.
        That happens when there are not positive eigenvalues in matrix B(g)
        (--> if leading, descending ordered eigenvalue is 0, then the subgraph
        is indivisible)
        
        
        Thus our algorithm is as follows. We construct the
        modularity matrix for our network and find its leading
        (most positive) eigenvalue and eigenvector. We divide
        the network into two parts according to the signs of the
        elements of this vector, and then repeat for each of the
        parts. If at any stage we find that the proposed split
        makes a zero or negative contribution to the total modularity, 
        we leave the corresponding subgraph undivided.
        When the entire network has been decomposed into indivisible subgraphs in this way, the algorithm ends.
        One immediate corollary of this approach is that all
        “communities” in the network are, by definition, indivisible subgraphs
        */
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    private static double countEdges(SimpleLinkedListNode[] adjList) {
        
        int n = 0;
        
        for (int u = 0; u < adjList.length; ++u) {
            SimpleLinkedListNode vNode = adjList[u];
            while (vNode != null && vNode.getKey() != -1) {
                int v = vNode.getKey();
                n++;
                vNode = vNode.getNext();
            }
        }
        
        return n;
    }
    
}
