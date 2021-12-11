package algorithms.graphs;

import algorithms.graphs.ApproxGraphSearchZeng.Graph;
import algorithms.matrix.MatrixUtil;
import algorithms.sort.MultiArrayMergeSort;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import java.util.Arrays;

/**
implementing a star structure and related methods following the
 algorithms presented as pseudocode in 
 <pre> Comparing Stars: On Approximating Graph Edit Distance
   Zeng, Tung, Wang, Feng, & Zhou 2009, 
   Proceedings of the VLDB Endowment, Volume 2, Issue 1,
   August 2009 pp 25–36
   https://doi.org/10.14778/1687627.1687631
</pre>
* Additionally have added edge labels and modifications to the edit distance 
* costs for them.
 * @author nichole
 */
public class StarStructure {
    
    /**
     * ε is a special label indicating that the vertex is virtual
     */
    public static int eps = Integer.MAX_VALUE;

    /** index of root; label of root*/
    public int rootIdx;  
    public int rootLabel;
    /** the adjacent vertexes of root in V[i] treated as leaves.  
    vLabels values are labels.  vLabels is sorted by value already.
    */
    public int[] vLabels; 
    /** root -> V[i] edge labels, sorted by the vLabel sorting*/
    public int[] eLabels;
    /** the original indexes in V[i] of vLabels before sorting*/
    public int[] origVIndexes;
    /*
    vLabels[0] refers to graph g.vLabels[ origVIndexes[0] ]
    if have an original V index and want to find the star structure index for
    its vLabels or eLabels, use this mapping.
    key = original vIndex, value = local star structure vLabel or eLabel index
    */
    public TIntIntMap reverseOrigVIndexes;
    
    private StarStructure(){};
    public StarStructure(Graph g, int rootIndex) {
        init(g, rootIndex);
    }
    
    private void init(Graph g, int rootIndex) {
        TIntList vL = new TIntArrayList();
        TIntList eL = new TIntArrayList();
        TIntList vLIdxs = new TIntArrayList();
        TIntSet rNeighbors = g.adjMap.get(rootIndex);
        TIntIterator iter = rNeighbors.iterator();
        int nhbrIdx;
        while (iter.hasNext()) {
            nhbrIdx = iter.next();
            vLIdxs.add(nhbrIdx);
            vL.add(g.vLabels.get(nhbrIdx));
            eL.add(g.eLabels.get(nhbrIdx));
        }
        sortAndSet(rootIndex, g.vLabels.get(rootIndex), vL.toArray(), 
            eL.toArray(), vLIdxs.toArray());
    }
    
    public StarStructure(int rIdx, int rLabel, int[] vLabels, 
        int[] eLabels, int[] vIndexes) {
        sortAndSet(rIdx, rLabel, vLabels, eLabels, vIndexes);
    }
    
    /**
    node and edge labeling functions which can be overridden with 
    specialization.   default is to sort but not to alter the values in vLabels 
    and eLabels.  This method sorts
    vLabels, eLabels, and  the array of original indexes [0, 1, 2, ..vLabels.length - 1]
    by ascending sort on vLabels.  All results are stored in this instance.
    NOTE: the design uses specialization of labeling, but could be refactored to use 
    composition.
     * @param rIdx
     * @param rLabel
     * @param vLabels
     * @param eLabels
     * @param vIndexes indexes of the vertexes of vLabels in the context of the original graph g.
    */
    protected void sortAndSet(int rIdx, int rLabel, int[] vLabels, int[] eLabels, int[] vIndexes) {
        MultiArrayMergeSort.sortBy1stArgThen2nd(vLabels, eLabels, vIndexes);
        this.rootIdx = rIdx;
        this.rootLabel = rLabel;
        this.vLabels = vLabels;
        this.eLabels = eLabels;
        this.origVIndexes = vIndexes;
        this.reverseOrigVIndexes = new TIntIntHashMap();
        for (int i = 0; i < origVIndexes.length; ++i) {
            reverseOrigVIndexes.put(origVIndexes[i], i);
        }
    }
    
    /**
     * Lemma 4.1
     * @param s1
     * @param s2
     * @return edit distance for transforming s1 into s2.  the cost includes
     * vertex insert, delete, and substitutions, and edge substitutions.
     */
    public static int calculateEditDistance(StarStructure s1, StarStructure s2) {
        /*
        modifying Star Edit Distance to include edge labels:
            lambda(s1, s2) = T(r1, t2) + d(L1, L2)
              where T(r1, r2) = 0 if label(r1)==label(r2), else = 1
                       d(L1, L2) = ||L1|-|L2|| + M(L1, L2)
                       where psi(L) is the multiset of vertex labels on L
                       where M(L1,L2) = max( |psi(L1)|, |psi(L2)| ) - | intersection of psi(L1) with psi(L2) |
                       (M(L1,L2) is the cost of substitution of vertices).
            edge labels: node labelling difference does not use distance, so edge will not either.
            will compare edge label values after they have been sorted by the vertex label sorting.  a cost of 1 when edge labels are not equal, else the cost is 0.
            call this function d(L1E, L2E) and modified Star Edit Distance:
               lambda(s1, s2) = T(r1, t2) + d(L1, L2) + d(L1E, L2E)
            where d(L1E, L2E) = max( |psi(L1E)|, |psi(L2E)| ) - | intersection of psi(L1E) with psi(L2E) |        
        */
        int t = 0;
        //T(r1, t2)
        if (s1.rootLabel != s2.rootLabel) {
            t++;
        }
        
        int dLs = calculateEditDistanceL1L2(s1, s2);
        //lambda(s1, s2) = T(r1, t2) + d(L1, L2) + d(L1E, L2E)
        int editDist = t + dLs;

        return editDist;        
    }
    
    /**
     * calculating d(L1, L2) + d(L1E, L2E)
     * @param s1
     * @param s2
     * @return edit distance for transforming s1 into s2.  the cost includes
     * vertex insert, delete, and substitutions, and edge substitutions.
     */
    private static int calculateEditDistanceL1L2(StarStructure s1, StarStructure s2) {
        /*
        modifying Star Edit Distance to include edge labels:
            lambda(s1, s2) = T(r1, t2) + d(L1, L2)
              where T(r1, r2) = 0 if label(r1)==label(r2), else = 1
                       d(L1, L2) = ||L1|-|L2|| + M(L1, L2)
                       where psi(L) is the multiset of vertex labels on L
                       where M(L1,L2) = max( |psi(L1)|, |psi(L2)| ) - | intersection of psi(L1) with psi(L2) |
                       (M(L1,L2) is the cost of substitution of vertices).
            edge labels: node labelling difference does not use distance, so edge will not either.
            will compare edge label values after they have been sorted by the vertex label sorting.  a cost of 1 when edge labels are not equal, else the cost is 0.
            call this function d(L1E, L2E) and modified Star Edit Distance:
               lambda(s1, s2) = T(r1, t2) + d(L1, L2) + d(L1E, L2E)
            where d(L1E, L2E) = max( |psi(L1E)|, |psi(L2E)| ) - | intersection of psi(L1E) with psi(L2E) |        
        */
       
        /*editing here to perform edge value comparison at same time that
        vertex intersection is calculated.
        the intersection of the vertex labels finds vertices with same values 
            and there may be more than one vertex with same value,
            so we need to further assign the vertexes with same labels 
            by matching edge labels.
            to do so, will store the edge labels for intersection of same labeled vertexes. then
            make assignments by matching edges.
            cost for unequal edges is 1 when calculating d(L1E, L2E).
        */
        
        // map key is intersecting vertex label values
        // map value is list of the s1 edge values for the index i
        TIntObjectMap<TIntList> s1InterMap = new TIntObjectHashMap<TIntList>();
        
        // map key is intersecting vertex label values
        // map value is list of the s2 edge values for the index j
        TIntObjectMap<TIntList> s2InterMap = new TIntObjectHashMap<TIntList>();
        
        int i = 0, j = 0, nVIntersect = 0;
        while (i < s1.vLabels.length && j < s2.vLabels.length) {
            if (s1.vLabels[i] == s2.vLabels[j]) {
                if (!s1InterMap.containsKey(s1.vLabels[i])) {
                    s1InterMap.put(s1.vLabels[i], new TIntArrayList());
                    s2InterMap.put(s1.vLabels[i], new TIntArrayList());
                }
                s1InterMap.get(s1.vLabels[i]).add(s1.eLabels[i]);
                s2InterMap.get(s2.vLabels[i]).add(s2.eLabels[j]);
                //vIntersect.add(s1.vLabels[i]);
                nVIntersect++;
                i++; 
                j++;
            } else if (s1.vLabels[i] < s2.vLabels[j]) {
                i++;
            } else {
                j++;
            }
        }
        
        int pL1 = s1.vLabels.length;
        int pL2 = s2.vLabels.length;
        
        //M(L1,L2) = max( |psi(L1)|, |psi(L2)| ) - | intersection of psi(L1) with psi(L2) |
        int mL1L2 = Math.max(pL1, pL2) - nVIntersect;
        //d(L1, L2) = ||L1|-|L2|| + M(L1, L2)
        int dL1L2 = Math.abs(pL1 - pL2) + mL1L2;
        
        // further assign edges and count those that match
        TIntObjectIterator<TIntList> iter1 = s1InterMap.iterator();
        int[] s1InterE, s2InterE;
        int[] interE;
        int nEIntersect = 0;
        for (i = 0; i < s1InterMap.size(); ++i) {
            iter1.advance();
            s1InterE = iter1.value().toArray();
            s2InterE = s2InterMap.get(iter1.key()).toArray();
            // sort then count the intersection
            Arrays.sort(s1InterE);
            Arrays.sort(s2InterE);
            interE = MatrixUtil.multisetIntersection(s1InterE, s2InterE);
            nEIntersect += interE.length;
        }
        
        //M(L1E,L2E) = max( |psi(L1)|, |psi(L2)| ) - | intersection of psi(L1) with psi(L2) |
        int mL1EL2E = Math.max(pL1, pL2) - nEIntersect;
        //d(L1E, L2E) = ||L1|-|L2|| + M(L1, L2)
        int dL1EL2E = Math.abs(pL1 - pL2) + mL1EL2E;
        
        return dL1L2 + dL1EL2E;        
    }
    
    /**
     * the edit distance for use with sub-graph search.
     * section 5.2.1 of Zeng et al. 2009.
     * @param s1
     * @param s2
     * @return 
     */
    public static int calculateEditDistanceNoRelabeling(StarStructure s1, StarStructure s2) {
        int t = 0;
        //T'(r1, t2)
        if (s1.rootLabel != s2.rootLabel) {
            t = 2 + s1.vLabels.length + s2.vLabels.length;
        }
        
        int dLs = calculateEditDistanceL1L2(s1, s2);
        //lambda'(s1, s2) = T'(r1, t2) + d(L1, L2) + d(L1E, L2E)
        int editDist = t + dLs;

        return editDist;        
    }

    /**create edit distance matrix for S(g_1) to S(g_2) from 
    StarStructure.calculateEditDistance
     * @param sg1
     * @param sg2
     * @return 
    */
    public static double[][] createDistanceMatrix(StarStructure[] sg1, StarStructure[] sg2) {
        int m = sg1.length;
        int n = sg2.length;
        double[][] dist = new double[m][];
        int i, j;
        StarStructure s1, s2;
        for (i = 0; i < m; ++i) {
            dist[i] = new double[n];
            s1 = sg1[i];
            assert(s1.rootIdx == i);
            for (j = 0; j < n; ++j) {
                s2 = sg2[j];
                assert(s2.rootIdx == j);
                dist[i][j] = calculateEditDistance(s1, s2);
            }
        }
        return dist;
    }

    /**create multi-set of |V| star structures from graph g with |V| vertices
     * @param g graph holding adjacency map, vertex labels, and edge labels.
     * Note that g must use vertices 0 through g.vLabels.size() - 1, inclusive.
     * @return star structures for each vertex
     */
    public static StarStructure[] createStarStructureMultiset(Graph g) {
        int n = g.vLabels.size();
        StarStructure[] out = new StarStructure[n];
        int i, vertexIdx;
        StarStructure s;
        TIntObjectIterator<TIntSet> iter = g.adjMap.iterator();
        for (i = 0; i < n; ++i) {
            iter.advance();
            vertexIdx = iter.key();
            if (vertexIdx >= n || vertexIdx < 0) {
                throw new IllegalArgumentException("g.adjMap is inconsistent "
                        + " with g.vLabels.  g must use vertices 0 through "
                        + "g.vLabels.size() - 1, inclusive.");
            }
            s = new StarStructure(g, vertexIdx);
            out[vertexIdx] = s;
        }
        return out;
    }
    
    public static StarStructure[] copy(StarStructure[] s) {
        int n = s.length;
        StarStructure[] c = new StarStructure[n];
        int i;
        for (i = 0; i < n; ++i) {
            c[i] = s[i].copy();
        }
        return c;
    }
    
    public StarStructure copy() {
        return new StarStructure(rootIdx, rootLabel, 
            Arrays.copyOf(vLabels, vLabels.length), 
            Arrays.copyOf(eLabels, eLabels.length), 
            Arrays.copyOf(origVIndexes, origVIndexes.length));
    }
}
