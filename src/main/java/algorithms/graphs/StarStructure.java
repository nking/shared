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
   Zeng, Tung, Wang, Feng, and Zhou 2009, 
   Proceedings of the VLDB Endowment, Volume 2, Issue 1,
   August 2009 pp 25–36
   https://doi.org/10.14778/1687627.1687631
</pre>
* Additionally have added edge labels and modifications to the edit distance 
* costs for them.
* 
* For use with attributed graphs.
* 
* Fernandez and Valiente 2001:
* Attributed relational graphs, introduced in 􏶵Tsai and Fu, 1979), extend the
* notion of labeled graph by allowing the value of several attributes as labels 
* of vertices and edges.
* 
 * @author nichole
 */
public class StarStructure {
    
    /**
     * ε is a special label indicating that the vertex is virtual
     */
    public final static int eps = Integer.MAX_VALUE;

    /** index of root; label of root*/
    public int rootIdx;  
    public int rootLabel;
    /** the adjacent vertexes of root in V[i] treated as leaves.  
    vLabels values are labels.  vLabels is sorted by value already.
    */
    public int[] vLabels; 
    /** root to  V[i] edge labels, sorted by the vLabel sorting*/
    public int[] eLabels;
    /** the indexes w.r.t. the graph's vertexes, which is the indexes of the
     * array of star structures. the original indexes in V[i] are 
     * the indexes of vLabels before re-ordered by sorting */
    public int[] origVIndexes;
   
    private StarStructure(){};
    /**
     * convert the graph vertex rootIndex and it's immediate neighbors into
     * a star structure.
     * @param g
     * @param rootIndex the vertex index in g that will be the root of the star structure.
     */
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
        
        // see Zeng et al. 4.2.2, points 1-3
        
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
     * for the graph model with vertex labeling, but no edge labeling,
     * calculate the edit distance between star structures s1 and s2.
     * Lemma 4.1
     * @param s1
     * @param s2
     * @return edit distance for transforming s1 into s2.  the cost includes
     * vertex insert, delete, and substitutions, and edge substitutions.
     */
    public static int calculateEditDistanceV(StarStructure s1, StarStructure s2) {
        /*
            lambda(s1, s2) = T(r1, t2) + d(L1, L2)
              where T(r1, r2) = 0 if label(r1)==label(r2), else = 1
                       d(L1, L2) = ||L1|-|L2|| + M(L1, L2)
                       where psi(L) is the multiset of vertex labels on L
                       where M(L1,L2) = max( |psi(L1)|, |psi(L2)| ) - | intersection of psi(L1) with psi(L2) |
                       (M(L1,L2) is the cost of substitution of vertices).
               lambda(s1, s2) = T(r1, t2) + d(L1, L2)
        */
        
        int t = 0;
        //T(r1, t2)
        if (s1.rootLabel != s2.rootLabel) {
            t++;
        }
        
        int[] inters = MatrixUtil.multisetIntersection(s1.vLabels, s2.vLabels);
       
        int pL1 = s1.eLabels.length;
        int pL2 = s2.eLabels.length;
        
        //M(L1,L2) = max( |psi(L1)|, |psi(L2)| ) - | intersection of psi(L1) with psi(L2) |
        int mL1L2 = Math.max(pL1, pL2) - inters.length;
        //d(L1, L2) = ||L1|-|L2|| + M(L1, L2)
        int dL1L2 = /*Math.abs(pL1 - pL2)*/ + mL1L2; 
        int editDist = t + dL1L2;
        
        return editDist;        
    }
    
    /**
     * calculating d(L1, L2) + d(L1E, L2E)
     * @param s1
     * @param s2
     * @return edit distance for transforming s1 into s2.  the cost includes
     * vertex insert, delete, and substitutions, and edge substitutions where
     * all operations cost +1.
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
        
        // 2 more maps to count multiplicity of edge labels
        TIntIntMap s1InterCountMap = new TIntIntHashMap();
        TIntIntMap s2InterCountMap = new TIntIntHashMap();
        
        int i = 0, j = 0, nVIntersect = 0, c;
        while (i < s1.vLabels.length && j < s2.vLabels.length) {
            if (s1.vLabels[i] == s2.vLabels[j]) {
                if (!s1InterMap.containsKey(s1.vLabels[i])) {
                    s1InterMap.put(s1.vLabels[i], new TIntArrayList());
                    s2InterMap.put(s2.vLabels[j], new TIntArrayList());
                }
                s1InterMap.get(s1.vLabels[i]).add(s1.eLabels[i]);
                s2InterMap.get(s2.vLabels[j]).add(s2.eLabels[j]);
                
                // count multiplicity of edges:
                c = s1InterCountMap.containsKey(s1.eLabels[i]) ? 
                    s1InterCountMap.get(s1.eLabels[i]) : 0;
                s1InterCountMap.put(s1.eLabels[i], c + 1);
                c = s2InterCountMap.containsKey(s2.eLabels[j]) ? 
                    s2InterCountMap.get(s2.eLabels[j]) : 0;
                s2InterCountMap.put(s2.eLabels[j], c + 1);
                
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
        int dL1L2 = /*Math.abs(pL1 - pL2)*/ + mL1L2;
        
        /*
        for the edges that are in the vertex intersection (in s*InterMap),
            will calc dist using the intersection of those edges amoung themseleves.
        for the edges not in the vertex intersection,
            will calc dist using edge label intersection as is done with vertexes.
        the total dist of those two methods will be <= the number of edges.
        */
        
        // assign edges int the vertex intersection and count those that match
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
        
        int n2 = s1.eLabels.length - nVIntersect;
        if (n2 > 0) {
            j = 0;
            int[] s1NotInterE = new int[n2];
            for (i = 0; i < s1.eLabels.length; ++i) {
                if (s1InterCountMap.containsKey(s1.eLabels[i])) {
                    c = s1InterCountMap.get(s1.eLabels[i]);
                    c--;
                    if (c == 0) {
                        s1InterCountMap.remove(s1.eLabels[i]);
                    } else {
                        s1InterCountMap.put(s1.eLabels[i], c);
                    }
                    continue;
                }
                s1NotInterE[j] = s1.eLabels[i];
                j++;
            }
            j = 0;
            int[] s2NotInterE = new int[s2.eLabels.length - nVIntersect];
            for (i = 0; i < s2.eLabels.length; ++i) {
                if (s2InterCountMap.containsKey(s2.eLabels[i])) {
                    c = s2InterCountMap.get(s2.eLabels[i]);
                    c--;
                    if (c == 0) {
                        s2InterCountMap.remove(s2.eLabels[i]);
                    } else {
                        s2InterCountMap.put(s2.eLabels[i], c);
                    }
                    continue;
                }
                s2NotInterE[j] = s2.eLabels[i];
                j++;
            }
            int[] intersE2 = MatrixUtil.multisetUnorderedIntersection(s1NotInterE, s2NotInterE);

            nEIntersect += intersE2.length;
        }
        
        int mL1L2E = Math.max(pL1, pL2) - nEIntersect;

        int dL1EL2E = /*Math.abs(pL1 - pL2)*/ + mL1L2E;
        
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
    
    /**
     * for the graph model with vertex labeling, but no edge labeling, calculate
     * the edit distance for use with sub-graph search.
     * section 5.2.1 of Zeng et al. 2009.
     * @param s1
     * @param s2
     * @return 
     */
    public static int calculateEditDistanceNoRelabelingV(StarStructure s1, StarStructure s2) {
        int t = 0;
        //T'(r1, t2)
        if (s1.rootLabel != s2.rootLabel) {
            t = 2 + s1.vLabels.length + s2.vLabels.length;
        }
        
        int[] inters = MatrixUtil.multisetIntersection(s1.vLabels, s2.vLabels);
       
        int pL1 = s1.vLabels.length;
        int pL2 = s2.vLabels.length;
        
        //M(L1,L2) = max( |psi(L1)|, |psi(L2)| ) - | intersection of psi(L1) with psi(L2) |
        int mL1L2 = Math.max(pL1, pL2) - inters.length;
        //d(L1, L2) = ||L1|-|L2|| + M(L1, L2)
        int dL1L2 = /*Math.abs(pL1 - pL2)*/ + mL1L2;
        int editDist = t + dL1L2;

        return editDist;                
    }

    /**create edit distance matrix for S(g_1) to S(g_2) from 
    StarStructure.calculateEditDistance.
    * for use with approx full-graph search and graphs with vertex and edge labels.
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
    
    /**for the graph model with vertex labeling, but no edge labeling,
     * create edit distance matrix for S(g_1) to S(g_2) from 
    StarStructure.calculateEditDistance.
    * for use with approx full-graph search and graphs with vertex labels and no edge labels.
     * @param sg1
     * @param sg2
     * @return 
    */
    public static double[][] createDistanceMatrixV(StarStructure[] sg1, StarStructure[] sg2) {
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
                dist[i][j] = calculateEditDistanceV(s1, s2);
            }
        }
        return dist;
    }

    /**create edit distance matrix for S(g_1) to S(g_2) for use with
     * approx sub-graph search and graphs with vertex and edge labels.
     * @param sg1
     * @param sg2
     * @return 
    */
    public static double[][] createDistanceMatrixNoRelabeling(StarStructure[] sg1, StarStructure[] sg2) {
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
                dist[i][j] = calculateEditDistanceNoRelabeling(s1, s2);
            }
        }
        return dist;
    }
    
    /**for the graph model with vertex labeling, but no edge labeling,
     * create edit distance matrix for S(g_1) to S(g_2) from 
    StarStructure.calculateEditDistance.
    * for use with approx sub-graph search and graphs with vertex labels, but no edge labels.
     * @param sg1
     * @param sg2
     * @return 
    */
    public static double[][] createDistanceMatrixNoRelabelingV(StarStructure[] sg1, StarStructure[] sg2) {
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
                dist[i][j] = calculateEditDistanceNoRelabelingV(s1, s2);
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
        TIntObjectIterator<TIntSet> iter = g.adjMap.iterator();
        for (i = 0; i < n; ++i) {
            iter.advance();
            vertexIdx = iter.key();
            if (vertexIdx >= n || vertexIdx < 0) {
                throw new IllegalArgumentException("g.adjMap is inconsistent "
                        + " with g.vLabels.  g must use vertices 0 through "
                        + "g.vLabels.size() - 1, inclusive.");
            }
            out[vertexIdx] = new StarStructure(g, vertexIdx);
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
    
    /**
     * calculate the sum of the degree of the vertices which are adjacent to 
     * vertex v, i.e., support (v) = s(v) = ∑ d(u).
     * @param sg the graph as an array of star structures.
     * @param vIdx the index of vertex v in the star structure array sg.
     * @return 
     */
    public static int calculateSupport(StarStructure[] sg, int vIdx) {
        int s = 0, jIdx;
        for (int j = 0; j < sg[vIdx].vLabels.length; ++j) {
            jIdx = sg[vIdx].origVIndexes[j];
            s += sg[jIdx].vLabels.length;
        }
        return s;
    }
}
