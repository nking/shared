package algorithms.graphs;

import algorithms.matrix.MatrixUtil;
import algorithms.util.PairInt;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.set.TIntSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import thirdparty.graphMatchingToolkit.algorithms.VolgenantJonker;

/**
 implementing subgraph and full-graph query methods following the
 algorithms presented as pseudocode in 
 <pre> Comparing Stars: On Approximating Graph Edit Distance
   Zeng, Tung, Wang, Feng, & Zhou 2009, 
   Proceedings of the VLDB Endowment, Volume 2, Issue 1,
   August 2009 pp 25–36
   https://doi.org/10.14778/1687627.1687631
</pre>
subgraph search: retrieve all the super-graphs of a given query graph Q from a 
   graph database D = (D_1,D_2,...D_n).
   It finds all the graphs D_i(i=1,2,…,m.lt.n) such that Q is a subgraph 
   of D_i.  It also finds D_i that are subgraphs of Q.
 
full graph search: find all graphs D_i in D s.t. D_i is the same as the query graph Q.

similarity search: find all graphs D_i in D s.t. D_i is similar to the
   query graph Q within a user-specified threshold based on some similarity measures.
 
 * @author nichole
 */
public class ApproxGraphSearchZeng {
    
    /**
     * property used in calculating the edit distance.  if edgesAreLabeled = true,
     * the cost of edge insert, delete, and substitutions are added.
     */
    private boolean edgesAreLabeled = false;
    /**
     * set the property edgesAreLabeled to true or false (default is false),
     * to add the cost of edge insert, delete, and substitutions into edit distances.
     */
    public void setEdgesAreLabeled(boolean labeled) {
        this.edgesAreLabeled = labeled;
    }
    
    /**
     * find all graphs db_i in db s.t. db_i is the same as the query graph q
     * within a graph edit distance threshold w.
     * implementing full-graph query method following the
     algorithms presented as pseudocode in 
     <pre> Comparing Stars: On Approximating Graph Edit Distance
       Zeng, Tung, Wang, Feng, & Zhou 2009, 
       Proceedings of the VLDB Endowment, Volume 2, Issue 1,
       August 2009 pp 25–36
       https://doi.org/10.14778/1687627.1687631
    </pre>
    The method attempts to avoid the expensive graph edit distance computation 
    by filtering out graphs that definitely will not be in the answer set.
    The best-case runtime complexity is Θ(n^3), while the worse is exponential 
    * in terms of n for the exact search, where n is the number of vertices.
    * This method implements Algorithm 3 of the paper.
     * @param q the query graph
     * @param db list of graphs in a database
     * @param w graph edit distance threshold for the matches in the search
     * @return graphs all db_i in db s.t. db_i is the same as the query graph q
     * within a graph edit distance threshold w
     */
    public List<Graph> approxFullSearch(Graph q, List<Graph> db, int w) {
        
        /*
        for each graph g ∈ D do 
            if L_m(g,q) > w {  // Lm(g1, g2) is the lower bound of λ(g1, g2)
                continue;
            }
            if τ(g,q) ≤ w {   // τ(g1,g2) is the suboptimal value of λ(g1, g2)
                report g as a result;
                continue;
            }
            if ρ(g,q) ≤ w {  // ρ(g1, g2) is the refined suboptimal value λ(g1, g2)
                report g as a result;
                continue; 
            }
            if λ(g, q) ≤ w { 
                report g as a result;
            }
        */
        
        List<Graph> results = new ArrayList<Graph>();
        
        Graph dbi;
        StarStructure[] sQ = StarStructure.createStarStructureMultiset(q);
        StarStructure[] sg1, sg2;
        int[][] a1, a2;
        double[][] distM;
        double lM, tau, rho, lambda;
        int[] refinedAssign;
        StarStructure s;
        Set<PairInt> e1, e2;
        Graph g;
        boolean swappedSG;
        int i, k, rIdx;
        for (int ii = 0; ii < db.size(); ++ii) {
            dbi = db.get(ii);
            
            sg1 = StarStructure.copy(sQ);
            sg2 = StarStructure.createStarStructureMultiset(dbi);
            
            a1 = createAdjacencyMatrix(sg1);
            a2 = createAdjacencyMatrix(sg2);
            
            // normalize sq1 and sg2 to have same cardinality for bipartite vertex assignments
            
            // order so that sg1.length >= sg2.length
            if (sg1.length < sg2.length) {
                StarStructure[] tmp = sg1;
                sg1 = sg2;
                sg2 = tmp;
                swappedSG = true;
            } else {
                swappedSG = false;
            }
            if (sg1.length > sg2.length) {
                k = sg1.length - sg2.length;
                //insert k vertices to sg2 and set their labels to eps
                StarStructure[] _sg2 = new StarStructure[sg1.length];
                System.arraycopy(sg2, 0, _sg2, 0, sg2.length);
                
                for (i = 0; i < k; ++i) {
                    rIdx = sg1.length + i;
                    s = new StarStructure(rIdx, StarStructure.eps,
                            new int[0], new int[0], new int[0]);
                    _sg2[sg2.length + i] = s;
                }
                sg2 = _sg2;
            }
            assert (sg1.length == sg2.length);
            
            e1 = getEdges(sg1);
            e2 = getEdges(sg2);
            
            // create cost matrix for bipartite assignments of vertexes in sg1 to sg2
            if (this.edgesAreLabeled) {
                distM = StarStructure.createDistanceMatrixV(sg1, sg2);
            } else {
                distM = StarStructure.createDistanceMatrix(sg1, sg2);
            }
            VolgenantJonker vj = new VolgenantJonker();
            double cost = vj.computeAssignment(distM);
            int[] assign = vj.getAssignment();
   
            int mappingDist = mappingDistance(sg1, sg2, assign);
            
            lM = lowerBoundEditDistance(sg1, sg2, mappingDist);
            
            if (lM > w) {
                continue;
            }
            
            if (this.edgesAreLabeled) {
                tau = suboptimalEditDistance(sg1, sg2, e1, e2, assign);
            } else {
                tau = suboptimalEditDistanceV(sg1, sg2, a1, a2, assign);
            }
            if (tau <= w) {
                results.add(dbi);
                continue;
            }
            
            refinedAssign = Arrays.copyOf(assign, assign.length);
            rho = refinedSuboptimalEditDistance(sg1, sg2, e1, e2, a1, a2, refinedAssign, tau, distM);
            if (rho <= w) {
                results.add(dbi);
                continue;
            }
            
            /*
            if λ(g, q) ≤ w { 
                report g as a result;
            }
            */
                        
        } // end loop over db graphs
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    /**
     retrieve all the super-graphs of a given query graph Q from a 
       graph database D = (D_1,D_2,...D_n).
       It finds all the graphs db_i(i=1,2,…,m.lt.n) such that Q is a subgraph 
       of db_i.  It also finds db_i that are subgraphs of Q.
     * @param q the query graph
     * @param db list of graphs in a database
     * @param w graph edit distance threshold
     * @return all graphs db_i in db s.t. db_i is the same as the query graph Q
     */
    public List<Graph> approxSubSearch(Graph q, List<Graph> db, int w) {
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    /**
     * calculate compute τ(g,h) as C(g, h, P)
      <pre>
      references:
      Section 4.3 of Zeng et al. 2009.
      and
      Chapter 3., pg 45, Algorithm 1 of Feng 2017, PHD Thesis in CSE, UNSW, AU,
      "Efficiently Computing Graph Similarity and Graph Connectivity"
     </pre>
     * @param sg1 star structures for graph g1
     * @param sg2 star structures for graph g2
     * @param e1 edges in graph g1. each PairInt is ordered so that u .lt. v.
     * @param e2 edges in graph g2. each PairInt is ordered so that u .lt. v.
     @param assignments array of bipartite assignments.
     * assignments[0] is the matching sg1[0] to sg2[assignments[0]];
     * @return 
     */
    protected double suboptimalEditDistance(StarStructure[] sg1, StarStructure[] sg2,
        Set<PairInt> e1, Set<PairInt> e2, int[] assignments) {
        
        /*
        // this section commented out implements Zeng et al. 2009 which has vertex edits but not edge edits
        int[][] p = createP(assignments);
        
        int[][] c = createLabelMatrix(sg1, sg2, assignments);
        assert(c.length == p.length);
        assert(c[0].length == p[0].length);
        
        //C(g, h, P') = sum_(i|0:n-1) sum_(j|0:n-1) ( c[i][j]*p[i][j]  
        //              + (1/2) || a1 - P*a2*P^T ||_1
        //
        //Assuming that the L1-norm here is the same convention as MatLab:
        //    For p-norm = 1, the L1-norm is the maximum absolute column sum of the matrix.
        //    ||X||_1 = max sum for an arg j where (0<=j<=n-1) sum_(i=0 to n-1) ( |a[i][j] )
        int[][] pA2PT = MatrixUtil.multiply(p,
            MatrixUtil.multiply(a2, MatrixUtil.transpose(p)));
        
        double term2 = 0.5*MatrixUtil.lp1Norm(
            MatrixUtil.elementwiseSubtract(a1, pA2PT));
        
        double term1 = 0;
        int i, j;
        for (i = 0; i < c.length; ++i) {
            for (j = 0; j < c[i].length; ++j) {
                term1 += c[i][j] * p[i][j];
            }
        }
        
        return term1 + term2;
        */
        
        int n = sg1.length;
        int cost = 0;
        
        // vertex relabeling
        int i, j;
        for (i = 0; i < assignments.length; ++i) {
            j = assignments[i];
            if (sg1[i].rootLabel == sg2[j].rootLabel) {
                cost++;
            }
        }
        int i2, j2;
        PairInt edge2;
        // Edge deletion or relabeling
        for (PairInt edge1 : e1) {
            /*
            for each edge (v,v′) in q do
                if edge (f(v), f(v′)) is not in g or l(v,v′) != l(f(v), f(v′)) cost++;
            */
            i = edge1.getX();
            i2 = edge1.getY();
            j = assignments[i];
            j2 = assignments[j];
            if (j < j2) {
                edge2 = new PairInt(j, j2);
            } else {
                edge2 = new PairInt(j2, j);
            }
            int edgeLabel1 = sg1[i].eLabels[sg1[i].reverseOrigVIndexes.get(i2)];
            int edgeLabel2 = sg2[j].eLabels[sg2[j].reverseOrigVIndexes.get(j2)];
            if (!e2.contains(edge2)) {
                cost++;
            } else if (edgeLabel1 != edgeLabel2) {
                cost++;
            }
        }
 //TODO: revisit this.  the only edge inserts are the eps normalization inserts which have no cost
        // Edge insertion
        /*
        for each edge (u,u′) in g
            if edge(f'(u),f'(u′)) is not in q cost++;
        */
        return cost;
    }
    
    /**
     * calculate compute τ(g,h) as C(g, h, P) for the graphs with vertex labeling
     * but no edge labeling
      <pre>
      references:
      Section 4.3 of Zeng et al. 2009.
      and
      Chapter 3., pg 45, Algorithm 1 of Feng 2017, PHD Thesis in CSE, UNSW, AU,
      "Efficiently Computing Graph Similarity and Graph Connectivity"
     </pre>
     * @param sg1 star structures for graph g1
     * @param sg2 star structures for graph g2
     * @param a1 adjacency matrix for graph g1
     * @param a2 adjacency matrix for graph g2
     @param assignments array of bipartite assignments.
     * assignments[0] is the matching sg1[0] to sg2[assignments[0]];
     * @return 
     */
    protected double suboptimalEditDistanceV(StarStructure[] sg1, StarStructure[] sg2,
        int[][] a1, int[][] a2, int[] assignments) {
        
        // this section commented out implements Zeng et al. 2009 which has vertex edits but not edge edits
        int[][] p = createP(assignments);
        
        int[][] c = createLabelMatrix(sg1, sg2, assignments);
        assert(c.length == p.length);
        assert(c[0].length == p[0].length);
        
        //C(g, h, P') = sum_(i|0:n-1) sum_(j|0:n-1) ( c[i][j]*p[i][j]  
        //              + (1/2) || a1 - P*a2*P^T ||_1
        //
        //Assuming that the L1-norm here is the same convention as MatLab:
        //    For p-norm = 1, the L1-norm is the maximum absolute column sum of the matrix.
        //    ||X||_1 = max sum for an arg j where (0<=j<=n-1) sum_(i=0 to n-1) ( |a[i][j] )
        int[][] pA2PT = MatrixUtil.multiply(p,
            MatrixUtil.multiply(a2, MatrixUtil.transpose(p)));
        
        double term2 = 0.5*MatrixUtil.lp1Norm(
            MatrixUtil.elementwiseSubtract(a1, pA2PT));
        
        double term1 = 0;
        int i, j;
        for (i = 0; i < c.length; ++i) {
            for (j = 0; j < c[i].length; ++j) {
                term1 += c[i][j] * p[i][j];
            }
        }
        
        return term1 + term2;
    }
    
    /**
     * calculate Lm(g1, g2) = μ(g1, g2) / max{4, [max{δ(g1), δ(g2)} + 1]}
     * following Sect 4.2.2 of Zeng et al. 2009.
     * @param sg1 star structures for graph g1
     * @param sg2 star structures for graph g2
     * @param mappingDist mapping distance calculated using Definition 4.3 in
     * Section 4.2.1 in Zeng et al. 2009.
     * @return 
     */
    protected double lowerBoundEditDistance(StarStructure[] sg1, StarStructure[] sg2,
        int mappingDist) {
        
        /*
        μ(g1, g2) = mapping distance
        
        δ(g) = max_{v∈V(g)}deg(v)
        
        Lm(g1, g2) = μ(g1, g2) / max{4, [max{δ(g1), δ(g2)} + 1]}
        */
        int maxDeg1 = maxDegree(sg1);
        int maxDeg2 = maxDegree(sg2);
        double denom = Math.max(maxDeg1, maxDeg2) + 1;
        
        double lM = mappingDist/denom;
        
        return lM;
    }
    
    /**
     * Given two multi-sets of star structures S1 and S2,
     * normalize them to create the same cardinality, 
     * and assume P :S1→S2 is a bijection. 
     * The distance ζ between S_1 and S_2 is the summation of the edit distance
     * over an assignment of vertexes solved by bipartite matching of the
     * vertex labels.
     * @param sg1 star structures for graph g1
     * @param sg2 star structures for graph g2
     * @param assignments array of bipartite assignments.
     * assignments[0] is the matching sg1[0] to sg2[assignments[0]];
     * @return 
     */
    protected int mappingDistance(StarStructure[] sg1, StarStructure[] sg2,
        int[] assignments) {
                               
        // usign the assignments, sum the edit distances.
        int sum = 0, i;
        StarStructure s1, s2;
        for (i = 0; i < assignments.length; ++i) {
            s1 = sg1[i];
            s2 = sg2[assignments[i]];
            if (this.edgesAreLabeled) {
                sum += StarStructure.calculateEditDistanceV(s1, s2);
            } else {
                sum += StarStructure.calculateEditDistance(s1, s2);
            }
        }
        
        return sum;
    }

    /**
     * find the maximum degree of a vertex for the graph g.
     * @param sg
     * @return 
     */
    private int maxDegree(StarStructure[] sg) {
        int max = 0, deg;
        for (int i = 0; i < sg.length; ++i) {
            deg = sg[i].vLabels.length;
            if (deg > max) {
                max = deg;
            }
        }
        return max;
    }

    private int[][] createAdjacencyMatrix(StarStructure[] s) {
        int nV = s.length;
        int[][] a = new int[nV][];
        StarStructure si;
        int uIdx, vIdx, j;
        for (int i = 0; i < nV; ++i) {
            a[i] = new int[nV];
            si = s[i];
            uIdx = si.rootIdx;
            for (j = 0; j < si.vLabels.length; ++j) {
                vIdx = si.origVIndexes[j];
                a[uIdx][vIdx] = 1;
            }
        }
        return a;
    }
    
    private Set<PairInt> getEdges(StarStructure[] s) {
        int nV = s.length;
        Set<PairInt> edges = new HashSet<PairInt>();
        StarStructure si;
        int uIdx, vIdx, j;
        for (int i = 0; i < nV; ++i) {
            si = s[i];
            uIdx = si.rootIdx;
            for (j = 0; j < si.vLabels.length; ++j) {
                vIdx = si.origVIndexes[j];
                if (uIdx < vIdx) {
                    edges.add(new PairInt(uIdx, vIdx));
                } else {
                    edges.add(new PairInt(vIdx, uIdx));
                }
            }
        }
        return edges;
    }

    private int[][] createP(int[] assignments) {
        int n = assignments.length;
        int[][] p = new int[n][];
        int i;
        for (i = 0; i < n; ++i) {
            p[i] = new int[n];
        }
        for (i = 0; i < n; ++i) {
            p[i][assignments[i]] = 1;
        }
        assert(MatrixUtil.isAPermutationMatrix(p));
        return p;
    }

    /**
     * create the label matrix C used in Section 3.1 in eqn (1) of 
     * Zeng et al. 2009.
     * C_i,j = 1 if l_g(vi) = l_h(uj)(vi ∈ V(g), uj ∈ V(h)), 
     *       otherwise C_i,j = 0.
     * TODO: revise to include consideration for edge labels.
     * @param sg1
     * @param sg2
     * @param assignments
     * @return 
     */
    private int[][] createLabelMatrix(StarStructure[] sg1, StarStructure[] sg2, 
        int[] assignments) {
        //Ci,j = 1 if l_g(vi) = l_h(uj)(vi ∈ V(g),uj ∈ V(h)), otherwise Ci,j = 0
        int[][] c = new int[sg1.length][];
        int i, j;
        for (i = 0; i < sg1.length; ++i) {
            c[i] = new int[sg2.length];
        }
        for (i = 0; i < assignments.length; ++i) {
            j = assignments[i];
            if (sg1[i].rootLabel == sg2[j].rootLabel) {
                c[i][j] = 1;
            }
        }
        return c;
    }

    /**
     * calculates the refined suboptimal edit distance by swapping the 
     * vertex assignments and calculating the
     * suboptimal distance, keeping the permuted assignments that result in
     * the smallest edit distance.
     * @param sg1 star structures for graph g1
     * @param sg2 star structures for graph g2
     * @param e1
     * @param e2
     * @param a1 adjacency matrix for graph g1
     * @param a2 adjacency matrix for graph g2
     * @param refinedAssign input initial vertex assignments and output
     * refined vertex assignments.
     * @param tau the sub-optimal edit distance C(g,h,P) where
     * P is formed from the given assignments in refinedAssign.
     * @param distM cost matrix for bipartite assignments of vertexes in sg1 to sg2
     * @return
     */
    protected double refinedSuboptimalEditDistance(StarStructure[] sg1, StarStructure[] sg2,
        Set<PairInt> e1, Set<PairInt> e2, int[][] a1, int[][] a2,
        int[] refinedAssign, double tau, double[][] distM) {
        /*
        dist ← C(g,h,P);
        min ← dist;
        for any pair (ui, uj) ∈ V (g) {
            get P′ based on ui and uj; 
            if min > C(g,h,P′) {
                min ← C(g,h,P′);
                Pmin ←P′;
            }
        }
        if min < dist then
            min ← Refine(g, h, Pmin )
        return min;
        */
        
        /*
        How to choose the pair for the statement "for any pair (ui, uj) ∈ V (g)"?
        
        distM holds the cost matrix for bipartite assignments of vertexes in sg1 to sg2.
        
        what swap in assign has best chance of decreasing the graph edit cost?
        one case is an edge in sg1 which is mapped to 2 vertexes in sg2 which 
           are not an edge.
           v1_i in sg1 maps to v2_j in sg2.
           v1_i2 in sg1 maps to v2_j2 in sg2.
           v1_i is adjacent to v1_i2.
           v1_j is not adjacent to v2_j2.
           goal is to find a vertex in sg2 adjacent to v2_j.
           that v2_j_adj vertex in sg2 maps in reverse to v1_i3.
           new assignments:
             V1_i2 <—>  V2_j_adj
             V1_i3 <—> V2_j2
           consistent hanges need to be made to the reverse assign array also.
        */
                 
        int[] assign = Arrays.copyOf(refinedAssign, refinedAssign.length);
        int[] revAssign = reverseAssignment(assign);
        double min = tau;
        
        int iV1, iV2, jV1, jV2, iV3, jV1Adj, jj;
        PairInt pV;
        int[] jV1AdjIdxs;
        for (PairInt edge1 : e1) {
            iV1 = edge1.getX();
            iV2 = edge1.getY();
            jV1 = assign[iV1];
            jV2 = assign[iV2];
            if (jV1 < jV2) {
                pV = new PairInt(jV1, jV2);
            } else {
                pV = new PairInt(jV2, jV1);
            }
            if (e2.contains(pV)) {
                continue;
            }
            jV1AdjIdxs = sg2[jV1].origVIndexes;
            for (jj = 0; jj < jV1AdjIdxs.length; ++jj) {
                jV1Adj = jV1AdjIdxs[jj];
                iV3 = revAssign[jV1Adj];
                
                // tentative changes to assign
                //V1_i2 <—>  V2_j_adj, V1_i3 <—> V2_j2
                assign[iV2] = jV1Adj;
                assign[iV3] = jV2;
                
                if (this.edgesAreLabeled) {
                    tau = suboptimalEditDistance(sg1, sg2, e1, e2, assign);
                } else {
                    tau = suboptimalEditDistanceV(sg1, sg2, a1, a2, assign);
                }
                
                if (tau < min) {
                   min = tau;
                   System.arraycopy(assign, 0, refinedAssign, 0, assign.length);
                   revAssign[jV1Adj] = iV2;
                   revAssign[jV2] = iV3;
                   break;
                } else {
                   //restore assign to latest refineAssign
                   System.arraycopy(refinedAssign, 0, assign, 0, assign.length);
                }
            }
        }
        return tau;
    }

    protected int[] reverseAssignment(int[] assign) {
        int[] r = new int[assign.length];
        for (int i = 0; i < assign.length; ++i) {
            r[assign[i]] = i;
        }
        return r;
    }

    /*
    let C be an n×n label matrix, 
        where Ci,j = 1 if l_g(vi) = l_h(uj) (vi ∈ V(g), uj ∈ V(h)),
         otherwiseCi,j =0. 
    Next, let P be an orthogonal matrix having the property P P T =P T P =I 
    where I is the identity matrix. For each per-
    */
    
    /**
     * an undirected attributed graph
     */
    public static class Graph {
        public final TIntIntMap vLabels;
        public final TObjectIntMap<PairInt> eLabels;
        public final TIntObjectMap<TIntSet> adjMap;
        /**
         * construct a graph instance using given data structures.  Note that this
         * method copies by reference and does not copy by value into new
         * data structures (can use MatrixUtil.copy() methods).
         * @param adjMap
         * @param vLabels
         * @param eLabels 
         */
        public Graph(TIntObjectMap<TIntSet> adjMap, TIntIntMap vLabels, TObjectIntMap<PairInt> eLabels) {
            this.vLabels = vLabels;
            this.eLabels = eLabels;
            this.adjMap = adjMap;
        }        
    }
}
