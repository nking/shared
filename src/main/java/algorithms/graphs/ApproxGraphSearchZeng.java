package algorithms.graphs;

import algorithms.PermutationsWithAwait;
import algorithms.matrix.MatrixUtil;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import thirdparty.HungarianAlgorithm;

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
 
 in pattern recognition, GED is called graph matching.
in information theory, GED is seeking the matched configuration of vertices
   that has maximum a posteriori probability w.r.t. the available vertex
   attribute information.
   * 
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
    * in terms of n for the optimal search, where n is the number of vertices.
    * This method implements Algorithm 3 of the paper.
     * @param q the query graph
     * @param db list of graphs in a database
     * @param w graph edit distance threshold for the matches in the search.
     * the lower bound, sub-optimal, refined sub-optimal, and optimal
     * costs are less than or equal to w.
     * @param useAsFilterWithoutOptimal if true, the algorithm will return
     * db graphs that passed the bounds of graph edit distances within the 
     * given threshold w and the algorithm will not execute the exponential
     * optimal algorithm.  if false, the algorithm will run the filters
     * and then the optimal graph search and return the results.
     * @return graphs all db_i in db s.t. db_i is the same as the query graph q
     * within a graph edit distance threshold w
     * @throws java.lang.InterruptedException
     */
    public List<Result> approxFullSearch(Graph q, List<Graph> db, double w,
         boolean useAsFilterWithoutOptimal) throws InterruptedException {
        
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
        
        List<Result> results = new ArrayList<Result>();
        
        Graph dbi;
        StarStructure[] sQ = StarStructure.createStarStructureMultiset(q);
        StarStructure[] sg1, sg2;
        int[][] a1, a2;
        double[][] distM;
        double lM, tau, rho, lambda;
        int[] refinedAssign;
        StarStructure s;
        Graph g;
        boolean swappedSG;
        Norm norm;
        int i, k, rIdx;
        for (int ii = 0; ii < db.size(); ++ii) {
            dbi = db.get(ii);
            
            sg1 = StarStructure.copy(sQ);
            sg2 = StarStructure.createStarStructureMultiset(dbi);
            
            // normalize sq1 and sg2 to have same cardinality for bipartite vertex assignments
            norm = normalize(sg1, sg2);
            sg1 = norm.sg1;
            sg2 = norm.sg2;
            swappedSG = norm.swapped;
            
            a1 = createAdjacencyMatrix(sg1);
            a2 = createAdjacencyMatrix(sg2);
            
            // create cost matrix for bipartite assignments of vertexes in sg1 to sg2
            if (this.edgesAreLabeled) {
                distM = StarStructure.createDistanceMatrix(sg1, sg2);
            } else {
                distM = StarStructure.createDistanceMatrixV(sg1, sg2);
            }
            int[] assign = balancedBipartiteAssignment(distM);
   
            int mappingDist = mappingDistance(sg1, sg2, assign);
            
            lM = lowerBoundEditDistance(sg1, sg2, mappingDist);
            
            if (lM > w) {
                continue;
            }
            
            if (this.edgesAreLabeled) {
                if (swappedSG) {
                    tau = suboptimalEditDistance(sg1, sg2, dbi.eLabels, q.eLabels, assign);
                } else {
                    tau = suboptimalEditDistance(sg1, sg2, q.eLabels, dbi.eLabels, assign);
                }
            } else {
                tau = suboptimalEditDistanceV(sg1, sg2, a1, a2, assign);
            }
            if (tau <= w) {
                results.add(new Result(ii, Result.BOUND.SUBOPTIMAL, assign, tau));
                continue;
            }
            
            refinedAssign = Arrays.copyOf(assign, assign.length);
            if (swappedSG) {
                rho = refinedSuboptimalEditDistance(sg1, sg2, dbi.eLabels, q.eLabels, a1, a2, refinedAssign, tau, distM);
            } else {
                rho = refinedSuboptimalEditDistance(sg1, sg2, q.eLabels, dbi.eLabels, a1, a2, refinedAssign, tau, distM);
            }
            if (rho <= w) {
                results.add(new Result(ii, Result.BOUND.REFINED_SUBOPTIMAL, assign, rho));
                continue;
            }
           
            if (!useAsFilterWithoutOptimal){
                // exponential runtime complexity:
                if (swappedSG) {
                    lambda = optimalEditDistance(sg1, sg2, dbi.eLabels, q.eLabels, a1, a2, refinedAssign, tau);
                } else {
                    lambda = optimalEditDistance(sg1, sg2, q.eLabels, dbi.eLabels, a1, a2, refinedAssign, tau);
                }
                if (lambda <= w) {
                    results.add(new Result(ii, Result.BOUND.OPTIMAL, refinedAssign, lambda));
                }
            } else {
                results.add(new Result(ii, Result.BOUND.LOWER, assign, lM));
            }
        } // end loop over db graphs
        
        return results;
    }
    
    /**
     * perform a lower-bound filter for sub-graph search, then an approx similarity
     * sub-graph search if requested.
     * sub-graph search retrieves all the super-graphs of a given query graph Q from a 
       graph database db = (db_0, db_1, ...db_n-1).
       It finds all the graphs db_i(i=0,1,2,…,m-1.lt.n-1) such that Q is a subgraph 
       of db_i.  It also finds db_i that are subgraphs of Q.
     * The method implements the description of the APPSUB algorithm from 
     * Zeng et al. 2009.
      <pre> 
      Comparing Stars: On Approximating Graph Edit Distance
       Zeng, Tung, Wang, Feng, & Zhou 2009, 
       Proceedings of the VLDB Endowment, Volume 2, Issue 1,
       August 2009 pp 25–36
       https://doi.org/10.14778/1687627.1687631
    </pre>
     from https://en.wikipedia.org/wiki/Subgraph_isomorphism_problem
     Subgraph isomorphism is a generalization of both the maximum clique problem and the problem of testing whether a graph contains a Hamiltonian cycle, and is therefore NP-complete. However certain other cases of subgraph isomorphism may be solved in polynomial time.
     still browsing these for the sub-graph search after the lower bound filter.
     * <pre>
     * 
     * inexact graph matching method with subgraph indexing:
     * T. E. Choe, H. Deng, F. Guo, M. W. Lee and N. Haering, 
     * "Semantic Video-to-Video Search Using Sub-graph Grouping and Matching," 
     * 2013 IEEE International Conference on Computer Vision Workshops, 
     * 2013, pp. 787-794, doi: 10.1109/ICCVW.2013.108.
     * 
     * Graph Pattern Matching: From Intractable to Polynomial Time, 
     * Fan et al. 2010
     * The 36th International Conference on Very Large Data Bases, September 13-17, 
     * 2010, Singapore.  Proceedings of the VLDB Endowment, Vol. 3, No. 1
     * In a variety of applications one wants to inspect the connectivity of a 
     * pair of nodes via a path of an arbitrary length [10, 16, 29] or with a 
     * bound on the number of hops.
     *      * 
     * Tu et al. 2020, "Inexact Attributed Subgraph Matching"
     * 
     * https://github.com/lihuiliullh/GFinder-Proj
     * 
     * Fard, A., Nisar, M.U., Ramaswamy, L., Miller, J.A., Saltz, M.: A distributed vertex-centric approach for pattern matching in massive graphs. In: BigData, pp 403–411 (2013)

       Foggia, P., Percannella, G., Vento, M.: Graph matching and learning in pattern recognition in the last 10 years. Int. J. Pattern Recognit. Artif. Intell. 28(1), 1450001 (2014)

       Gallagher, B.: Matching structure and semantics: A survey on graph-based pattern matching. In: AAAI FS-06-02, pp 45–53 (2006)

       Khan, A., Wu, Y., Aggarwal, C.C., Yan, X.: NeMa: Fast graph search with label similarity. PVLDB
6(3), 181–192 (2013)
       
       Liu,G.,Zheng,K.,Wang,Y.,Orgun,M.A.,Liu,A.,Zhao,L.,Zhou,X.:Multi-constrainedgraphpattern matching in large-scale contextual social graphs. In: ICDE, pp 351–362 (2015)
       
       Ma, S., Cao, Y., Fan, W., Huai, J., Wo, T.: Capturing topology in graph pattern matching. PVLDB 5(4),
310–321 (2011)
* 
*      Tabei, Y., Tsuda, K.: Kernel-based similarity search in massive graph databases with wavelet trees. In:
SDM, pp 154–163 (2011)
* 
*      Wang, X., Smalter, A.M., Huan, J., Lushington, G.H.: G-hash: towards fast kernel-based similarity search in large graph databases. In: EDBT, pp 472–480 (2009)
* 
*      Yan, X., Yu, P.S., Han, J.: Substructure similarity search in graph databases. In: SIGMOD, pp 766–777
(2005)
* 
*      Zhang, S., Yang, J., Jin, W.: SAPPER: Subgraph indexing and approximate matching in large graphs. PVLDB 3(1), 1185–1194 (2010
* 
*      Zou,L.,Chen,L.,O ̈zsu,M.T.:Distance-join:Patternmatchqueryinalargegraphdatabase.PVLDB 2(1), 886–897 (2009)
* 
     * </pre>
  
     * @param q the query graph
     * @param db list of graphs in a database
     * @param w graph edit distance threshold used in the search when comparing
     * the query with each database graph.
     * the total threshold for the lower bound is abs(|V_q|-|V_db_i|) + abs(|E_q|-|E_db_i|) + 2*w.
     * @param useAsFilterWithoutOptimal if true, the algorithm will return
     * db graphs that passed the bounds of graph edit distances within the 
     * given threshold w and the algorithm will not execute the exponential
     * optimal algorithm.  if false, the algorithm will run the filters
     * and then the optimal graph search and return the results.
     * @return all graphs db_i in db s.t. db_i is the same as the query graph Q
     */
    public List<Result> approxSubSearch(Graph q, List<Graph> db, double w,
        boolean useAsFilterWithoutOptimal) {
        
        /*
        AppSub inherently supports both two kinds of subgraph search, i.e., 
        traditional subgraph search[37] and containment search[9].
        
        Lemma 2.2 g1 is subgraph isomorphic to g2 iff λ(g1, g2) = (|E2| − |E1|) + (|V2| − |V1|).
        
        λ′(s1, s2 ) = T′(s1, s2 ) + d(L1, L2) (+ d(L1E´, L2E) if there are edge labels are in the graph model)
           where T′(s1,s2)= {2+|L1|+|L2| if l(r1)̸=l(r2) }
                            {0, otherwise               }
           This is implemented in StarStructure.calculateEditDistanceNoRelabeling()
        
        A graph g1 is said to be θ-subgraph isomorphic to g2 if there exists a 
            graph g3 s.t. g3 ⊑ g2 and λ′(g1, g3) ≤ θ.
        
        if g1 is a θ-subgraph of g2, λ′(g1, g2) ≤ L + 2θ
            where L = |E2| − |E1| + |V2| − |V1|,
        
        Therefore if L′_m(g1, g2) > L + 2θ, g2 can be safely filtered.        
        */
        
        List<Result> results = new ArrayList<Result>();
        
        Graph dbi;
        StarStructure[] sQ = StarStructure.createStarStructureMultiset(q);
        StarStructure[] sg1, sg2;
        double[][] distM;
        double lM, w2;
        Norm norm;
        boolean swappedSG;
        for (int ii = 0; ii < db.size(); ++ii) {
            dbi = db.get(ii);
                        
            sg1 = StarStructure.copy(sQ);
            sg2 = StarStructure.createStarStructureMultiset(dbi);
            
            // normalize sq1 and sg2 to have same cardinality for bipartite vertex assignments
            norm = normalize(sg1, sg2);
            sg1 = norm.sg1;
            sg2 = norm.sg2;
            swappedSG = norm.swapped;
                
            // create cost matrix for bipartite assignments of vertexes in sg1 to sg2
            if (this.edgesAreLabeled) {
                distM = StarStructure.createDistanceMatrixNoRelabeling(sg1, sg2);
            } else {
                distM = StarStructure.createDistanceMatrixNoRelabelingV(sg1, sg2);
            }
            int[] assign = balancedBipartiteAssignment(distM);
   
            int mappingDist = mappingDistance(sg1, sg2, assign);
            
            lM = lowerBoundEditDistance(sg1, sg2, mappingDist);
            
            int l = Math.abs(q.eLabels.size() - dbi.eLabels.size()) 
                + Math.abs(q.vLabels.size() - dbi.vLabels.size());
            
            w2 = l + 2*w;
            System.out.printf("i=%d lM=%.3f w2=%.3f%n", ii, lM, w2);
            
            //L′_m(g1, g2) > L + 2θ can filter out g2
            if (lM > w2) {
                continue;
            }
            
            if (useAsFilterWithoutOptimal) {
                results.add(new Result(ii, Result.BOUND.LOWER, assign, lM));
                continue;
            }
            
            // an exponential exact solution would be to use subsets of the
            //   graph dbi of size q then solve using full graph (which uses filters and then optimal solution which is all permutations).
            throw new UnsupportedOperationException("the filter is implemented, "
                + "but the approx subgraph search is not yet implemented.  "
                + " meanwhile, considering the following implemented java projects: "
                + " https://github.com/jgrapht/jgrapht/blob/master/docs/guide-templates/VertexAndEdgeTypes.md"
                + " https://github.com/pmoris/miles-subgraph-miner"
                + " There are C++ approx subgraph search implementations for large graphs also.");
            
        } // end loop over db graphs
        
        return results;        
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
     * @param e1Labels map having key = edge in graph g1, value = edge label
     * @param e2Labels map having key = edge in graph g2, value = edge label
     @param assignments array of bipartite assignments.
     * assignments[0] is the matching sg1[0] to sg2[assignments[0]];
     * @return 
     */
    protected double suboptimalEditDistance(StarStructure[] sg1, StarStructure[] sg2,
        TObjectIntMap<PairInt> e1Labels, TObjectIntMap<PairInt> e2Labels,
        int[] assignments) {
        
        /*
        Efficiently Computing Graph Similarity and Graph Connectivity
        by Xing Feng, 2017, PhD Thesis, CSE, UNSW, AU
        Chapter 3. Graph Edit Distance Computation, pg 45
        Algorithm 1: EditorialCost
        Input: Graphs q and g, and a mapping f from V(q) to V(g)
        Output: Editorial cost δ_f(q, g)
        int cost = 0;
        // vertex relabeling
        for each vertex v in q do
            if l(v)!=l(f(v)) cost++;
        // Edge deletion or relabeling
        for each edge (v,v′) in q do
            if edge (f(v), f(v′)) is not in g or l(v,v′) != l(f(v), f(v′)) cost++;
        // Edge insertion
        for each edge (u,u′) in g
            if edge(f'(u),f'(u′)) is not in q cost++;
        */
        
        int n = sg1.length;
        int costVSubst = 0;
        int costEDelSubst = 0;
        int costEIns = 0;
                
        // vertex relabeling
        int i, j;
        for (i = 0; i < assignments.length; ++i) {
            j = assignments[i];
            if (sg1[i].rootLabel != sg2[j].rootLabel) {
                costVSubst++;
            }
        }
        int i2, j2, i1r, j1r;
        PairInt edge2;
        // Edge deletion or relabeling
        int e2Label;
        for (PairInt edge1 : e1Labels.keySet()) {
            /*
            for each edge (v,v′) in q do
                if edge (f(v), f(v′)) is not in g or l(v,v′) != l(f(v), f(v′)) cost++;
            */
            i = edge1.getX();
            i2 = edge1.getY();
            j = assignments[i];
            j2 = assignments[i2];
            edge2 = new PairInt(j, j2);
            if (e2Labels.containsKey(edge2)) {
                e2Label = e2Labels.get(edge2);
            } else if (e2Labels.containsKey(new PairInt(j2, j))) {
                edge2 = new PairInt(j2, j);
                e2Label = e2Labels.get(edge2);
            } else {
                costEDelSubst++;
                continue;
            }
            if (e1Labels.get(edge1) != e2Label) {
                costEDelSubst++;
            }
        }
        
        TIntIntMap revAssign = reverseAssignment(assignments);

        // Edge insertion
        //for each edge (u,u′) in g
        //    if edge(f'(u),f'(u′)) is not in q cost++;
        TObjectIntIterator<PairInt> iter = e2Labels.iterator();
        PairInt edge1;
        int e1Label;
        while (iter.hasNext()) {
            iter.advance();
            edge2 = iter.key();
            e2Label = iter.value();
            j = edge2.getX();
            j2 = edge2.getY();
            i = revAssign.get(j);
            i2 = revAssign.get(j2);
            
            edge1 = new PairInt(i, i2);
            if (e1Labels.containsKey(edge1)) {
                e1Label = e1Labels.get(edge1);
            } else if (e1Labels.containsKey(new PairInt(i2, i))) {
                edge1 = new PairInt(i2, i);
                e1Label = e1Labels.get(edge1);
            } else {
                costEIns++;
                continue;
            }
            if (e1Label != e2Label) {
                costEIns++;
            }
        }
        
        //int cost = costVSubst + costEDelSubst + costEIns;
        // or use the maximum for edges to not count them twice
        int cost = costVSubst + Math.max(costEDelSubst, costEIns);
        
        //System.out.printf("costs vSubst=%d, eDelSubst=%d, eIns=%d ==> c=%d%n",
        //    costVSubst, costEDelSubst, costEIns, cost);
        
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
        
        // Zeng et al. 2009 which has vertex edits but not edge edits
        int[][] p = MatrixUtil.createPermutationMatrix(assignments);
        
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
            MatrixUtil.pointwiseSubtract(a1, pA2PT));
        
        //this is a large term, summed from labels being equal.
        //origins are  BLP paper  by Justice & Hero, eqn(18).
        double term1 = 0;
        int i, j;
        for (i = 0; i < c.length; ++i) {
            for (j = 0; j < c[i].length; ++j) {
                term1 += c[i][j] * p[i][j];
            }
        }
        
        //System.out.printf("term1=%.2f term2=%.2f%n", term1, term2);
        
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
        int sum = 0, i, d;
        StarStructure s1, s2;
        for (i = 0; i < assignments.length; ++i) {
            s1 = sg1[i];
            s2 = sg2[assignments[i]];
            if (this.edgesAreLabeled) {
                d = StarStructure.calculateEditDistance(s1, s2);
            } else {
                d = StarStructure.calculateEditDistanceV(s1, s2);
            }
            sum += d;
            //System.out.printf("i=%d ed=%d sum=%d%n", i, d, sum);
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

    static int[][] createAdjacencyMatrix(StarStructure[] s) {
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
    
    static Set<PairInt> getEdges(StarStructure[] s) {
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
            //NLK, changing the fill values from 0 to 1 because the non-matches will cost 1
            Arrays.fill(c[i], 1);
        }
        for (i = 0; i < assignments.length; ++i) {
            j = assignments[i];
            if (sg1[i].rootLabel == sg2[j].rootLabel) {
                //NLK, changing to not penalize if vertex labels match
                //c[i][j] = 1;
                c[i][j] = 0;
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
     * @param e1Labels map having key = edge in graph g1, value = edge label
     * @param e2Labels map having key = edge in graph g2, value = edge label
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
        TObjectIntMap<PairInt> e1Labels, TObjectIntMap<PairInt> e2Labels, int[][] a1, int[][] a2,
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
        TIntIntMap revAssign = reverseAssignment(assign);
        double min = tau;
        
        int iV1, iV2, jV1, jV2, iV3, jV1Adj, jj;
        int[] jV1AdjIdxs;
        for (PairInt edge1 : e1Labels.keySet()) {
            iV1 = edge1.getX();
            iV2 = edge1.getY();
            jV1 = assign[iV1];
            jV2 = assign[iV2];
            
            if (e2Labels.containsKey(new PairInt(jV1, jV2)) ||
                e2Labels.containsKey(new PairInt(jV2, jV1))) {
                continue;
            }
            
            jV1AdjIdxs = sg2[jV1].origVIndexes;
            for (jj = 0; jj < jV1AdjIdxs.length; ++jj) {
                jV1Adj = jV1AdjIdxs[jj];
                iV3 = revAssign.get(jV1Adj);
                
                // tentative changes to assign
                //V1_i2 <—>  V2_j_adj, V1_i3 <—> V2_j2
                assign[iV2] = jV1Adj;
                assign[iV3] = jV2;
                
                if (this.edgesAreLabeled) {
                    tau = suboptimalEditDistance(sg1, sg2, e1Labels, e2Labels, assign);
                } else {
                    tau = suboptimalEditDistanceV(sg1, sg2, a1, a2, assign);
                }
                
                if (tau < min) {
                   min = tau;
                   System.arraycopy(assign, 0, refinedAssign, 0, assign.length);
                   revAssign.put(jV1Adj, iV2);
                   revAssign.put(jV2, iV3);
                   break;
                } else {
                   //restore assign to latest refineAssign
                   System.arraycopy(refinedAssign, 0, assign, 0, assign.length);
                }
            }
        }
        return min;
    }
    
    /**
     * calculates the optimal edit distance by permuting the vertex assignments
     * to find the one which results in the minimum C(g, h, P).  
     * (see section 4.3 of Zeng et al. 2009).
     * the runtime complexity is n!.
     * @param sg1 star structures for graph g1
     * @param sg2 star structures for graph g2
     * @param e1Labels map having key = edge in graph g1, value = edge label
     * @param e2Labels map having key = edge in graph g2, value = edge label
     * @param a1 adjacency matrix for graph g1
     * @param a2 adjacency matrix for graph g2
     * @param refinedAssign input initial vertex assignments and output
     * refined vertex assignments.
     * @param tau the sub-optimal edit distance C(g,h,P) where
     * P is formed from the given assignments in refinedAssign.
     * @return
     * @throws java.lang.InterruptedException exception thrown if thread is 
     * interrupted.  The permutation code is running in a separate thread using
     * a semaphore model to pause and continue execution.
     */
    protected double optimalEditDistance(StarStructure[] sg1, StarStructure[] sg2,
        TObjectIntMap<PairInt> e1Labels, TObjectIntMap<PairInt> e2Labels, int[][] a1, int[][] a2,
        int[] refinedAssign, double tau) throws InterruptedException {
                
        int[] assign = new int[refinedAssign.length];
        double min = tau;
        
        PermutationsWithAwait perm = new PermutationsWithAwait(Arrays.copyOf(refinedAssign, refinedAssign.length));
        
        while (perm.hasNext()) {
            
            perm.getNext(assign);
                        
            if (this.edgesAreLabeled) {
                tau = suboptimalEditDistance(sg1, sg2, e1Labels, e2Labels, assign);
            } else {
                tau = suboptimalEditDistanceV(sg1, sg2, a1, a2, assign);
            }
            
            if (tau < min) {
                min = tau;
                System.arraycopy(assign, 0, refinedAssign, 0, assign.length);
            }            
        }
        
        return min;
    }

    protected TIntIntMap reverseAssignment(int[] assign) {
        TIntIntMap r = new TIntIntHashMap();
        for (int i = 0; i < assign.length; ++i) {
            r.put(assign[i], i);
        }
        return r;
    }

    static Norm normalize(StarStructure[] sg1, StarStructure[] sg2) {
        Norm norm = new Norm();
        // order so that sg1.length >= sg2.length
        if (sg1.length < sg2.length) {
            StarStructure[] tmp = sg1;
            sg1 = sg2;
            sg2 = tmp;
            norm.swapped = true;
        } else {
            norm.swapped = false;
        }
        StarStructure s;
        int k, i, rIdx;
        if (sg1.length > sg2.length) {
            k = sg1.length - sg2.length;
            //insert k vertices to sg2 and set their labels to eps
            StarStructure[] _sg2 = new StarStructure[sg1.length];
            System.arraycopy(sg2, 0, _sg2, 0, sg2.length);

            for (i = 0; i < k; ++i) {
                rIdx = sg2.length + i;
                s = new StarStructure(rIdx, StarStructure.eps,
                    new int[0], new int[0], new int[0]);
                _sg2[rIdx] = s;
            }
            sg2 = _sg2;
        }
        assert (sg1.length == sg2.length);
        norm.sg1 = sg1;
        norm.sg2 = sg2;
        return norm;
    }

    static int[] balancedBipartiteAssignment(double[][] distM) {
        HungarianAlgorithm ha = new HungarianAlgorithm();
        int[][] hAssign = ha.computeAssignments(MatrixUtil.convertToFloat(distM));
        assert(hAssign.length == distM.length);
        int[] assign = new int[distM.length];
        int i1, i2, i; 
        for (i = 0; i < hAssign.length; ++i) {
            i1 = hAssign[i][0];
            i2 = hAssign[i][1];
            assign[i1] = i2;
        }
        return assign;
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
        
        private static TObjectIntMap<PairInt> copy(TObjectIntMap<PairInt> a) {
            TObjectIntMap<PairInt> c = new TObjectIntHashMap<PairInt>();

            TObjectIntIterator<PairInt> iter = a.iterator();
            PairInt p;
            int v;
            while (iter.hasNext()) {
                iter.advance();
                p = iter.key();
                v = iter.value();
                c.put(p.copy(), v);
            }

            return c;
        }
        
        public static Graph copy(Graph g) {
            return new Graph(Graph.copy(g.adjMap), Graph.copy(g.vLabels),
                Graph.copy(g.eLabels));
        }
        
        private static TIntObjectMap<TIntSet> copy(TIntObjectMap<TIntSet> a) {
            TIntObjectMap<TIntSet> c = new TIntObjectHashMap<TIntSet>();

            TIntObjectIterator<TIntSet> iter = a.iterator();
            int k;
            TIntSet v;
            while (iter.hasNext()) {
                iter.advance();
                k = iter.key();
                v = iter.value();
                c.put(k, new TIntHashSet(v));
            }

            return c;
        }

        private static TIntIntMap copy(TIntIntMap a) {
            TIntIntMap c = new TIntIntHashMap();

            TIntIntIterator iter = a.iterator();
            while (iter.hasNext()) {
                iter.advance();
                c.put(iter.key(), iter.value());
            }
            return c;
        }
    }
    
    public static class Result {
        
        public static enum BOUND {
            LOWER, SUBOPTIMAL, REFINED_SUBOPTIMAL, OPTIMAL
        }
        
        /**
         * the graph node assignment w.r.t. the query graph
         */
        public final int[] assignment;
        
        public final int dbGraphIndex;
        
        public final double editCost;
        
        public final BOUND bound;
        
        public Result(int dbGraphIndex, BOUND bound, int[] assign, double editCost) {
            this.bound = bound;
            this.dbGraphIndex = dbGraphIndex;
            this.assignment = Arrays.copyOf(assign, assign.length);
            this.editCost = editCost;
        }
        
    }

    static class Norm {
        StarStructure[] sg1;
        StarStructure[] sg2;
        boolean swapped;
    }
}
