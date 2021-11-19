package algorithms.graphs;

import algorithms.TreeTraversal;
import algorithms.DoublyLinkedList;
import algorithms.NAryTreeNode;
import algorithms.optimization.LinearProgramming;
import algorithms.optimization.LinearProgramming.SlackForm;
import algorithms.optimization.LinearProgramming.StandardForm;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 A vertex cover is a subset of a graph's vertices which represents at 
 least one vertex from every edge in the full graph.
 
 optimal vertex cover is np-hard (non-deterministic polynomial class
   problems at least as hard as the hardest problems in NP.  No known polynomial
   time algorithm, but one can guess a single solution and verify it.)    

 * @author nichole
 */
public class VertexCover {
    
    /**
     * for a tree, find an exact vertex cover which does not include any leaves.
     * runtime complexity is O(|V|).
     <pre>
     implemented from pseudocode in lecture slides of Principal lecturer: Dr Thomas Sauerwald
        Advanced Algorithms, University of Cambridge.
        VII. Approximation Algorithms: Covering Problems
      https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/vertexcover.pdf
      who reference Cormen et al. "Introduction to Algorithms"
     </pre>
     Note that a copy of the tree is made internally (excluding prev and next attributes)
     * and the copied nodes are returned in a Set.  At this time, the identity
     * of each NAryTreeNode in the original root and copied root will not be
     * the same, so only the data properties should be used for identity in the returned set.
     * @param root a n-ary tree root.
     * @return an exact vertex cover of the tree, excluding leaves.
     * Note that a copy of the tree is made internally (excluding prev and next attributes)
     * and the copied nodes are returned in a Set.  At this time, the identity
     * of each NAryTreeNode in the original root and copied root will not be
     * the same, so only the data properties should be used for identity in the returned set.
     */
    public Set<NAryTreeNode> exact(NAryTreeNode root) {
        
        root = NAryTreeNode.copyTree(root);
        
        /*
        There exists an optimal vertex cover which does not include any leaves.
        VERTEX-COVER-TREES(G)
            C=∅ 
            while there exists leaves in G 
                Add all parents of leaves to C
                Remove all leaves and their parents from G 
            return C
        */
                
        TreeTraversal tt = new TreeTraversal();
        DoublyLinkedList<NAryTreeNode> revLevelOrder = tt.getReverseLevelOrderIterative2(root);
        Set<NAryTreeNode> c = new HashSet<NAryTreeNode>();
        if (revLevelOrder.size() < 1) {
            return c;
        }
        Set<NAryTreeNode> children;
        NAryTreeNode node;
        int i = 0;
        while (!revLevelOrder.isEmpty()) {
            node = revLevelOrder.peekFirst();
            //System.out.println("visiting " + node.getData());
            if (node.getParent() != null) {
                if (node.getChildren().isEmpty() && !c.contains(node.getParent())) {
                    // this is a leaf node.  
                    // add it's parent to set C
                    c.add(node.getParent());
                    //System.out.println("  *C.add==" + node.getParent().getData());
                    //  remove node and parent from LinkedList
                    //System.out.println("  unlink parent=" + node.getParent().getData());
                    revLevelOrder.unlink(node.getParent());
                    children = node.getParent().getChildren();
                    // remove children of parent
                    for (NAryTreeNode child : children) {
                        //System.out.println("  unlink child=" + child.getData());
                        if (!c.contains(child)) {
                            revLevelOrder.unlink(child);
                        }
                    }
                    // remove parent from it's parent's children
                    if (node.getParent().getParent() != null) {
                        children = node.getParent().getParent().getChildren();
                        children.remove(node.getParent());
                    }
                } else {
                    //System.out.println("  unlink " + node.getData());
                    revLevelOrder.unlink(node);
                }
            } else if (revLevelOrder.size() == 1) {
                c.add(node);
                //System.out.println("  unlink root " + node.getData());
                revLevelOrder.unlink(node);
            }
        }
        return c;
    }
    
    /**
     * find a vertex cover that is 2-approximate, that is no more than 2 times
     * as large as the optimal vertex cover.
     *
     * A vertex cover is a subset of a graph's vertices which represents at least one vertex 
     * from every edge in the full graph.
     * 
     * The runtime complexity is O(|E| + |V|).

     * <pre>
     * The algorithm implements pseudocode from
     *  from https://www.ics.uci.edu/~goodrich/teach/graph/notes/Approximation.pdf
     * Also see Section 35.1 of Cormen et al. "Introduction to Algorithms".
     * </pre>
     * @param adjMap adjacency map for an undirected graph.
     * @return the vertex cover for undirected graph adjMap, no larger than twice that of optimal.
     */
    public TIntSet approx2(TIntObjectMap<TIntSet> adjMap) {
        
        // make a copy of the graph to edit it
        TIntObjectMap<TIntSet> g = copy(adjMap);
        
        // make a reverse mapping of the graph in order to find edges incident
        //    on to a given vertex.
        TIntObjectMap<TIntSet> reverseMap = reverse(adjMap);
        
        //TIntIntMap c = new TIntIntHashMap();
        TIntSet c = new TIntHashSet();
        
        int u, v, incident;
        TIntSet uEdges, revEdgesU, revEdgesV;
        TIntIterator iter;
        while (!g.isEmpty()) {
            //(1) select an edge e = (u,v) of G
            u = g.keySet().iterator().next();
            uEdges = g.get(u);
            if (uEdges == null || uEdges.isEmpty()) {
                g.remove(u);
                continue;
            }
            v = uEdges.iterator().next();
            
            //(2) add vertices u and v to C
            c.add(u);
            c.add(v);
            
            // remove all edges u->other
            g.remove(u);
            // remove all edges v->other
            g.remove(v);

            // remove all edges other->u
            //   can use reverse map to find those links                        
            revEdgesU = reverseMap.get(u);
            if (revEdgesU != null) {
                iter = revEdgesU.iterator();
                while (iter.hasNext()) {
                    incident = iter.next();
                    if (g.containsKey(incident) && g.get(incident).contains(u)) {
                        g.get(incident).remove(u);
                        if (g.get(incident).isEmpty()) {
                            g.remove(incident);
                        }
                    }
                }
            }
            
            // remove all edges other->v
            //   can use reverse map to find those links
            revEdgesV = reverseMap.get(v);
            if (revEdgesV != null) {
                iter = revEdgesV.iterator();
                while (iter.hasNext()) {
                    incident = iter.next();
                    if (g.containsKey(incident) && g.get(incident).contains(v)) {
                        g.get(incident).remove(v);
                        if (g.get(incident).isEmpty()) {
                            g.remove(incident);
                        }
                    }
                }
            }
        }       
        return c;
    }
    
    /**
     * find a minimum weighted vertex cover for a weighted undirected graph  
     * using a 2-approximate algorithm
     *(though the linear programming implementation used does not have polynomial runtime).
     *
     * A vertex cover is a subset of a graph's vertices which represents at least one vertex 
     * from every edge in the full graph.  The given vertexes have weights
     * associated with them.
     * 
     * <pre>
     * The algorithm implements pseudocode from
     * Section 35.4 of Cormen et al. "Introduction to Algorithms".
     * </pre>
     * @param adjMap adjacency map for an undirected graph.
     * @param weights the weights of each vertex
     * @return the minimum weighted vertex cover for graph G represented by adjMap with vertex weights.
     */
    public TIntSet approx2Weighted(TIntObjectMap<TIntSet> adjMap, double[] weights) {
        
        /*
        for the linear program:
            minimize: 
                summation_v_in_V( w(v)*x(v)
            subject to:
                x(u) + x(v) >= 1 for each (u,v) in E
                x(v) <= 1 for each v in V
            non-negativity constraints:
                x(v) >= 0 for each v in V
        
        for the weighted vc:
            compute optimal x from linear programming.
            C = empty set
            for each v in V
                if (x(v) >= 0.5)
                    add v to C
        */
        
        StandardForm standForm = createLinearProgramInStandardForm(adjMap, weights);
        LinearProgramming lp = new LinearProgramming();
        SlackForm soln = lp.solveUsingSimplexMethod(standForm);
        double[] x = soln.computeBasicSolution();
        
        //System.out.printf("x=%s\n", FormatArray.toString(x, "%.3f"));
        
        TIntSet c = new TIntHashSet();
        int i;
        for (i = 0; i < weights.length; ++i) {
            if (x[i] >= 0.5) {
                c.add(i);
            }
        }
        return c;
    }

    protected TIntObjectMap<TIntSet> copy(TIntObjectMap<TIntSet> adjMap) {
        
        TIntObjectMap<TIntSet> c = new TIntObjectHashMap<TIntSet>();
        TIntObjectIterator<TIntSet> iter = adjMap.iterator();
        TIntSet v;
        int k;
        while (iter.hasNext()) {
            iter.advance();
            k = iter.key();
            v = iter.value();
            c.put(k, new TIntHashSet(v.toArray()));
        }
        return c;
    }

    protected TIntObjectMap<TIntSet> reverse(TIntObjectMap<TIntSet> adjMap) {
        
        TIntObjectMap<TIntSet> r = new TIntObjectHashMap<TIntSet>();
        
        TIntObjectIterator<TIntSet> iter = adjMap.iterator();
        TIntIterator iterV;
        TIntSet vSet, rVSet;
        int u, v;
        while (iter.hasNext()) {
            iter.advance();
            u = iter.key();
            vSet = iter.value();
            
            if (vSet != null && !vSet.isEmpty()) {
                iterV = vSet.iterator();
                while (iterV.hasNext()) {
                    v = iterV.next();
                    rVSet = r.get(v);
                    if (rVSet == null) {
                        rVSet = new TIntHashSet();
                        r.put(v, rVSet);
                    }
                    rVSet.add(u);
                }
            }
        }
        return r;
    }

    protected StandardForm createLinearProgramInStandardForm(
        TIntObjectMap<TIntSet> adjMap, double[] weights) {
        
        /*
         minimize: 
                summation_v_in_V( w(v)*x(v) )
            subject to:
                x(u) + x(v) >= 1 for each (u,v) in E
                x(v) <= 1 for each v in V
            non-negativity constraints:
                x(v) >= 0 for each v in V
        */
        double[] c = Arrays.copyOf(weights, weights.length);
        
        // rows are pairs of u, v, e.g. edges[0][0] = u, edges[0][1] = v
        int[][] edges = extractEdges(adjMap);
        int nE = edges.length;
        int nV = weights.length;
        
        double[][] a = new double[nE + nV][nV];
        double[] b = new double[nE + nV];
        Arrays.fill(b, 1);
        
        int i, u, v;
        for (i = 0; i < nE; ++i) {
            u = edges[i][0];
            v = edges[i][1];
            a[i] = new double[nV];
            a[i][u] = 1;
            a[i][v] = 1;
        }
        int i2;
        for (i = 0, i2=nE; i < nV; ++i, ++i2) {
            a[i2] = new double[nV];
            a[i2][i] = 1;
        }
        
        boolean isMaximization = false;
        int[] constraintComparisons = new int[nE + nV];
        Arrays.fill(constraintComparisons, 0, nE, 1);
        Arrays.fill(constraintComparisons, nE, nV, -1);
        boolean[] nonnegativityConstraints = new boolean[nV];
        Arrays.fill(nonnegativityConstraints, true);
        
        StandardForm standForm = LinearProgramming
            .convertLinearProgramToStandardForm(isMaximization, a, b, c, 
            constraintComparisons, nonnegativityConstraints);
        
        //System.out.printf("graph as Linear Program in standard form=\n%s\n", standForm.toString());

        return standForm;
    }

    /**
     * extract edges from the adjacency map and put them in a 2 dimensional array.
     * @param adjMap adjacency map of an undirected graph.  edge u->v and v->u are
     * considered the same edge.
     * @return two dimensional array of the edges.
     */
    protected int[][] extractEdges(TIntObjectMap<TIntSet> adjMap) {
        Set<PairInt> edges = new HashSet<PairInt>();
        int u, v, i, k, k2;
        PairInt p;
        
        TIntObjectIterator<TIntSet> iter = adjMap.iterator();
        TIntIterator iter2;
        TIntSet adj;
        for (i = 0; i < adjMap.size(); ++i) {
            iter.advance();
            adj = iter.value();
            if (adj == null || adj.isEmpty()) {
                continue;
            }
            k = iter.key();
            iter2 = adj.iterator();
            while (iter2.hasNext()) {
                k2 = iter2.next();
                if (k2 < k) {
                    p = new PairInt(k2, k);
                } else {
                    p = new PairInt(k, k2);
                }
                edges.add(p);
            }
        }
        int[][] out = new int[edges.size()][];
        i = 0;
        for (PairInt p2 : edges) {
            out[i] = new int[]{p2.getX(), p2.getY()};
            i++;
        }
        return out;
    }
}
