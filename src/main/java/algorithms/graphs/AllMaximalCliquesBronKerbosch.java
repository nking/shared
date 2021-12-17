package algorithms.graphs;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Stack;

/**
 * In computer science, the Bron–Kerbosch algorithm is an enumeration algorithm 
 * for finding all maximal cliques in an undirected graph. That is, it lists all 
 * subsets of vertices with the two properties that each pair of vertices in 
 * one of the listed subsets is connected by an edge, and no listed subset can 
 * have any additional vertices added to it while preserving its complete 
 * connectivity. The Bron–Kerbosch algorithm was designed by Dutch scientists 
 * Coenraad Bron and Joep Kerbosch, who published its description in 1973.
 * https://en.m.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
 * 
 * 
 * @author nichole
 */
public class AllMaximalCliquesBronKerbosch {
    
    private final TIntObjectMap<TIntSet> adjMap;
    
    /**
     * 
     * @param adjMap adjacency map for an undirected graph. 
     */
    public AllMaximalCliquesBronKerbosch(TIntObjectMap<TIntSet> adjMap) {
        this.adjMap = adjMap;
    }
    
    /*
    algorithm BronKerbosch2(R, P, X) is
    if P and X are both empty then
        report R as a maximal clique
    choose a pivot vertex u in P ⋃ X
    for each vertex v in P \ N(u) do
        BronKerbosch2(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
        P := P \ {v}
        X := X ⋃ {v}
    
    algorithm BronKerbosch3(G) is
    P = V(G)
    R = X = empty
    for each vertex v in a degeneracy ordering of G do
        BronKerbosch2({v}, P ⋂ N(v), X ⋂ N(v))
        P := P \ {v}
        X := X ⋃ {v}
    */
    
    /**
     * find all maximal cliques in an undirected graph.
     * The worse case runtime complexity is O(3^(n/3)) which is from
     *  Moon & Moser (1965) who found that any n-vertex graph has at most 3^(n/3) maximal cliques.
     * <pre>
     * references https://en.m.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
     * algorithm BronKerbosch3(G)
     * </pre>
     * 
     * @return 
     */
    public List<TIntSet> bronKerBosch() {
        
        /*
        P = V(G)
        R = X = empty
        for each vertex v in a degeneracy ordering of G do
            BronKerbosch2({v}, P ⋂ N(v), X ⋂ N(v))
            P := P \ {v}
            X := X ⋃ {v}
        */
        
        // maximal cliques
        List<TIntSet> mc = new ArrayList<TIntSet>();
                
        TIntSet p = extractVertexes(adjMap);
        TIntSet r = new TIntHashSet();
        TIntSet x = new TIntHashSet();
        
        TIntSet p2, r2, x2;
        
        int[] dOrdering = new int[p.size()];
        DegeneracyOrderingMatulaBeck.findDegeneracyOrder(adjMap, dOrdering);
        for (int v : dOrdering) {

            r2 = new TIntHashSet(1);
            r2.add(v);
            p2 = intersection(p, adjMap.get(v));
            x2 = intersection(x, adjMap.get(v));
                    
            bronKerBosch(r2, p2, x2, mc);
            p.remove(v);
            x.add(v);
        }
        
        return mc;
    }
        
    /**
     * find all maximal cliques in an undirected graph.
     * The worse case runtime complexity is O(3^(n/3)) which is from
     *  Moon & Moser (1965) who found that any n-vertex graph has at most 3^(n/3) maximal cliques.
     * 
     * @param r
     * @param p
     * @param x
     * @param results array to add results to 
     */
    public void bronKerBosch(TIntSet r, TIntSet p, TIntSet x, List<TIntSet> results) {
        if (p.isEmpty() && x.isEmpty()) {
            // r is a maximal clique
            results.add(r);
            return;
        }
        TIntSet px = new TIntHashSet(p);
        px.addAll(x);

        //choose a pivot vertex u in P ⋃ X
        int u = px.iterator().next();
        
        //P \ N(u)
        TIntSet pDiffNU = difference(p, adjMap.get(u));

        TIntIterator iter = pDiffNU.iterator();
        int v;
        TIntSet r2, p2, x2;
        while (iter.hasNext()) {
            v = iter.next();
            r2 = new TIntHashSet(r);
            r2.add(v);
            p2 = intersection(p, adjMap.get(v));
            x2 = intersection(x, adjMap.get(v));
            bronKerBosch(r2, p2, x2, results);
            p.remove(v);
            x.add(v);
        }
    }
    
    protected TIntSet intersection(TIntSet a, TIntSet b) {
        //  (0,1,2,3)  {1,2,4}  C=A-B=0,3  A-C=1,2
        TIntSet c = new TIntHashSet(a);
        TIntSet out = new TIntHashSet(a);
        c.removeAll(b);
        out.removeAll(c);
        return out;
    }
    
    private TIntSet extractVertexes(TIntObjectMap<TIntSet> adjMap) {
        TIntSet vs = new TIntHashSet(adjMap.size());
        TIntObjectIterator<TIntSet> iter = adjMap.iterator();
        TIntIterator iter2;
        while (iter.hasNext()) {
            iter.advance();
            vs.add(iter.key());
            // all vertexes should be a key too, so no need to visit adjacent
            /*iter2 = iter.value().iterator();
            while (iter2.hasNext()) {
                vs.add(iter2.next());
            }*/
        }
        return vs;
    }

    //P \ N(u)
    private TIntSet difference(TIntSet p, TIntSet minus) {
        if (p == null) {
            throw new IllegalArgumentException("p cannot be null");
        }
        TIntSet d = new TIntHashSet(p);
        if (minus != null) {
            d.removeAll(minus);
        }
        return d;
    }

}
