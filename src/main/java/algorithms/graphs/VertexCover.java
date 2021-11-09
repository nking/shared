package algorithms.graphs;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.List;

/**
 A vertex cover is a subset of a graph's vertices which represents at 
 * least one vertex from every edge in the full graph.
 * 
 * optimal vertex cover is np-hard (non-deterministic polynomial class
   problems at least as hard as the hardest problems in NP.  No known polynomial
   time algorithm, but one can guess a single solution and verify it.)    

 * @author nichole
 */
public class VertexCover {
    
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
     * @return the set cover for undirected graph adjMap, no larger than twice that of optimal.
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
}
