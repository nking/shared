package algorithms.graphs;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class Complement {
    
    /**
       Given an undirected graph G=(V,E), we define the complement 
       of G as G_bar=(V,E_bar), where E_bar=((u,v): u,v in V, u!=v, and (u,v) 
       is not in E).  In other words, G is the graph containing exactly those 
       edges that are not in G.
       <pre>
       Chap 34 Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms".
       </pre>
       runtime complexity is O(|V|^2).
     @param g an undirected graph
     @return 
     */
    public static Set<PairInt> graphComplement(TIntObjectMap<TIntSet> g) {
        Set<PairInt> edges = extractEdges(g);
        return graphComplement(edges);
    }
    
    /**
     * Given an undirected graph G=(V,E), we define the complement 
       of G as G_bar=(V,E_bar), where E_bar=((u,v): u,v in V, u!=v, and (u,v) is not in E).
       In other words, G is the graph containing exactly those edges that are not in G.
       <pre>
       Chap 34 Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms"
       </pre>
       runtime complexity is O(|V|^2).
     @param edges edges of an undirected graph
     @return 
     */
    public static Set<PairInt> graphComplement(Set<PairInt> edges) {
        TIntSet vertices = extractVertices(edges);
        int[] vs = vertices.toArray();
        // since the algorithm runtime is O(|V|^2), adding an O{|V|lg_2(|V|))
        Arrays.sort(vs);
        
        Set<PairInt> c = new HashSet<PairInt>();
        int i, j, u, v;
        PairInt p;
        for (i = 0; i < vs.length; ++i) {
            u = vs[i];
            for (j = i+1; j < vs.length; ++j) {
                v = vs[j];
                p = new PairInt(u, v);
                if (!edges.contains(p)) {
                    c.add(p);
                }
            }
        }
        return c;
    }
    
    /**
     * given undirected graph g=(V,E) and a subset s of the vertices V,
     * return the complement of s.
     * runtime complexity is O(|V| + |E|)... between best case O(|V| + |E|) and worse case O(|V| * |E|).
     @param g undirected graph g=(V,E)
     @param s subset of vertices in graph g.
     @return the complement of s, that is, V - {s}.
     */
    public static TIntSet setComplement(TIntObjectMap<TIntSet> g, TIntSet s) {
        TIntSet vertices = extractVertices(g);
        vertices.removeAll(s);
        return vertices;
    }

    /**
     *
     @param g
     @return
     */
    public static Set<PairInt> extractEdges(TIntObjectMap<TIntSet> g) {
        TIntObjectIterator<TIntSet> iter = g.iterator();
        Set<PairInt> edges = new HashSet<PairInt>();
        int i, u, v;
        TIntSet set;
        TIntIterator iter2;
        for (i = 0; i < g.size(); ++i) {
            iter.advance();
            u = iter.key();
            set = iter.value();
            iter2 = set.iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                if (v > u) {
                    edges.add(new PairInt(u, v));
                } else {
                    edges.add(new PairInt(v, u));
                }
            }
        }
        return edges;
    }

    /**
     *
     @param edges
     @return
     */
    public static TIntSet extractVertices(Set<PairInt> edges) {
        TIntSet vs = new TIntHashSet();
        for (PairInt uv : edges){
            vs.add(uv.getX());
            vs.add(uv.getY());
        }
        return vs;
    }
    
    /**
     *
     @param g
     @return
     */
    public static TIntSet extractVertices(TIntObjectMap<TIntSet> g) {
        TIntSet vs = new TIntHashSet();
        TIntObjectIterator<TIntSet> iter = g.iterator();
        int i, u, v;
        TIntSet set;
        TIntIterator iter2;
        for (i = 0; i < g.size(); ++i) {
            iter.advance();
            u = iter.key();
            set = iter.value();
            iter2 = set.iterator();
            vs.add(u);
            while (iter2.hasNext()) {
                v = iter2.next();
                vs.add(v);
            }
        }
        return vs;
    }
}
