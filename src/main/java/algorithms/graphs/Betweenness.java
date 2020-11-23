package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TObjectFloatIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Queue;

/**
 *
 * @author nichole
 */
public class Betweenness {
    
    /**
     * implementation of unweighted graph edge scoring from Girvan-Newman 
     * algorithm, accepting a DAG, but the DAG can 
     * Reference is 2004 Newman and Girvan,
     * "Finding and evaluating community structure in networks".
     * 
     * The runtime complexity is <em>1 + the number of graph roots (i.e. nodes without predecessors)
     *   times  O(|V| + |E|)</em>.
     * 
     * For more information and other graph socring and distance algorithms and 
     * cluster finding (a.k.a. community finding) see also
     * <pre>
     * Chapter 10 from "Mining of Massive Datasets"
       by Leskovec, Rajaraman, and Ullman
     * http://infolab.stanford.edu/~ullman/mmds/ch10n.pdf
     * 
     * and
     * 
     * 2005 paper "Complex networks: Structure and dynamics" by
                Boccalettia, Latorab, Morenod, Chavezf, and  Hwanga
                Physics Reports
               
       and
       
       2010 WWW2010 conference paper "Empirical Comparison of Algorithms for
            Network Community Detection" by Leskovec, Lang, and Mahoney
     
     * </pre>
     * @param adjacencyList
     * @param s
     * @return 
     */
    public Results girvanNewmanDistances(SimpleLinkedListNode[] adjacencyList, final int s) {
        
        //avg O(|E|);  worst: O(|V| + |E|)
        int[] rootIndexes = findRoots(adjacencyList, s);
        
        final int nV = adjacencyList.length;
        // init
        final int[] d = new int[nV];
        final int[] dBest = new int[nV];
        Arrays.fill(dBest, Integer.MAX_VALUE);
        final TIntList[] p = new TIntList[nV];
        final int[] color = new int[nV];
        
        // index is vertex index
        final int[] w = new int[nV];
        final int[] wG = new int[nV];
        
        // key = edge vertexes (u, v) where u is closer to the root, value = edge weight
        final TObjectFloatMap<PairInt> wEdges = new TObjectFloatHashMap<PairInt>();
        
        // key = source vertex of current root, value = hashmap w/ key=vertex, value=distance
        //final TIntObjectMap<TIntIntMap> d = new TIntObjectHashMap<TIntIntMap>();
        
        final TIntSet members = new TIntHashSet();
        
        // runtime is # of roots * O(|V| + |E|)
        
        for (int src : rootIndexes) {
        
            Arrays.fill(w, 0);
            Arrays.fill(color, 0);
            Arrays.fill(d, Integer.MAX_VALUE);
            for (int i = 0; i < nV; ++i) {
                p[i] = new TIntArrayList();
            }
            
            //d.put(src, new TIntIntHashMap());
            //d.get(src).put(src, 0);
            
            color[src] = 1;
            d[src] = 0;
            dBest[src] = 0;
            w[src] = 1;
            
            //System.out.println("root=" + src);
      
            // calc vertex weights
            final TIntList leaf = new TIntArrayList();
            final Queue<Integer> queue = new ArrayDeque<Integer>();
            queue.add(src);
            int u;
            
            while (!queue.isEmpty()) {
                //System.out.println("w=" + Arrays.toString(w));
                //System.out.println("d=" + Arrays.toString(d));
                //System.out.println("dBest=" + Arrays.toString(dBest));
                //System.out.println("color=" + Arrays.toString(color));
                u = queue.remove().intValue();
                members.add(u);
                //System.out.printf("u=%d\n", u);
                SimpleLinkedListNode vNode = adjacencyList[u];
                if (vNode == null || vNode.getKey() == -1) {
                    leaf.add(u);
                    color[u] = 2;
                    //System.out.printf("    LEAF\n");
                    continue;
                }
                while (vNode != null && vNode.getKey() != -1) {
                    int v = vNode.getKey();
                    if (color[v] == 0) {
                        color[v] = 1;
                        d[v] = d[u] + 1;
                        w[v] = w[u];
                        queue.add(v);
                        //System.out.printf("  v=%d\n", v);
                        if (d[v] < dBest[v]) {
                            dBest[v] = d[v];
                        }
                    } else if (d[v] == (d[u] + 1)) {
                        w[v] += w[u];
                    } else {
                        assert (d[v] < (d[u] + 1));
                    }
                    p[v].add(u);
                    vNode = vNode.getNext();
                }
                color[u] = 2;
            }
            assert(queue.isEmpty());
         
            //add weights from previous root traversals
            if (rootIndexes.length > 1) {
                //System.out.println("\n  w=" + Arrays.toString(w));
                for (int i = 0; i < w.length; ++i) {
                    w[i] += wG[i];
                }
                //System.out.println("  wG=" + Arrays.toString(wG));
                //System.out.println("->w=" + Arrays.toString(w));
            }        
            
            // calc edge weights
            TIntSet enqd = new TIntHashSet();
            float e;
            TIntIterator tIter = leaf.iterator();
            int t, i, ip;
            TIntIterator pIter;
            PairInt uv;
            while (tIter.hasNext()) {
                t = tIter.next();
                //System.out.printf("t leaf=%d\n", t);
                pIter = p[t].iterator();
                while (pIter.hasNext()) {
                    i = pIter.next();
                    //System.out.printf("  p=%d\n", i);
                    if (!enqd.contains(i)) {
                        queue.add(i);
                        enqd.add(i);
                    }
                    e = (float) w[i] / (float) w[t];
                    assert(d[i] < d[t]);
                    uv = new PairInt(i, t);
                    wEdges.put(uv, e);
                    //System.out.printf("  edge=(%d, %d) w=%.3e\n", i, t, e);
                }
            }
            //System.out.println("w=" + Arrays.toString(w));
            //System.out.println("d=" + Arrays.toString(d));
            //System.out.println("q="+ queue.toString());
            float e2;
            PairInt ij;
            while (!queue.isEmpty()) {
                i = queue.remove().intValue();
                e = 1;
                //System.out.printf("  e0[%d]=%.3f", i, e);
                SimpleLinkedListNode jNode = adjacencyList[i];
                while (jNode != null && jNode.getKey() != -1) {
                    int j = jNode.getKey();
                    ij = new PairInt(i, j);
                    assert(d[i] < d[j]);
                    assert(wEdges.containsKey(ij));
                    e += wEdges.get(ij);
                    //System.out.printf("  + (w[%d]/w[%d]=%.3f)", i, j, wEdges.get(ij));
                    jNode = jNode.getNext();
                }
                //System.out.printf("  \n   =>%.3f\n", e);

                pIter = p[i].iterator();
                while (pIter.hasNext()) {
                    ip = pIter.next();
                    if (!enqd.contains(ip)) {
                        queue.add(ip);
                        enqd.add(ip);
                    }
                    //System.out.printf("  d[%d]=%d,  d[%d]=%d\n", ip, d[ip], i, d[i]);
                    assert(d[ip] < d[i]);
                    e2 = (float) w[ip] / (float) w[i];
                    e2 *= e;
                    //System.out.printf("    e=(%.3f)*(w[%d]/w[%d]=%.3f)\n", e, ip, i, (float) w[ip] / (float) w[i]);
                    uv = new PairInt(ip, i);
                    //System.out.printf("    edge=(%d, %d) w=%.3e\n", ip, i, e2);
                    wEdges.put(uv, e2);
                }
            }
            
            //update the total tree weights with current
            for (int ii = 0; ii < w.length; ++ii) {
                if (w[ii] > 0) {
                    assert(wG[ii] <= w[ii]);
                    wG[ii] = w[ii];
                }
            }
        }
        
        Results results = new Results();
        results.edges = wEdges;
        results.rootIndexes = rootIndexes;
        results.vertexes = members;
        return results;
    }

    private int[] findRoots(SimpleLinkedListNode[] adjacencyList, int s) {
         
        DFS dfs = new DFS(adjacencyList);
        dfs.walk();
        
        int[] p = dfs.getPredecessorIndexes();
        int nRoots = 0;
        for (int i = 0; i < p.length; ++i) {
            if (p[i] == -1) {
                nRoots++;
            }
        }
        int[] roots = new int[nRoots];
        nRoots = 0;
        for (int i = 0; i < p.length; ++i) {
            if (p[i] == -1) {
                roots[nRoots] = i;
                nRoots++;
            }
        }
        return roots;
    }
    
    public static class Results {
        
        private TObjectFloatMap<PairInt> edges = null;
        
        private int[] rootIndexes = null;
        
        private TIntSet vertexes = null;
       
        /**
         * @return the betweenness scores for the edges
         */
        public TObjectFloatMap<PairInt> getEdges() {
            return edges;
        }

        /**
         * @return the src
         */
        public int[] getRootIndexes() {
            return rootIndexes;
        }

        /**
         * @return the vertexes visited by the algorithm.  any vertex given
         * to algorithm that was not connected to source s is not present in this.
         */
        public TIntSet getVertexes() {
            return vertexes;
        }
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("root indexes=");
            if (rootIndexes != null) {
                sb.append(Arrays.toString(rootIndexes));
            }
            sb.append("\nedges=");
            if (edges != null) {
                TObjectFloatIterator<PairInt> iter = edges.iterator();
                for (int i = 0; i < edges.size(); ++i) {
                    iter.advance();
                    PairInt e = iter.key();
                    float w = iter.value();
                    sb.append(String.format("\n  (%d,%d)=%.3e", e.getX(), e.getY(), w));
                }
            }
            sb.append("\nvertex indexes=");
            if (vertexes != null) {
                sb.append(Arrays.toString(vertexes.toArray()));
            }
            return sb.toString();
        }
    }
}
