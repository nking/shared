package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.Queue;

/**
 *
 * @author nichole
 */
public class Betweenness {
    
    /**
     * implementation of unweighted graph distance metric of Girvan-Newman 
     * algorithm.
     * Reference is 2004 Newman and Girvan,
     * "Finding and evaluating community structure in networks".
     * 
     * For more information and other graph distance algorithms and 
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
        final int nV = adjacencyList.length;
        // init
        final int[] d = new int[nV];
        final TIntList[] p = new TIntList[nV];
        final int[] color = new int[nV];
        final int[] w = new int[nV];
        for (int i = 0; i < nV; ++i) {
            p[i] = new TIntArrayList();
            d[i] = Integer.MAX_VALUE;
        }
        color[s] = 1;
        d[s] = 0;
        w[s] = 1;
        
        // calc vertex weights
        int nE = 0;
        final TIntList leaf = new TIntArrayList();
        final Queue<Integer> queue = new ArrayDeque<Integer>();
        queue.add(s);
        int u;
        final TIntSet members = new TIntHashSet();
        while (!queue.isEmpty()) {
            u = queue.remove().intValue();
            members.add(u);
            SimpleLinkedListNode vNode = adjacencyList[u];
            if (vNode == null || vNode.getKey() == -1) {
                leaf.add(u);
                color[u] = 2;
                continue;
            }
            while (vNode != null && vNode.getKey() != -1) {
                int v = vNode.getKey();
                if (color[v] == 0) {
                    color[v] = 1;
                    d[v] = d[u] + 1;
                    w[v] = w[u];
                    queue.add(v);
                    nE++;
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
        
        // calc edge weights
        TObjectFloatMap edges = new TObjectFloatHashMap<PairInt>();
        TIntSet enqd = new TIntHashSet();
        float e;
        TIntIterator tIter = leaf.iterator();
        int t, i, ip;
        TIntIterator pIter;
        while (tIter.hasNext()) {
            t = tIter.next();
            pIter = p[t].iterator();
            while (pIter.hasNext()) {
                i = pIter.next();
                if (!enqd.contains(i)) {
                    queue.add(i);
                    enqd.add(i);
                }
                e = (float) w[i] / (float) w[t];
                edges.put(new PairInt(i, t), e);
            }
        }
        float e2;
        while (!queue.isEmpty()) {
            i = queue.remove().intValue();
            e = 1;
            SimpleLinkedListNode jNode = adjacencyList[i];
            while (jNode != null && jNode.getKey() != -1) {
                int j = jNode.getKey();
                e += (float) w[i] / (float) w[j];
                jNode = jNode.getNext();
            }
            
            pIter = p[i].iterator();
            while (pIter.hasNext()) {
                ip = pIter.next();
                if (!enqd.contains(ip)) {
                    queue.add(ip);
                    enqd.add(ip);
                }
                e2 = (float) w[ip] / (float) w[i];
                e2 *= e;
                edges.put(new PairInt(ip, i), e2);
            }
        }
        
        Results results = new Results();
        results.edges = edges;
        results.src = s;
        results.vertexes = members;
        return results;
    }
    
    public static class Results {
        private TObjectFloatMap<PairInt> edges;
        private int src;
        private TIntSet vertexes;

        /**
         * @return the betweenness scores for the edges
         */
        public TObjectFloatMap<PairInt> getEdges() {
            return edges;
        }

        /**
         * @return the src
         */
        public int getSrc() {
            return src;
        }

        /**
         * @return the vertexes visited by the algorithm.  any vertex given
         * to algorithm that was not connected to source s is not present in this.
         */
        public TIntSet getVertexes() {
            return vertexes;
        }
    }
}
