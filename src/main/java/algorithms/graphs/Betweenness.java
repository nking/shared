package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TObjectFloatIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Queue;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class Betweenness {
    
    private Logger log = Logger.getLogger(getClass().getSimpleName());
    private Level logLevel = Level.FINE;
    /**
     * implementation of unweighted graph edge scoring from Girvan-Newman 
     * algorithm, accepting a DAG.   Each graph node without a predecessor
     * is iterated over to calculated node and edge weights.  
     * The choice of the src node as a method argument does not affect the 
     * final scored edge weights,
       but is present in case information about the graph provides a best src
       node to start with (e.g. when the graph only has one parent node).
     
     * Reference is 2004 Newman and Girvan,
     * "Finding and evaluating community structure in networks".
     * 
     * The runtime complexity is <em>1 + the number of graph roots (i.e. nodes without predecessors)
     *   times  O(|V| + |E|)</em>.
     * 
     * For more information and other graph scoring and distance algorithms and 
     * cluster finding (a.k.a. community finding) see also
      <pre>
       Chapter 10 from "Mining of Massive Datasets"
       by Leskovec, Rajaraman, and Ullman
       http://infolab.stanford.edu/~ullman/mmds/ch10n.pdf
       
       and
       
       2005 paper "Complex networks: Structure and dynamics" by
                Boccalettia, Latorab, Morenod, Chavezf, and  Hwanga
                Physics Reports
               
       and
       
       2010 WWW2010 conference paper "Empirical Comparison of Algorithms for
            Network Community Detection" by Leskovec, Lang, and Mahoney
     
     </pre>
     * @param adjacencyList an adjacency list for an unweighted DAG.
     * @param s source node to begin the graph traversal
     * @return scores for the edges of the graph, calculated using the Girvan-Newman algorithm.
     */
    public Results girvanNewmanDAG(SimpleLinkedListNode[] adjacencyList, final int s) {
        
        //avg O(|E|);  worst: O(|V| + |E|)
        int[] rootIndexes = findRoots(adjacencyList, s);
        
        if (rootIndexes.length == 0) {
            throw new IllegalArgumentException("this algorithm operates on a " 
            + " DAG and has been given a graph without a true root node.  " +
            " If have a directionless adjacency list (i.e. u:v and v:u exist)" +
            " then use girvanNewmanDirectionless with a choise of a source node");
        }
        
        //print(adjacencyList);
        //System.out.println("roots=" + Arrays.toString(rootIndexes));
        
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
        
        final TIntSet members = new TIntHashSet();
        
        // runtime is # of roots * O(|V| + |E|)
        
        for (int src : rootIndexes) {
        
            Arrays.fill(w, 0);
            Arrays.fill(color, 0);
            Arrays.fill(d, Integer.MAX_VALUE);
            for (int i = 0; i < nV; ++i) {
                p[i] = new TIntArrayList();
            }
            
            color[src] = 1;
            d[src] = 0;
            dBest[src] = 0;
            w[src] = 1;
            
            log.log(logLevel, "root=" + src);
      
            // calc vertex weights
            final TIntList leaf = new TIntArrayList();
            final Queue<Integer> queue = new ArrayDeque<Integer>();
            queue.add(src);
            int u;
            
            while (!queue.isEmpty()) {
            //    log.log(logLevel, "w=" + Arrays.toString(w));
            //    log.log(logLevel, "d=" + Arrays.toString(d));
            //    log.log(logLevel, "dBest=" + Arrays.toString(dBest));
            //    log.log(logLevel, "color=" + Arrays.toString(color));
                u = queue.remove().intValue();
                members.add(u);
                log.log(logLevel, String.format("u=%d\n", u));
                SimpleLinkedListNode vNode = adjacencyList[u];
                if (vNode == null || vNode.getKey() == -1) {
                    leaf.add(u);
                    color[u] = 2;
                    log.log(logLevel, "    LEAF\n");
                    continue;
                }
                while (vNode != null && vNode.getKey() != -1) {
                    int v = vNode.getKey();
                    if (color[v] == 0) {
                        color[v] = 1;
                        d[v] = d[u] + 1;
                        w[v] = w[u];
                        queue.add(v);
                        log.log(logLevel, String.format("  v=%d\n", v));
                        if (d[v] < dBest[v]) {
                            dBest[v] = d[v];
                        }
                    } else if (d[v] == (d[u] + 1)) {
                        w[v] += w[u];
                    }
                    if (d[u] < d[v]) {
                        p[v].add(u);
                        log.log(logLevel, String.format("  p[%d]=%d\n", v, u));
                    }
                    vNode = vNode.getNext();
                }
                color[u] = 2;
            }
            assert(queue.isEmpty());
         
            //add weights from previous root traversals
            if (rootIndexes.length > 1) {
           //     log.log(logLevel, "\n  w=" + Arrays.toString(w));
                for (int i = 0; i < w.length; ++i) {
                    w[i] += wG[i];
                }
           //     log.log(logLevel, "  wG=" + Arrays.toString(wG));
           //     log.log(logLevel, "->w=" + Arrays.toString(w));
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
                log.log(logLevel, String.format("t leaf=%d\n", t));
                pIter = p[t].iterator();
                while (pIter.hasNext()) {
                    i = pIter.next();
                    log.log(logLevel, String.format("  p=%d\n", i));
                    if (!enqd.contains(i)) {
                        queue.add(i);
                        enqd.add(i);
                    }
                    /*
                    e = (float) w[i] / (float) w[t];
                    assert(d[i] < d[t]);
                    uv = new PairInt(i, t);
                    wEdges.put(uv, e);
                    log.log(logLevel, String.format("  edge=(%d, %d) w=%.3e\n", i, t, e));
                    */
                    // calc edges for all direct children of i instead of only child t
                    SimpleLinkedListNode jNode = adjacencyList[i];
                    while (jNode != null && jNode.getKey() != -1) {
                        int j = jNode.getKey();
                        e = (float) w[i] / (float) w[j];
                        assert(d[i] < d[j]);
                        uv = new PairInt(i, j);
                        wEdges.put(uv, e);
                        log.log(logLevel, String.format("  edge=(%d, %d) w=%.3e\n", i, j, e));
                        jNode = jNode.getNext();
                    }
                }
            }
         //   log.log(logLevel, "w={0}", Arrays.toString(w));
         //   log.log(logLevel, "d=" + Arrays.toString(d));
            log.log(logLevel, "q="+ queue.toString());
            float e2;
            PairInt ij;
            while (!queue.isEmpty()) {
                i = queue.remove();
                e = 1;
                log.log(logLevel, String.format("  e0[%d]=%.3f", i, e));
                SimpleLinkedListNode jNode = adjacencyList[i];
                while (jNode != null && jNode.getKey() != -1) {
                    int j = jNode.getKey();
                    ij = new PairInt(i, j);
                    log.log(logLevel, String.format("  i=%d j=%d", i, j));
///ERROR: look for why no edge from 8->2 in graph karate.gml
                    assert(d[i] < d[j]);
                    assert(wEdges.containsKey(ij));
                    e += wEdges.get(ij);
                    log.log(logLevel, String.format("  + (w[%d]/w[%d]=%.3f)", i, j, wEdges.get(ij)));
                    jNode = jNode.getNext();
                }
                log.log(logLevel, String.format("  \n   =>%.3f\n", e));

                pIter = p[i].iterator();
                while (pIter.hasNext()) {
                    ip = pIter.next();
                    if (!enqd.contains(ip)) {
                        queue.add(ip);
                        enqd.add(ip);
                    }
                    log.log(logLevel, String.format("  d[%d]=%d,  d[%d]=%d\n", ip, d[ip], i, d[i]));
                 //   assert(d[ip] < d[i]);
                    e2 = (float) w[ip] / (float) w[i];
                    e2 *= e;
                    log.log(logLevel, String.format("    e=(%.3f)*(w[%d]/w[%d]=%.3f)\n", e, ip, i, (float) w[ip] / (float) w[i]));
                    uv = new PairInt(ip, i);
                    log.log(logLevel, String.format("    edge=(%d, %d) w=%.3e\n", ip, i, e2));
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

    /**
     * implementation of unweighted graph edge scoring from Girvan-Newman 
     * algorithm, accepting a directionless graph.   The given source node is
     * the root of the search.  Each graph node without a predecessor
     * is iterated over also to calculate node and edge weights.  
    
    * The directionless method here has a longer runtime complexity than the
    * DAG method because it has to check for leaf nodes as nodes without 
    * links further from source than it is, so this adds an iteration over edges.
     
     * Reference is 2004 Newman and Girvan,
     * "Finding and evaluating community structure in networks".
     * 
     * The runtime complexity is <em>1 + the number of graph roots (i.e. nodes without predecessors)
     *   times (add here)</em>.
     * 
     * For more information and other graph scoring and distance algorithms and 
     * cluster finding (a.k.a. community finding) see also
      <pre>
       Chapter 10 from "Mining of Massive Datasets"
       by Leskovec, Rajaraman, and Ullman
       http://infolab.stanford.edu/~ullman/mmds/ch10n.pdf
       
       and
       
       2005 paper "Complex networks: Structure and dynamics" by
                Boccalettia, Latorab, Morenod, Chavezf, and  Hwanga
                Physics Reports
               
       and
       
       2010 WWW2010 conference paper "Empirical Comparison of Algorithms for
            Network Community Detection" by Leskovec, Lang, and Mahoney
     
     </pre>
     * @param adjacencyList an adjacency list for an unweighted DAG.
     * @param s source node to begin the graph traversal
     * @return scores for the edges of the graph, calculated using the Girvan-Newman algorithm.
     */
    public Results girvanNewmanUndirected(SimpleLinkedListNode[] adjacencyList, final int s) {
        
        //avg O(|E|);  worst: O(|V| + |E|)
        int[] rootIndexes = findRoots(adjacencyList, s);
        
        if (rootIndexes.length == 0) {
            rootIndexes = new int[]{s};
        }
        
        final int nV = adjacencyList.length;
        // init
        final int[] d = new int[nV];
        final int[] dBest = new int[nV];
        Arrays.fill(dBest, Integer.MAX_VALUE);
        
        // holds parents of index.  e.g. for root node 0, p[0] is empty.
        final TIntList[] p = new TIntList[nV];
        final int[] color = new int[nV];
        
        // index is vertex index
        final int[] w = new int[nV];
        final int[] wG = new int[nV];
        
        // key = edge vertexes (u, v) where u is closer to the root, value = edge weight
        final TObjectFloatMap<PairInt> wEdges = new TObjectFloatHashMap<PairInt>();
        
        final TIntSet members = new TIntHashSet();
        
        // runtime is # of roots * O(|V| + |E|)
        
        for (int src : rootIndexes) {
        
            Arrays.fill(w, 0);
            Arrays.fill(color, 0);
            Arrays.fill(d, Integer.MAX_VALUE);
            for (int i = 0; i < nV; ++i) {
                p[i] = new TIntArrayList();
            }
            
            color[src] = 1;
            d[src] = 0;
            dBest[src] = 0;
            w[src] = 1;
            
            log.log(logLevel, "root=" + src);
      
            // calc vertex weights
            final TIntList leaf = new TIntArrayList();
            final Queue<Integer> queue = new ArrayDeque<Integer>();
            queue.add(src);
            int u;
            
            while (!queue.isEmpty()) {
                //log.log(logLevel, "w=" + Arrays.toString(w));
                //log.log(logLevel, "d=" + Arrays.toString(d));
                //log.log(logLevel, "dBest=" + Arrays.toString(dBest));
                //log.log(logLevel, "color=" + Arrays.toString(color));
                u = queue.remove().intValue();
                members.add(u);
                log.log(logLevel, String.format("u=%d\n", u));
                SimpleLinkedListNode vNode = adjacencyList[u];
                // for adjacencyLists that are directionless, edges from u to v and v to u exist,
                // so need to make one additional check for a pseudo-leaf:
                //     a pseudo-leaf has no descendants that are not also predecessors.
                boolean isALeafNode = isALeafNode(u, adjacencyList, d);
                if (isALeafNode) {
                    leaf.add(u);
                    color[u] = 2;
                    log.log(logLevel, "    LEAF\n");
                    continue;
                }
                while (vNode != null && vNode.getKey() != -1) {
                    int v = vNode.getKey();
                    if (color[v] == 0) {
                        color[v] = 1;
                        d[v] = d[u] + 1;
                        w[v] = w[u];
                        queue.add(v);
                        log.log(logLevel, String.format("  v=%d\n", v));
                        if (d[v] < dBest[v]) {
                            dBest[v] = d[v];
                        }
                    } else if (d[v] == (d[u] + 1)) {
                        w[v] += w[u];
                    }
                    if (d[u] < d[v]) {
                        log.log(logLevel, String.format("adding parent=%d of node=%d (d[%d]=%d, d[%d]=%d])\n",
                            u, v, u, d[u], v, d[v]));
                        p[v].add(u);
                    }
                    vNode = vNode.getNext();
                }
                color[u] = 2;
            }
            assert(queue.isEmpty());
         
            //add weights from previous root traversals
            if (rootIndexes.length > 1) {
            //    log.log(logLevel, "\n  w=" + Arrays.toString(w));
                for (int i = 0; i < w.length; ++i) {
                    w[i] += wG[i];
                }
            //    log.log(logLevel, "  wG=" + Arrays.toString(wG));
            //    log.log(logLevel, "->w=" + Arrays.toString(w));
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
                log.log(logLevel, String.format("t leaf=%d\n", t));
                pIter = p[t].iterator();
                while (pIter.hasNext()) {
                    i = pIter.next();
                    log.log(logLevel, String.format("  p=%d\n", i));
                    if (!enqd.contains(i)) {
                        queue.add(i);
                        enqd.add(i);
                    }
                    e = (float) w[i] / (float) w[t];
                    assert(d[i] < d[t]);
                    uv = new PairInt(i, t);
                    wEdges.put(uv, e);
                    log.log(logLevel, String.format("  edge=(%d, %d) w=%.3e\n", i, t, e));
                }
            }
          //  log.log(logLevel, "w=" + Arrays.toString(w));
          //  log.log(logLevel, "d=" + Arrays.toString(d));
          //  log.log(logLevel, "q="+ queue.toString());
            if (logLevel.equals(Level.INFO)) {
                for (int iv=0; iv<nV;++iv) {
                    TIntList ps = p[iv];
                    if (!ps.isEmpty()) {
                        log.log(logLevel, String.format("   *p[%d]=%s", iv, ps.toString()));
                    }
                }
            }
            float e2;
            PairInt ij;
            while (!queue.isEmpty()) {
                i = queue.remove().intValue();
                e = 1;
                log.log(logLevel, String.format("  e0[%d]=%.3f", i, e));
                SimpleLinkedListNode jNode = adjacencyList[i];
                while (jNode != null && jNode.getKey() != -1) {
                    int j = jNode.getKey();
                    if (d[i] < d[j]) {
                        ij = new PairInt(i, j);
                        assert (d[i] < d[j]);
                        assert (wEdges.containsKey(ij));
                        e += wEdges.get(ij);
                        log.log(logLevel, String.format("  + (w[%d]/w[%d]=%.3f)", i, j, wEdges.get(ij)));
                    }
                    jNode = jNode.getNext();
                }
                log.log(logLevel, String.format("  \n   =>%.3f\n", e));

                pIter = p[i].iterator();
                while (pIter.hasNext()) {
                    ip = pIter.next();
                    if (!enqd.contains(ip)) {
                        queue.add(ip);
                        enqd.add(ip);
                    }
                    log.log(logLevel, String.format("  parent=d[%d]=%d,  child=d[%d]=%d\n", 
                        ip, d[ip], i, d[i]));
                    assert(d[ip] < d[i]);
                    e2 = (float) w[ip] / (float) w[i];
                    e2 *= e;
                    log.log(logLevel, String.format("    e=(%.3f)*(w[%d]/w[%d]=%.3f)\n", e, ip, i, (float) w[ip] / (float) w[i]));
                    uv = new PairInt(ip, i);
                    log.log(logLevel, String.format("    edge=(%d, %d) w=%.3e\n", ip, i, e2));
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

    private boolean isALeafNode(int u, SimpleLinkedListNode[] adjacencyList,
        int[] d) {
        
        // for directionless graphs, presumably the input adjacency list has been
        // created such that all U->V edges have V->U edges. 
        // Looking for the equivalent of a leaf w.r.t. a source node:
        //    would expect no u-> pairs where d[u] < d[v].
        
        //Set<PairInt> uvs = new HashSet<PairInt>();
        SimpleLinkedListNode vNode = adjacencyList[u];
        while (vNode != null && vNode.getKey() != -1) {
            int v = vNode.getKey();
            if (d[v] > d[u]) {
                return false;
            }
            vNode = vNode.getNext();
        }
        return true;
    }

    /*private void print(SimpleLinkedListNode[] adjacencyList) {
        
        StringBuffer sb = new StringBuffer();
        
        for (int u = 0; u < adjacencyList.length; ++u) {
            sb.append(String.format("u=%d v=[", u));
            SimpleLinkedListNode vNode = adjacencyList[u];

            while (vNode != null && vNode.getKey() != -1) {
                int v = vNode.getKey();
                sb.append(String.format("%d, ", v));
                vNode = vNode.getNext();
            }
            sb.append(String.format("]\n"));
        }
        System.out.println(sb.toString());
    }*/
    
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
