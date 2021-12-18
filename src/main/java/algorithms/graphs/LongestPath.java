package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TObjectDoubleIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;

/**
 * In graph theory and theoretical computer science, the longest path problem 
 is the problem of finding a simple path of maximum length in a given graph. 
 A path is called simple if it does not have any repeated vertices; the length 
 of a path may either be measured by its number of edges, or (in weighted graphs) 
 by the sum of the weights of its edges.

 The problem is NP-hard and the decision version of the problem, which asks 
 whether a path exists of at least some given length, is NP-complete.

there is  a linear time approx solution for directed acyclic graphs, which has 
important applications in finding the critical path in scheduling problems.

proof of NP-hardness:
      reduction from Hamiltonian path problem (i.e. a cycle that includes all vertices once).
      Decision version: if input graph G has a path of k or more edges.
      
 * @author nichole
 */
public class LongestPath {
    
    /**
     * given a DAG, find the longest simple path of maximum length where lengths
     * are given by the edge weights in graph g.
     * <pre>
     * https://en.m.wikipedia.org/wiki/Longest_path_problem
     * </pre>
     * runtime complexity is O(|V| + |E|).
     *
     * @param g input weighted DAG where key = a directed edge of vertexes (u,v)
     * and key = weight of the key edge where the weights are non-negative).
     * @return vertex indexes of longest path through the DAG.
     */
    public static int[] solve(TObjectDoubleMap<PairInt> g) {
        // convert g to an adjacency list where every vertex has an entry even if empty.
        // key = vertex
        SimpleLinkedListNode[] nodes = convert(g);
       
        //O(V + E)
        TopologicalSort ts = new TopologicalSort(nodes);
                
        int[] tsIdxs = ts.sort();
        
        return solve(g, nodes, tsIdxs, tsIdxs[0]);
    }
    
    /**
     * given a DAG, find the longest simple path of maximum length from the
     * given start vertex to any other vertex.
     * <pre>
     * https://en.m.wikipedia.org/wiki/Longest_path_problem
     * </pre>
     * runtime complexity is O(|V| + |E|)
     *
     * @param g input weighted DAG where key = a directed edge of vertexes (u,v)
     * and key = weight of the key edge where the weights are non-negative).
     * @param srcIdx the graph vertex index for the start (a.k.a. source) node.
     * @return vertex indexes of longest path through the DAG.
     */
    public static int[] solve(TObjectDoubleMap<PairInt> g, int srcIdx) {
        
        // convert g to an adjacency list where every vertex has an entry even if empty.
        // key = vertex
        SimpleLinkedListNode[] nodes = convert(g);
        if (nodes.length <= srcIdx) {
            throw new IllegalArgumentException("index srcIdx was not found in g");
        }
        
        //O(V + E)
        TopologicalSort ts = new TopologicalSort(nodes);
                
        int[] tsIdxs = ts.sort();
        
        return solve(g, nodes, tsIdxs, srcIdx);
    }
    
    /**
     * given a DAG, find the longest simple path of maximum length from the
     * given start vertex to any other vertex.
     * <pre>
     * https://en.m.wikipedia.org/wiki/Longest_path_problem
     * </pre>
     *
     * @param g input weighted DAG where key = a directed edge of vertexes (u,v)
     * and key = weight of the key edge where the weights are non-negative).
     * @param nodes an adjacency list extracted from graph g.
     * @param tsIdxs the topologically sorted vertexes of graph g.
     * @param srcIdx the graph vertex index for the start (a.k.a. source) node.
     * @return vertex indexes of longest path through the DAG.
     */
    static int[] solve(TObjectDoubleMap<PairInt> g, 
        SimpleLinkedListNode[] nodes, int[] tsIdxs, int srcIdx) {
        
        /*
        Find a topological ordering of the given DAG.
          For each vertex v of the DAG, in the topological ordering,
             compute the length of the longest path ending at v by looking at its
             incoming neighbors and adding one to the maximum length recorded for
             those neighbors.
             If v has no incoming neighbors, set the length of the longest path
             ending at v to zero. In either case, record this number so that
             later steps of the algorithm can access it.
         Once this has been done, the longest path in the whole DAG may be obtained
             by starting at the vertex v with the largest recorded value,
             then repeatedly stepping backwards to its incoming neighbor with the
             largest recorded value, and reversing the sequence of vertices found in this way.
        */
        
        
        // find srcIdx in tsIdxs
        int sIdx = -1;
        int i;
        for (i = 0; i < tsIdxs.length; ++i) {
            if (tsIdxs[i] == srcIdx) {
                sIdx = i;
                break;
            }
        }
        //System.out.printf("srcIdx=%d, sIdx=%d, ts=%s\n", srcIdx, sIdx, Arrays.toString(tsIdxs));
        assert(sIdx > -1);
        
        double[] dist = new double[nodes.length];
        int[] prev = new int[nodes.length];
        Arrays.fill(prev, -1);
        
        // for pred, fill in the map as progress from v 
        TIntObjectMap<TIntSet> prevMap = new TIntObjectHashMap<TIntSet>();
        
        double maxDistS = Double.NEGATIVE_INFINITY, distV, wUV, d;
        int maxDistIdx = -1;
        SimpleLinkedListNode nhbr;
        int v, v2, u;
        TIntSet set;
        TIntIterator iter;
        for (i = sIdx; i < tsIdxs.length; ++i) {
            
            v = tsIdxs[i];
            
            // add u to prevMap for the next v in tsIdxs to use
            nhbr = nodes[v];
            while (nhbr != null && nhbr.getKey() != -1) {
                v2 = nhbr.getKey();
                set = prevMap.get(v2);
                if (set == null) {
                    set = new TIntHashSet();
                    prevMap.put(v2, set);
                }
                set.add(v);
                nhbr = nhbr.getNext();
            }
            
            distV = Double.NEGATIVE_INFINITY;
            set = prevMap.get(v);
            if (set != null) {
                iter = set.iterator();
                while (iter.hasNext()) {
                    u = iter.next();
                    wUV = g.get(new PairInt(u, v));
                    assert(g.containsKey(new PairInt(u, v)));
                    d = dist[u] + wUV;
                    if (d > distV) {
                        distV = d;
                        prev[v] = u;
                    }
                }
            }
            
            if (distV == Double.NEGATIVE_INFINITY) {
                distV = 0;
            }
            dist[v] = distV;
            
            if (distV > maxDistS) {
                maxDistS = distV;
                maxDistIdx = v;
            }
        }
        
        //System.out.printf("src=%d  maxDist=%.3f  maxDistIdx=%d\n", srcIdx, maxDistS, maxDistIdx);
        //System.out.printf("dist=%s\n", FormatArray.toString(dist, "%.3f"));
        //System.out.printf("pred=%s\n", Arrays.toString(prev));
        
        // back track from maxS to S then reverse the nodes and return that
        TIntList path = new TIntArrayList();
        int idx = maxDistIdx;//4->3->7
        while (idx > -1) {
            path.add(idx);
            idx = prev[idx];
        }
        path.reverse();
        return path.toArray();
    }

    static SimpleLinkedListNode[] convert(TObjectDoubleMap<PairInt> g) {
        
        PairInt p;
        
        // count vertexes
        TIntSet vs = new TIntHashSet();
        TObjectDoubleIterator<PairInt> iter = g.iterator();
        while (iter.hasNext()) {
            iter.advance();
            p = iter.key();
            vs.add(p.getX());
            vs.add(p.getY());
        }
        
        int n = vs.size();
        int i;
        // initialize output
        SimpleLinkedListNode[] nodes = new SimpleLinkedListNode[n];
        for (i = 0; i < n; ++i) {
            nodes[i] = new SimpleLinkedListNode();
        }
        
        SimpleLinkedListNode node;
        iter = g.iterator();
        while (iter.hasNext()) {
            iter.advance();
            p = iter.key();
            nodes[p.getX()].insert(p.getY());
        }
        
        return nodes;
    }
}
