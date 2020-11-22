package algorithms.graphs;

import algorithms.util.SimpleLinkedListNode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * A strongly connected (a.k.a. diconnected) graph is one in which every vertex 
 * is reachable from every other vertex.
 * The strongly connected components of an arbitrary directed graph form a 
 * partition into subgraphs that are themselves strongly connected.
 * 
 * It is possible to test the strong connectivity of a graph, or to find 
 * its strongly connected components, in linear time (that is, Θ(V+E)).
 * 
 * A directed graph is called strongly connected if there is a path in each 
 * direction between each pair of vertices of the graph. That is, a path exists 
 * from the first vertex in the pair to the second, and another path exists 
 * from the second vertex to the first. 
 * 
 * A pair of vertices u and v in a directed graph G that may not itself be 
 * strongly connected, are said to be strongly connected to each other if there 
 * is a path in each direction between them.
 * 
 * Tarjan's strongly connected components algorithm traverses the graph once, 
 * using a stack in which the vertices are not removed unless found to be
 * strongly connected.
 * 
 * Runtime complexity is O(|E| + |V|) (average?), and the worse case runtime
 * complexity is one in which the graph is completely connected, strongly,
 * O(|V|*(2 + 5w)) where w is the word size in bits.
 * 
 * https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
 * 
 * Further information on use in processing very large datasets such as 
 * in transitive closure is in the book "Mining of Massive Datasets" 
 * by Jure Leskovec, Anand Rajaraman, Jeff Ullman
 * http://infolab.stanford.edu/~ullman/mmds/ch10n.pdf
 * A strongly connected component (SCC) of a graph is a set of nodes S such that:
      1. Every node in S reaches every other node in S.
      2. S is maximal, in the sense that there is no node outside S that both
         reaches every node in S and is reached by every node in S.
* 
 * 
 * first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright

 * @author nichole
 */
public class StronglyConnectedComponents {
    
    private SimpleLinkedListNode[] g;
    
    /**
     * index numbers the nodes consecutively in the order in which they are discovered
     */
    private int[] index;
    
    /**
    lowLink represents the smallest index of any node known to be reachable 
       from v through v's DFS subtree, including v itself.
    Therefore v must be left on the stack if v.lowlink less than v.index, 
    whereas v must be removed as the root of a strongly connected component 
    if v.lowlink equals v.index. 
     
    The value v.lowlink is computed during the depth-first search from v,
    as this finds the nodes that are reachable from v.
     */
    private int[] lowLink;
    private int[] onStack;
    
    private int time;
    private List<SimpleLinkedListNode> scc;
    private int[] inSCC;
    private Stack<Integer> stack;
    
    private Logger log = Logger.getLogger(getClass().getSimpleName());
    
    private Level logLevel = Level.FINEST;
    
    /**
     *  find the strongly connected components.
     * @param graph graph in the form of an adjacency list
     * @return strongly connected components as a list of linked-lists.  Note
     * that the return order of the vertices in the list and their 
     * linked-lists is the reverse of the topological sort.
     * 
     * following https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
     * 
     */
    public List<SimpleLinkedListNode> find(SimpleLinkedListNode[] graph) {
        
        check(graph);
        
        g = Arrays.copyOf(graph, graph.length);
        lowLink = new int[g.length];
        index = new int[g.length];
        Arrays.fill(index, -1);
        Arrays.fill(lowLink, -1);
        onStack = new int[g.length];
        scc = new ArrayList<SimpleLinkedListNode>();
        inSCC = new int[g.length];
        time = 0;
        stack = new Stack<Integer>();
        
        for (int i = 0; i < g.length; ++i) {
            strongConnect(i);
        }
        
        return scc;
    }
    
    private void strongConnect(int v) {
        
        index[v] = time;
        lowLink[v] = time;
        time++;
        
        stack.push(v);
        onStack[v] = 1;
        
        log.log(logLevel, "v:" + toString(v));
                
        SimpleLinkedListNode wNode = g[v];
        while (wNode != null && wNode.getKey() != -1) {
            int w = wNode.getKey();
            log.log(logLevel, "    w=" + toString(w));            
            if (index[w] == -1) {
                // Successor w has not yet been visited; recurse on it
                strongConnect(w);
                lowLink[v] = (lowLink[v] <= lowLink[w]) ? lowLink[v] : lowLink[w];
            } else if (onStack[w] == 1) {
                // w is in stack S and hence in the current SCC
                // If w is not on stack, then (v, w) is a cross-edge in the 
                // DFS tree and must be ignored.
                lowLink[v] = (lowLink[v] <= index[w]) ? lowLink[v] : index[w];
            }
            wNode = wNode.getNext();
        }
        
        // If v is a root node, pop the stack and generate an SCC
        if (lowLink[v] == index[v] && inSCC[v] == 0) {
            log.log(logLevel, "    START scc " + v);
            SimpleLinkedListNode sccNode = new SimpleLinkedListNode();
            scc.add(sccNode);
            inSCC[v] = 1;
            int w;
            do {
                if (stack.isEmpty()) {
                    break;
                }
                w = stack.pop();
                onStack[w] = 0;
                sccNode.insert(w);
                inSCC[w] = 1;
                log.log(logLevel, "    add " + w + " to scc");
            } while (w != v);
            log.log(logLevel, "stack size=" + stack.size());
        }
    }

    private void check(SimpleLinkedListNode[] graph) {

        if (graph == null) {
            throw new IllegalArgumentException("graph cannot be null");
        }
        
        for (SimpleLinkedListNode node : graph) {
            if (node == null) {
                throw new IllegalArgumentException("graph nodes cannot be null");
            }
        }
    }
    
    private String toString(int node) {
        
        StringBuilder sb = new StringBuilder();
        sb.append("node=").append(node).append(", lowLink=").append(lowLink[node])
            .append(", index=").append(index[node])
            .append(", onStack=").append(onStack[node])
            .append(" (time=").append(time).append(");");
        return sb.toString();
    }
}
