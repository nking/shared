package algorithms.shortestPaths;

import algorithms.bipartite.MinHeapForRT2012;
import algorithms.heapsAndPQs.HeapNode;
import algorithms.misc.MiscMath0;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * 
  An implementation of the A* algorithm using the "uniform cost search" pattern
  and heuristics given to the code.
  
  http://en.wikipedia.org/wiki/A*_search_algorithm
 *
 * "As A* traverses the graph, it follows a path of the lowest known cost,
 * keeping a sorted priority queue [minimum heap] of alternate path segments
 * along the way.   If, at any point, a segment of the path being traversed has
 * a higher cost than another encountered path segment, it abandons the
 * higher-cost path segment and traverses the lower-cost path segment instead.
 * This process continues until the goal is reached."
 *   -- maintains a minimum heap of nodes to be traversed = the open set.
 *   -- the implementation below uses a breadth first search w/ depth=1
 *
 * If heuristics aren't given to the code, it makes an assumption that all
 * nodes lie within a straight line distance to the destination node and hence
 * calculates the heuristic based upon that distance.
 *
 * * Variables:
 *     g[n] is the shortest distance traveled path from the sourceIndx to the
 *          node n
 *     h[n] is the smallest estimated cost from the node n to destinationIndx
 *     f[n] is the total search cost from sourceIndex to node n
 *          f(n) = g(n) + h(n)
 * Goal:
 *     find the path that creates the smallest f[destinationIndx]

*  <pre>
     for more information on heuristics, see:
        http://theory.stanford.edu/~amitp/GameProgramming/Heuristics.html
        https://en.wikipedia.org/wiki/Admissible_heuristic
        https://en.wikipedia.org/wiki/Consistent_heuristic
        http://ai.stanford.edu/~latombe/cs121/2011/slides/D-heuristic-search.pdf
        https://michael.kim/blog/puzzle
        
     "admissible" means roughly that the heuristic underestimates the cost to 
      the goal or formally that h(v_{i},v_{g}) .leq. optimal distv_{i},v_{g})
      for all (v_{i},v_{g})}(v_{i},v_{g}) where {i,g} are in [0,1,...,n].   
     
     from http://theory.stanford.edu/~amitp/GameProgramming/Heuristics.html
     
     h(n) is an estimate of the minimum cost from any vertex n to the goal.
      
     - if h(n) is 0, then f(n) = g(n) + 0, so A* turns into Dijkstra’s (UCS in this case) Algorithm, 
         which is guaranteed to find a shortest path.
     - If h(n) is always lower than (or equal to) the cost of moving from n to the goal, 
         then A* is guaranteed to find a shortest path. 
         The lower h(n) is, the more node A* expands, making it slower.
     - If h(n) is exactly equal to the cost of moving from n to the goal, 
         then A* will only follow the best path and never expand anything else, 
         making it very fast. Although you can’t make this happen in all cases, 
         you can make it exact in some special cases. 
     - If h(n) is sometimes greater than the cost of moving from n to the goal, 
         then A* is not guaranteed to find a shortest path, but it can run faster.
     - If h(n) is very high relative to g(n), then only h(n) plays a role, 
         and A* turns into Greedy Best-First-Search.
     </pre>

 Not that if your intention is motion planning and the system is imperfect in its motion
 results (e.g. a move to a cell has error and another cell is reached instead)
 then one should consider a Markov Decision Process (MDP) model if the state of the
 world is known perfectly, else a Partially Observable Markov Decision Process (POMDP) model
 if the state of the world is not perfectly known.

 * @author nichole
 */
public class AStar {
    
    /**
     * directed acyclic graph where the main index of the object array
     * is the vertex number, and each vertex has a linked list of outgoing 
     * edges connecting to the vertexes held in each linked list node key.
     */
    protected SimpleLinkedListNode[] graph = null;
    
    /* edge weights in format of outer array index being the
     * first vertex and the map within having a key being the 2nd vertex of the
     * edge with value being the edge weight.
    */

    /**
     *
     */

    protected TIntIntMap[] w = null;
    
    // f(n) = g(n) + h(n)

    /**
     *
     */
    protected int[] f = null;
    
    // g(n) = distance from source to node n

    /**
     *
     */
    protected int[] g = null;
    
    // h(n) = heuristic distance from node n to destination

    /**
     *
     */
    protected int[] h = null;
    
    /**
     *
     */
    protected int[] visited = null;

    /**
     *
     */
    protected int[] predecessor = null;   
    
    /**
     *
     */
    protected int src = -1;

    /**
     *
     */
    protected int dest = -1;
    
    private int sentinel = Integer.MAX_VALUE;
        
    // this is recalculated in constructor
    private int maxValue = sentinel - 1;

    // key is cost of path so far plus the edge weight

    /**
     *
     */
    protected MinHeapForRT2012 heap = null;

    // refs to nodes internal to heap for decrease key operations

    /**
     *
     */
    protected HeapNode[] nodes = null;
    
    private Logger log = Logger.getLogger(getClass().getSimpleName());
    
    private Level logLevel = Level.FINEST;
    
    /**
     *
     @param dAG directed acyclic graph where the main index of the object array
     * is the vertex number, and each vertex has a linked list of outgoing 
     * edges connecting to the vertexes held in each linked list node key.
     * Note that all vertexes, including edge vertexes, must be present as an
     * index of the array dAG, i.e. all vertexes must have numerical value 
     * less than dAG.length.
     @param weights the edge weights in format of outer array index being the
     * first vertex and the map within having a key being the 2nd vertex of the
     * edge with value being the edge weight.
     @param heuristics an array of feasible heuristics with the array index 
     * being the node index.
     <pre>
     for more information, see:
        http://theory.stanford.edu/~amitp/GameProgramming/Heuristics.html
        http://ai.stanford.edu/~latombe/cs121/2011/slides/D-heuristic-search.pdf
        https://michael.kim/blog/puzzle
     </pre>
     @param sourceVertex the source vertex index
     @param destVertex the destination vertex index
     */
    public AStar(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, 
        int[] heuristics,
        int sourceVertex, int destVertex) {
        
        if (dAG == null || dAG.length == 0) {
            throw new IllegalArgumentException("dAG cannot be null");
        }
        if (sourceVertex < 0 || sourceVertex >= dAG.length) {
            throw new IllegalArgumentException("sourceIndex cannot be null");
        }
        
        /*
        initialize single source (g, s)
          initial node is sourceVertex and has cost 0
          initialize the priority queue and add the source node
          initialize visited array to -1
          note: the cost function is g(n), the sum of the weights of the 
            edges from the source node to node n along the shortest currently 
            known path. g(v) = g(u) + w(u, v).
        loop do
          if priority queue empty, return failure
          node = queue.extractMin
          if (visited[goal] is completed) done
          set visited[node] to explored
          for each v adjacent to node
            g(v) = g(u) + w(u,v)
            if (v is not explored and is not in queue)
              insert v into queue
            else if child is in queue and dist[child] > path_cost
              decrease key of v to path_cost
        */ 

        // only source node is added to heap initially:
        init(dAG, weights, heuristics, sourceVertex, destVertex);        
    }
        
    /**
     * find the single shortest path in dAG with edge weights w starting from s
     * to destination vertex.  if destination supplied in constructor was
     * outside of graph index limits, all reachable nodes will be searched.
     */
    public void find() {
                               
        while (heap.getNumberOfNodes() > 0) {

            HeapNode uNode = heap.extractMin();
            
            int u = ((Integer)uNode.getData()).intValue();

            log.log(logLevel, "u: " + toString(u));
            
            if (u == dest) {
                log.log(logLevel, "exit heap.n=" + heap.getNumberOfNodes());
                return;
            }

            visited[u] = 2;
            
            // null the entry in nodes so it isn't used in decrease key
            nodes[u] = null;
            
            TIntIntMap uWeights = w[u];
            
            if (uWeights == null) {
                continue;
            }

            SimpleLinkedListNode vNode = graph[u];
            
            while (vNode != null && vNode.getNumberOfKeys() > 0) {
            
                int v = vNode.getKey();
               
                if (!uWeights.containsKey(v)) {
                    throw new IllegalStateException("no weight found for edge " 
                    + u + " to " + v);
                }

                //f(n) = g(n) + h(n)
                int gUPlusWUV = g[u];
                if (gUPlusWUV < sentinel) {
                    gUPlusWUV += uWeights.get(v);
                }
                
                if (gUPlusWUV >= g[v]) {
                    vNode = vNode.getNext();
                    continue;
                }
                                
                if (visited[v] == 0) {
                    
                    visited[v] = 1;
                    
                    g[v] = gUPlusWUV;
                    f[v] = g[v] + h[v];
                    
                    int key = f[v];
                    HeapNode node = new HeapNode(key);
                    node.setData(Integer.valueOf(v));
                    nodes[v] = node;
                    predecessor[v] = u;
                    heap.insert(node);
                } else if (nodes[v] != null) {
                    //assert(visited[v] == 1);
                    g[v] = gUPlusWUV;
                    f[v] = g[v] + h[v];
                    predecessor[v] = u;
                    
                    // if h(n) is exact, can use remove here instead of decreaseKey.
                    heap.decreaseKey(nodes[v], f[v]);
                }
                
                vNode = vNode.getNext();
            }
        }
    }
    
    private void init(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, 
        int[] heuristics, int sourceVertex,
        int destVertex) {
        
        graph = dAG.clone();
        for (int i = 0; i < dAG.length; ++i) {
            graph[i] = new SimpleLinkedListNode(dAG[i]);
        }
        w = weights.clone();
        for (int i = 0; i < weights.length; ++i) {
            if (weights[i] != null) {
                w[i] = new TIntIntHashMap(weights[i]);
            }
        }
        src = sourceVertex;
        dest = destVertex;
        
        // h(n)
        h = Arrays.copyOf(heuristics, heuristics.length);
        
        maxValue = calcUpperLimitKeyValue();
        
        sentinel = maxValue + 1;
    
        predecessor = new int[graph.length];
        visited = new int[graph.length];
        Arrays.fill(predecessor, -1);
        Arrays.fill(visited, 0);
        
        // f(n)
        f = new int[graph.length];
        Arrays.fill(f, sentinel);
                
        // g(n)
        g = new int[graph.length];
        Arrays.fill(g, sentinel);
                
        // presumably this is always correct:
        h[src] = sentinel;
        
        g[src] = 0;
        f[src] = h[src];
                
        initHeap();
    }
    
    private void initHeap() {
        
        int n = graph.length;
                
        int nBits = (int)Math.ceil(Math.log(maxValue/Math.log(2)));
        
        //int maxValue, int approxN, int maxNumberOfBits
        heap = new MinHeapForRT2012(sentinel, n, nBits);

        nodes = new HeapNode[n];

        int key = f[src];

        HeapNode node = new HeapNode(key);
        node.setData(Integer.valueOf(src));

        heap.insert(node);

        nodes[src] = node;
    }
    
    /**
     * get shortest path from source to idx
     @param idx
     @return 
     */
    public int getDistanceFromSrc(int idx) {
        if (idx < 0 || idx >= graph.length) {
            throw new IllegalArgumentException("idx must be within bounds of graph vertices");
        }
        if (predecessor == null) {
            throw new IllegalStateException("find must be run before this method can be used");
        }
        return g[idx];
    }
    
    /**
     * get shortest path from source to destIndex
     @param destVertex
     @return 
     */
    public int[] getShortestPathToVertex(int destVertex) {
        if (destVertex < 0 || destVertex >= graph.length) {
            throw new IllegalArgumentException("destIndex cannot be null");
        }
        if (predecessor == null) {
            throw new IllegalStateException("find must be run before this method can be used");
        }
                
        log.log(logLevel, "    f[]=" + Arrays.toString(f));
        log.log(logLevel, "    g[]=" + Arrays.toString(g));
        
        int[] p = new int[graph.length];
        p[p.length - 1] = destVertex;
                
        for (int i = p.length - 2; i > -1; --i) {
            if (destVertex == src) {
                int len = p.length - 1 - i;
                int[] t = new int[len];
                System.arraycopy(p, i + 1, t, 0, len);
                return t;
            } else if (destVertex == -1) {
                throw new IllegalStateException("path did not complete correctly");
            }
            p[i] = predecessor[destVertex];
            destVertex = p[i];
        }
        
        if (p[0] != src) {
            throw new IllegalStateException("path did not complete correctly for destIndex");
        }
        
        return p;
    }
    
    /**
     *
     @param vertexes
     @return
     */
    public int getSumOfPath(int[] vertexes) {
        if (vertexes == null) {
            throw new IllegalArgumentException("vertexes cannot be null");
        }
        if (f == null) {
            throw new IllegalStateException("find must be run before this method can be used");
        }
        int sum = 0;
        int u, v;
        for (int i = 1; i < vertexes.length; ++i) {
            u = vertexes[i - 1];
            v = vertexes[i];
            
            TIntIntMap uWeights = w[u];
            
            if (uWeights == null || !uWeights.containsKey(v)) {
                throw new IllegalStateException("no weight found for edge " 
                    + u + " to " + v);
            }
            
            sum += uWeights.get(v);
        }
        return sum;
    }

    private int calcUpperLimitKeyValue() {
        
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < w.length; ++i) {
            TIntIntMap map = w[i];
            if (map == null) {
                continue;
            }
            TIntIntIterator iter = map.iterator();
            int value;
            for (int j = 0; j < map.size(); ++j) {
                iter.advance();
                value = iter.value() + h[i];
                if (value > max) {
                    max = value;
                }
            }
        }
        
        // if # of bits in max * g.length is larger than max_integer,
        //   use max_integer, else use it
        int b1 = MiscMath0.numberOfBits(max);
        int b2 = MiscMath0.numberOfBits(graph.length);
        if ((b1 + b2) > 31) {
            max = Integer.MAX_VALUE - 1;
        } else {
            max *= graph.length;
        }
        
        return max;
    }
    
    private String toString(int u) {
        StringBuffer sb = new StringBuffer();
        sb.append("node=").append(u).append(": visited=").append(visited[u]).
               append(", dist=").append(f[u])
                .append(", prev=").append(predecessor[u])
                .append(" is in queue=").append(!(nodes[u] == null));
        return sb.toString();
    }
}
