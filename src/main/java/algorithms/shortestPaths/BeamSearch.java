package algorithms.shortestPaths;

import algorithms.util.SimpleLinkedListNode;

import algorithms.FixedSizeSortedVector; 
import algorithms.util.LinkedListCostNode;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.Arrays;

import gnu.trove.map.TIntIntMap;
import thirdparty.edu.princeton.cs.algs4.Queue;

        
/**
from https://www.wikiwand.com/en/Beam_search
Beam Search uses a greedy, breadth-first search to build its search tree,
but only keeps top k (beam size) nodes at each level in memory.
The next level will then be expanded from these k nodes.

Complete search variants of beam search are made by combining it with 
depth-first search resulting in beam stack search and depth-first beam search, 
or combining it with limited discrepancy search resulting in beam search 
using limited discrepancy backtracking (BULB). 
The resulting search algorithms are anytime algorithms that find good but 
likely sub-optimal solutions quickly, like beam search, then backtrack and 
continue to find improved solutions until convergence to an optimal solution.

breadth-first search:
 given a graph, it visits every node reachable from s and computes the
 distance to each of those nodes.  the predecessor tree is those nodes
 reachable from s.   the distance array can be used to find the
 shortest path.

 * @author nichole
 */
public class BeamSearch {

    /**
     *
     */
    protected int[] dist = null;

    /**
     *
     */
    protected int[] predecessor = null;

    /**
     *
     */
    protected int[] visited = null;

    /**
     *
     */
    protected SimpleLinkedListNode[] adjacencyList = null;

    /**
     * source index
     */
    protected final int src;

    /**
     *
     */
    protected final int beamSize;

    /**
     * constructor with adjacency list, with default equal cost edges.
     @param theAdjacencyList
     @param sourceNodeIndex start index for search
     @param beamSize in the breadth-first search, for each set of child nodes, 
     * only the smallest cost nodes are expanded to continue search and the
     * limited number of those smallest cost nodes is beamSize.
     */
    public BeamSearch(SimpleLinkedListNode[] theAdjacencyList, int sourceNodeIndex,
            final int beamSize) {

        this.adjacencyList = theAdjacencyList;
        this.dist = new int[adjacencyList.length];
        this.predecessor = new int[adjacencyList.length];
        this.visited = new int[adjacencyList.length];
        if (sourceNodeIndex < 0) {
            throw new IllegalArgumentException("sourceNodeIndex cannot be a negative number");
        }
        if (sourceNodeIndex >= adjacencyList.length) {
            throw new IllegalArgumentException("sourceNodeIndex must be an index within size limits of theAdjacencyList");
        }
        this.src = sourceNodeIndex;
        
        if (beamSize <= 0) {
            throw new IllegalArgumentException("beamSize must be a positive number greater than 0");
        }
        this.beamSize = beamSize;
    }
    
    /**
     * constructor with edge costs in adjacency list.
     @param theAdjacencyList adjacency list that includes edge costs
     @param sourceNodeIndex start index for search
     @param beamSize in the breadth-first search, for each set of child nodes, 
     * only the smallest cost nodes are expanded to continue search and the
     * limited number of those smallest cost nodes is beamSize.
     */
    public BeamSearch(LinkedListCostNode[] theAdjacencyList, int sourceNodeIndex,
            final int beamSize) {

        this.adjacencyList = theAdjacencyList;
        this.dist = new int[adjacencyList.length];
        this.predecessor = new int[adjacencyList.length];
        this.visited = new int[adjacencyList.length];
        if (sourceNodeIndex < 0) {
            throw new IllegalArgumentException("sourceNodeIndex cannot be a negative number");
        }
        if (sourceNodeIndex >= adjacencyList.length) {
            throw new IllegalArgumentException("sourceNodeIndex must be an index within size limits of theAdjacencyList");
        }
        this.src = sourceNodeIndex;
        
        if (beamSize <= 0) {
            throw new IllegalArgumentException("beamSize must be a positive number greater than 0");
        }
        this.beamSize = beamSize;
    }

    /**
     * add source index to the bfs tree.
     @return list of indexes searched
     */
    public TIntList search() {

        TIntList searched = new TIntArrayList();
                
        // initialize
        Arrays.fill(dist, Integer.MAX_VALUE);
        Arrays.fill(predecessor, -1);
        Arrays.fill(visited, 0);
        visited[src] = 1;
        dist[src] = 0;
        Queue queue = new Queue();
        queue.enqueue(src);

        int dUPlusWUV;

        while (!queue.isEmpty()) {

            final int u = (Integer)queue.dequeue();

            if (visited[u] == 2) {
                continue;
            }

            searched.add(u);

            SimpleLinkedListNode vNode = adjacencyList[u];

            FixedSizeSortedVector<Index> sorted = new FixedSizeSortedVector(beamSize, Index.class);
        
            while (vNode != null && vNode.getNumberOfKeys() > 0) {

                final int v = vNode.getKey();

                if (vNode instanceof LinkedListCostNode) {
                    dUPlusWUV = dist[u] + ((LinkedListCostNode)vNode).getCost();
                } else {
                    dUPlusWUV = dist[u] + 1;
                }
                if (dist[v] > dUPlusWUV) {
                    dist[v] = dUPlusWUV;
                    predecessor[v] = u;
                    if (visited[v] == 0) {
                        visited[v] = 1;
                        sorted.add(new Index(v, dist[v]));
                    }
                }

                vNode = vNode.getNext();
            }
            
            Index[] sortedArray = sorted.getArray();
            for (int z = 0; z < sorted.getNumberOfItems(); z++) {
                queue.enqueue(sortedArray[z].idx);
            }

            visited[u] = 2;
        }
        
        return searched;
    }

    private class Index implements Comparable<Index> {

        private final int idx;
        private final int dist;

        public Index(int index, int pathDist) {
            this.idx = index;
            this.dist = pathDist;
        }
        
        @Override
        public int compareTo(Index c) {
            // compare distances. this is within scope of d
            return Integer.compare(dist, c.dist);
        }
        
    }


    /**
     * get shortest path from source to destIndex
     @param destVertex
     @return
     */
    public int[] getShortestPathToVertex(int destVertex) {
        int n = adjacencyList.length;
        if (destVertex < 0 || destVertex >= n) {
            throw new IllegalArgumentException("destIndex cannot be null");
        }
        if (predecessor == null) {
            throw new IllegalStateException("find must be run before this method can be used");
        }

        int[] p = new int[n];
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
        if (dist == null) {
            throw new IllegalStateException("find must be run before this method can be used");
        }
        int sum = 0;
        int u, v;
        SimpleLinkedListNode vNode;
        int dUV;
        for (int i = 1; i < vertexes.length; ++i) {
            u = vertexes[i - 1];
            v = vertexes[i];

            vNode = adjacencyList[u];
            if (vNode == null) {
                throw new IllegalStateException("no edge from " + u + " to " + v);
            }

            if (vNode instanceof LinkedListCostNode) {
                dUV = ((LinkedListCostNode)vNode).getCost();
            } else {
                dUV = 1;
            }

            sum += dUV;
        }
        return sum;
    }

}
