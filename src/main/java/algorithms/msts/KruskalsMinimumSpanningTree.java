package algorithms.msts;

import algorithms.disjointSets.DisjointForest;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.disjointSets.UnionFind;
import algorithms.sort.MiscSorter;
import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TObjectDoubleIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TIntObjectHashMap;

import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.stream.IntStream;

/**
*
* minimum spanning tree is the subset of edges in a weighted undirected graph
 * that connect all vertexes for a total minimum cost (sum of edge weights).
 * 
  Kruskal's grows a forest by sorting the edges first and then adding edges that 
  are not yet connected to the tree.
  
  * Runtime Complexity is O(|E| lg_2|E|),
      which for sparse graphs having |E| .lt. |V|^2,
      gives O(|E| lg_2|V|).
    so for dense graphs, one should prefere Prim's.
      
  Best time Kruskal's is O(|E| lg_2|V|).
  Best time  Prim's w/ fib heaps is O(|E| + |V|lg_2|V|).
  
 * @author nichole
 */
public class KruskalsMinimumSpanningTree {
    
    /**
     * find a minimum spanning tree using Kruskal's algorithm.
     @param graph an adjacency list for the graph, where the index of the array
     * is the vertex number and the keys within each vertex's list node are the
     * vertex on the other end of an edge.
     @param edgeWeights map of edge weights with key being pairs of the
     * vertex numbers.
     @return a minimum spanning tree of the weighted given graph
     */
    public static TIntObjectMap<SimpleLinkedListNode> mst(
        SimpleLinkedListNode[] graph, TObjectDoubleMap<PairInt> edgeWeights) {
                
        /* pseudo code form Cormen, Leiserson, Rivest, and Stein "Introduction to Computer Science".
         *
         *  initialize set A to the empty set and create |V| trees, one containing each vertex
         *
         *      A <-0
         *      for each vertex v member of V[G]
         *          do Make-Set(v) 
         *      sort the edges of E into nondecreasing order by weight w
         *      for each edge (u, v) member of E, taken in nondecreasing order by weight
         *          do if Find-Set(u) != Find-Set(v)  
         *              then A <- A union {(u,v)}
         *                  Union(u,v)
         *      return A
        */
                
        DisjointForest<Integer> forest = new DisjointForest<>();
        TIntObjectMap<DisjointSet2Node<Integer>> vertexMap = new TIntObjectHashMap<>();
        
        int nV = graph.length;
        
        DisjointSet2Node<Integer> uVertex, vVertex;
        int u, v;
        for (u = 0; u < nV; ++u) {
            uVertex = new DisjointSet2Node<>(u);
            vertexMap.put(u, uVertex);
            forest.makeSet(uVertex);
        }
        
        PairInt[] w = sortWeightsNonDecreasing(edgeWeights);
        
        TIntObjectMap<SimpleLinkedListNode> a = new TIntObjectHashMap<>();
        
        DisjointSet2Node<Integer> aUV;
        long sum = 0;
        int nMSTEdges = 0;
        for (PairInt e : w) {
            u = e.getX();
            v = e.getY();
            uVertex = forest.findSet(vertexMap.get(u));
            vVertex = forest.findSet(vertexMap.get(v));
            
            if (!uVertex.equals(vVertex)) {
                
                aUV = forest.union(uVertex, vVertex);
                
                if (a.containsKey(u)) {
                    a.get(u).insert(v);
                } else {
                    SimpleLinkedListNode node = new SimpleLinkedListNode(v);
                    a.put(u, node);
                }
                sum += edgeWeights.get(new PairInt(u, v));
                nMSTEdges++;
            }
        }
        
        System.out.println("forest stats=" + forest.toString());
        System.out.printf("%d edges out of %d in minimumspanning tree.  sum=%d\n",
            nMSTEdges, edgeWeights.size(), sum);
        
        return a;
    }

    /**
     * given an edge weights matrix, return the edge indexes that make a minimum spanning tree.
     * @param edgeWeights edge weights given in format where each row contains
     *                    the 2 vertex numbers of the edge and the edge weight.
     *                    e.g. row 0 = [1,4, 100] states that edge 0 has
     *                    start vertex 1 and end vertex 4 and weight 100.
     * @return
     */
    public static List<Integer> mst(int[][] edgeWeights) {
        // [x,y, w]  sort by col2
        int n = edgeWeights.length;
        int[] sortedIdxs
                = IntStream.range(0, n).boxed()
                .sorted((i, j)-> Integer.compare(edgeWeights[i][2], edgeWeights[j][2]))
                .mapToInt(ele->ele).toArray();

        List<Integer> mstEdgeIndexes = new ArrayList<>();
        UnionFind uF = new UnionFind(n);
        for (int i = 0; i < n; ++i) {
            int idx = sortedIdxs[i];
            int u =edgeWeights[idx][0];
            int v = edgeWeights[idx][1];
            if (uF.find(u) != uF.find(v)) {
                mstEdgeIndexes.add(idx);
                uF.union(u, v);
            }
        }

        return mstEdgeIndexes;
    }

    static PairInt[] sortWeightsNonDecreasing(
        TObjectDoubleMap<PairInt> edgeWeights) {
        
        int n = edgeWeights.size();
        
        PairInt[] keys = new PairInt[n];
        double[] w = new double[n];
        
        TObjectDoubleIterator<PairInt> iter = edgeWeights.iterator();
        for (int i = 0; i < n; ++i) {
            iter.advance();
            keys[i] = iter.key();
            w[i] = iter.value();
        }
        
        int[] sortedIdxs = MiscSorter.mergeSortIncreasing(w);
        
        PairInt[] sortedKeys = new PairInt[n];
        for (int i = 0; i < n; ++i) {
            sortedKeys[i] = keys[sortedIdxs[i]];
        }
        
        return sortedKeys;
    }
    
    /**
     given an undirected weighted graph adjMap, find a minimum spanning tree, and return
     it as an adjacency map of the original graph.
     */
    public static Map<Integer, Map<Integer, Double>> 
        mst(Map<Integer, Map<Integer, Double>> adjMap, double[] outputSum) {

	/*
	 * sort the edges by weight and add them to the growing forest if not already
	 * present.*/
	int nEdges = 0;
	for (int u : adjMap.keySet()) {
	    nEdges += adjMap.get(u).size();
        }

	int[][] edges = new int[nEdges][];
	double[] weights = new double[nEdges];
	int i = 0;
	for (int u : adjMap.keySet()) {
            for (Map.Entry<Integer, Double> entry : adjMap.get(u).entrySet()) {
		edges[i] = new int[]{u, entry.getKey()};
		weights[i] = entry.getValue();
                ++i;
	    } 
	}

	List<int[]> tree = mst(edges, weights, outputSum);

	Map<Integer, Map<Integer, Double>> out = new HashMap<>();

        for (int[] edge : tree) {
	    if (adjMap.containsKey(edge[0]) && adjMap.get(edge[0]).containsKey(edge[1])) {
		out.putIfAbsent(edge[0], new HashMap<Integer, Double>());
		out.get(edge[0]).put(edge[1], adjMap.get(edge[0]).get(edge[1]));
            } else {
		out.putIfAbsent(edge[1], new HashMap<Integer, Double>());
		out.get(edge[1]).put(edge[0], adjMap.get(edge[1]).get(edge[0]));
            } 
	}	

	return out;
    }

    /**
     given an undirected weighted graph and weights, find a minimum spanning tree,
     and return it as edges of the original graph.
     */
    public static List<int[]> mst(int[][] edges, double[] weights, double[] outputSum) {

        int[] sortedIdxs = IntStream.range(0, edges.length) .boxed()
	    .sorted( (i, j) -> Double.compare(weights[i], weights[j]))
	    .mapToInt(ele -> ele).toArray();	

	// count number of vertices
	int nV = 0;
	Set<Integer> vS = new HashSet<>();
	for (int[] edge : edges) {
	    vS.add(edge[0]);
	    vS.add(edge[1]);
	}

        UnionFind uf = new UnionFind(vS.size());

	List<int[]> out = new ArrayList<>();
	
	if (outputSum != null && outputSum.length > 0) {
	    outputSum[0] = 0;
	}

	for (int i = 0; i < sortedIdxs.length; ++i) {

	    if (out.size() == (vS.size()-1)) break;

	    int idx = sortedIdxs[i];

	    if (uf.find(edges[idx][0]) != uf.find(edges[idx][1])) {
		uf.union(edges[idx][0], edges[idx][1]);

		out.add(Arrays.copyOf(edges[idx], 2));
		outputSum[0] += weights[idx];
	    }
	}

	return out;
    }
}
