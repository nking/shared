package algorithms.msts;

import algorithms.disjointSets.DisjointForest;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.sort.MiscSorter;
import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TObjectDoubleIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TIntObjectHashMap;

/**
*
  Minimum spanning tree is the minimal network that spans all nodes in a tree
  and has the smallest cost (sum of edges).
  
  Kruskal's grows a forest by sorting the edges first and then adding edges that 
  are not yet connected to the tree.
  
  * RuntimeComplexity is O(|E| lg_2|E|), 
      which for sparse graphs having |E| .lt. |V|^2,
      gives O(|E| lg_2|V|).
      
  Best time Kruskal's is O(|E| lg_2|V|).
  Best time  Prim's w/ fib heaps is O(|E| + |V|lg_2|V|).
  
 * @author nichole
 */
public class KruskalsMinimumSpanningTree {
    
    /**
     * find a minimum spanning tree using Kruskal's algorithm.
     * @param graph an adjacency list for the graph, where the index of the array
     * is the vertex number and the keys within each vertex's list node are the
     * vertex on the other end of an edge.
     * @param edgeWeights map of edge weights with key being pairs of the
     * vertex numbers.
     * @return a minimum spanning tree of the weighted given graph
     */
    public static TIntObjectMap<SimpleLinkedListNode> mst(
        SimpleLinkedListNode[] graph, TObjectDoubleMap<PairInt> edgeWeights) {
                
        /* pseudo code form Cormen et al. "Introduction to Computer Science".
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
    
}
