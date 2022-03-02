package algorithms.graphs;

import algorithms.heapsAndPQs.YFastTrie;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;

/**
 *In graph theory, a k-degenerate graph is an undirected graph in which every 
 * subgraph has a vertex of degree at most k: that is, some vertex in the 
 * subgraph touches k or fewer of the subgraph's edges. The degeneracy of a 
 * graph is the smallest value of k for which it is k-degenerate. 
 * The degeneracy of a graph is a measure of how sparse it is, and is within 
 * a constant factor of other sparsity measures such as the arboricity of a 
 * graph.

Degeneracy is also known as the k-core number, width, and linkage, and 
* is essentially the same as the coloring number or Szekeres-Wilf number 
* (named after Szekeres and Wilf (1968)). k-degenerate graphs have also been 
* called k-inductive graphs. The degeneracy of a graph may be computed in 
* linear time by an algorithm that repeatedly removes minimum-degree vertices.
* The connected components that are left after all vertices of degree less than 
* k have been removed are called the k-cores of the graph and the degeneracy of 
* a graph is the largest value k such that it has a k-core.
* 
 * As Matula & Beck (1983) describe, it is possible to find a vertex ordering 
 * of a finite graph G that optimizes the coloring number of the ordering, in 
 * linear time, by using a bucket queue to repeatedly find and remove the vertex 
 * of smallest degree. 
 * <pre>
 * https://en.wikipedia.org/wiki/Degeneracy_(graph_theory)#Relation_to_other_graph_parameters
 * </pre>
 * @author nichole
 */
public class DegeneracyOrderingMatulaBeck {
    
    /**
     * following
     * https://en.wikipedia.org/wiki/Degeneracy_(graph_theory)#Relation_to_other_graph_parameters
     * 
     * runtime complexity is O(|V| + |E|).
     * @param adjMap an undirected graph as an adjacency map with key=vertex index, value =
     * set of indexes of vertexes for the neighbors of the key.  note: the method
     * expects that the adjacency map has a key for every vertex listed in the values.
     * @param out output array of size |V|
     * @return the degeneracy of the graph.
     */
    public static int findDegeneracyOrder(TIntObjectMap<TIntSet> adjMap, int[] out) {
                
        int n = adjMap.size(), i;
        
        if (out.length != n) {
            throw new IllegalArgumentException("expecting the the number of vertexes"
                + " is adjMap.size() and that out.length==adjMap.size()");
        }
        
        // making a bucket queue with a YFastTrie and a HashMap        
        // where YFastTrie key = the degree of a vertex.  and the bucketMap holds key=degree, value=vertex indexes
        // (YFastTrie needs to know the max degree to set word size).
        
        // ---- initialize bucket ------
        BucketQueue bQ = new BucketQueue(adjMap);
        
        TIntSet bucket, adj;
        
        // initialize output vertex list to -1 for debugging
        Arrays.fill(out, -1);
        // keep a set of values of out
        TIntSet outSet = new TIntHashSet(n);
        
        TIntIterator iter2, iter3;
        
        // the degree k
        int v, k = 0, j = 0, w, nW, x;
        while (j < n && !bQ.isEmpty()) {
            
            // this handles check for null or empty bucket.
            //   note that the internal bucketMap size is out of sync w/ bucketQueue size until remove is invoked
            i = bQ.extractMinimum();
            
            k = Math.max(k, i);
            
            assert(!outSet.contains(i));
            
            bucket = bQ.getBucket(i);
            
            v = bucket.iterator().next();
            // add V to out and remove from bucket in bucketMap.  bucketQueue has been updated already.
            outSet.add(v);
            out[n-j-1] = v;
            bQ.remove(i, v, bucket);
            
            //System.out.printf("out=%s\n", Arrays.toString(out));
              
            // for each neighbor w of v that is not already in L
            adj = adjMap.get(v);
            if (adj == null || adj.isEmpty()) {
                j++;
                continue;
            }
            
            iter2 = adj.iterator();
            while (iter2.hasNext()) {
                w = iter2.next();
                if (outSet.contains(w)) {
                    continue;
                }
                // user should have created adjMap with a key for every vertex
                assert(adjMap.containsKey(w));
                // count neighbors of w not already in L
                nW = 0;
                if (adjMap.get(w).size() > 0) {
                    iter3 = adjMap.get(w).iterator();
                    while (iter3.hasNext()) {
                        x = iter3.next();
                        if (!outSet.contains(x)) {
                            nW++;
                        }
                    }
                }
                nW--;
                if (nW < 0) {
                    //System.out.printf("vertex w=%d has nW=%d.  bucket size=%d\n",
                    //    w, nW, bQ.size());
                    //w should remain in bQ until it becomes v
                } else {
                    // move w from current bucket to bucket dW if dW>=0
                    bQ.moveItem(w, nW);
                }
                
            }
            j++;
        }
        return k;
    }
    
    static class BucketQueue {
        
        // key = degree of a vertex, value = set of vertexes with degree given by key.
        protected final TIntObjectMap<TIntSet> bucketMap = new TIntObjectHashMap<TIntSet>();
        
        // needed in the move operation for neighbors of v
        // key = vertex index,  value = degree of vertex (== # of neighbors)
        protected final TIntIntMap reverseBucketMap = new TIntIntHashMap();
        
        // holds the degrees of vertexes, uniquely, and in a datastructure with
        // operations successor, predecessor, min, max, etc.
        // NOTE that the operations are O(log_2(w)) where w is the bitlength of the
        // largest value to be inserted into the trie.
        // If the number of items to be inserted, n, is small such that log_2(n) < log_2(w),
        // then a SortedSet would be a better datastructure here for the bucketQueue.
        // TODO: add that check and ability to use either data structure here.
        protected final YFastTrie bucketQueue;
        
        public BucketQueue(TIntObjectMap<TIntSet> adjMap) {
        
            int n = adjMap.size();
            
            TIntSet bucket, adj;
            TIntObjectIterator<TIntSet> iter = adjMap.iterator();
            int i, v, nV, maxDegree = 0;
            for (i = 0; i < n; ++i) {
                iter.advance();
                v = iter.key();
                adj = iter.value();
                if (adj == null) {
                    nV = 0;
                } else {
                    nV = adj.size();
                }
                bucket = bucketMap.get(nV);
                if (bucket == null) {
                    bucket = new TIntHashSet();
                    bucketMap.put(nV, bucket);
                }
                bucket.add(v);
                if (nV > maxDegree) {
                    maxDegree = nV;
                }
                reverseBucketMap.put(v, nV);
                
                //System.out.printf("bM key=%d, val=%d, max=%d\n", nV, v, maxDegree);
            }
            
            bucketQueue = new YFastTrie(1 + (int)Math.ceil(Math.log(maxDegree)/Math.log(2)) );
            iter = bucketMap.iterator();
            for (i = 0; i < bucketMap.size(); ++i) {
                iter.advance();
                nV = iter.key();
                bucketQueue.add(nV);
                //System.out.printf("bQ key=%d\n", nV);
            }
        }
        
        public boolean isEmpty() {
            return (bucketQueue.size() == 0);
        }

        public int minimum() {
            return bucketQueue.minimum();
        }
        
        public int size() {
            return bucketQueue.size();
        }
        
        /**
         * extract the minimum from the internal priority queue called bucketQueue.
         * note that the internal bucketMap size is out of sync w/ bucketQueue 
         * size until remove is invoked.
         * @return 
         */
        public int extractMinimum() {
            if (bucketQueue.size() == 0) {
                throw new IllegalStateException("bucket queue is empty");
            }
            int min = bucketQueue.minimum();
            TIntSet bucket = bucketMap.get(min);
            if (bucket == null || bucket.isEmpty()) {
                throw new IllegalStateException("programming error: bucket is null or empty");
            }
            if (bucket.size() > 1) {
                return min;
            }
            return bucketQueue.extractMinimum();
        }

        /**
         * find whether the bucket queue contains the degree key i.
         * @param i
         * @return 
         */
        public boolean contains(int i) {
            return bucketMap.containsKey(i);
        }
        
        /**
         * find whether the bucket queue contains the vertex v.
         * @param v
         * @return 
         */
        public boolean containsVertex(int v) {
            return reverseBucketMap.containsKey(v);
        }

        public TIntSet getBucket(int i) {            
            return bucketMap.get(i);
        }
       
        /**
         * remove value v from the bucket queue for key k.  this method also
         * removes key entries in bucketMap and bucketQueue if the bucket empties.
         * @param k
         * @param v 
         */
        public void remove(int k, int v) {
            TIntSet bucket = getBucket(k);
            remove(k, v, bucket);
        }
        
        public void remove(int k, int v, TIntSet bucket) {
            if (bucket == null || bucket.isEmpty()) {
                throw new IllegalStateException("programming error: bucket is null or empty");
            }
            bucket.remove(v);
            if (bucket.isEmpty()) {
                bucketMap.remove(k);
                bucketQueue.remove(k);
                reverseBucketMap.remove(v);
            }
        }

        /**
         * 
         * @param item
         * @param fromDBucket the bucket degree number to remove item from
         * @param toDBucket the bucket degree number to add item to
         */
        public void moveItem(int item, int fromDBucket, int toDBucket) {
            
            if (fromDBucket == toDBucket) {
                throw new IllegalArgumentException("fromDBucket cannot equal toDBucket");
            }
            
            TIntSet bucket = bucketMap.get(fromDBucket);
            
            if (bucket == null || !bucket.contains(item)) {
                throw new IllegalStateException("bucket " + fromDBucket + 
                    " is empty, so cannot remove " + item + " from it");
            }
            
            // remove from current bucket 
            bucket.remove(item);
            if (bucket.isEmpty()) {
                // remove bucket from bucketQueue
                bucketQueue.remove(fromDBucket);
                bucketMap.remove(fromDBucket);
            }
            
            // add to new bucket queue
            bucket = bucketMap.get(toDBucket);
            if (bucket == null) {
                assert(bucketQueue.find(toDBucket) == -1);
                if (bucket == null) {
                    bucket = new TIntHashSet();
                    bucketMap.put(toDBucket, bucket);
                }
                bucketQueue.add(toDBucket);
            }
            bucket.add(item);
            reverseBucketMap.put(item, toDBucket);
            assert(bucketQueue.find(toDBucket) > -1);
        }
        
        /**
         * remove item from current bucket that its in, and move it to toDBucket.
         * @param item
         * @param toDBucket the bucket degree number to add item to
         */
        public void moveItem(int item, int toDBucket) {
            int fromDBucket = reverseBucketMap.get(item);
            moveItem(item, fromDBucket, toDBucket);
        }
        
    }
}
