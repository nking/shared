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
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * In graph theory, a k-degenerate graph is an undirected graph in which every 
 * subgraph has a vertex of degree at most k: that is, some vertex in the 
 * subgraph touches k or fewer of the subgraph's edges. The degeneracy of a 
 * graph is the smallest value of k for which it is k-degenerate. 
 * The degeneracy of a graph is a measure of how sparse it is, and is within 
 * a constant factor of other sparsity measures such as the arboricity of a 
 * graph.
 <pre>
 https://en.wikipedia.org/wiki/Degeneracy_(graph_theory)#Relation_to_other_graph_parameters
</pre>

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
 * 
 * <pre>
     the runtime is O( (|V| + |E|) * log_2(w)) where w is maxDegree
     which is essentially O((|V| + |E|))
 * </pre>
 * 
 * NOTE:  for an undirected graph that is a simple undirected graph not a multigraph, 
        the maximum possible number of edges |E|_max = C(|V|, 2) 
        and the max possible degree of a vertex = |V|-1,
        and in that case |V| will not be less than a number whose maximum is |V|-1.
        
  This code does not at this time handle multigraphs, but the ability to store the 
   structure needed for it is present.  When the multigraph datastructure is used,
   the runtime complexity is O( (|V| + |E|) * log_2(|V|))
 
 * @author nichole
 */
public class DegeneracyOrderingMatulaBeck {
    
    /**
     * following
     * https://en.wikipedia.org/wiki/Degeneracy_(graph_theory)#Relation_to_other_graph_parameters
     * 
      <pre>
        the runtime is O( (|V| + |E|) * log_2(w)) where w is maxDegree
          which is essentially O((|V| + |E|))
        
      </pre>
 
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
        
        IBucketQueue bQ;
        
        int maxDegree = calculateMaxDegree(adjMap);
        
        System.out.printf("maxDegree=%d, n=%d\n", maxDegree, n);
        
        // use a YFastTrie and HashMap to create the bucket queue
        // runtime complexity is O(|V|*log_2(w)) where q=maxDegree
        bQ = new TBucketQueue(adjMap);
        
        TIntSet bucket, adj;
        
        // initialize output vertex list to -1 for debugging
        Arrays.fill(out, -1);
        // keep a set of values of out
        TIntSet outSet = new TIntHashSet(n);
        
        TIntIterator iter2, iter3;
        
        // while loop: runtime complexity of each iteration is O(log_2(maxDegree)) or O(log_2(n))
        // and the loop visits each vertex and it's neighbors which are not
        // yet in the output set = |V| + |E| iterations.
        // if using a SBucketQueue, the runtime complexity is O( (|V| + |E|) * log_2(|V|))
        // if using a TBucketQueue, the runtime complexity is O( (|V| + |E|) * log_2(w)) where w is maxDegree
        
        // the degree k
        int v, k = 0, j = 0, w, nW, x;
        while (j < n && !bQ.isEmpty()) {
            
            // this handles check for null or empty bucket.
            //   note that the internal bucketMap size is out of sync w/ bucketQueue size until remove is invoked
            i = bQ.extractMinimumPt1();
            
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
    
    protected static int calculateMaxDegree(TIntObjectMap<TIntSet> adjMap) {
        
        int n = adjMap.size();

        TIntSet adj;
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
            if (nV > maxDegree) {
                maxDegree = nV;
            }
        }
        return maxDegree;
    }
   
    static class TBucketQueue implements IBucketQueue {
        
        // key = degree of a vertex, value = set of vertexes with degree given by key.
        protected final TIntObjectMap<TIntSet> bucketMap = new TIntObjectHashMap<TIntSet>();
        
        // needed in the move operation for neighbors of v
        // key = vertex index,  value = degree of vertex (== # of neighbors)
        protected final TIntIntMap reverseBucketMap = new TIntIntHashMap();
        
        // holds the degrees of vertexes, uniquely, and in a datastructure with
        // operations successor, predecessor, min, max, etc.
        protected final YFastTrie bucketQueue;
        
        public TBucketQueue(TIntObjectMap<TIntSet> adjMap) {
        
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
        
        @Override
        public boolean isEmpty() {
            return (bucketQueue.size() == 0);
        }

        @Override
        public int minimum() {
            return bucketQueue.minimum();
        }
        
        @Override
        public int size() {
            return bucketQueue.size();
        }
        
        /**
         * extract the minimum from the internal priority queue called bucketQueue.
         * note that the internal bucketMap size is out of sync w/ bucketQueue 
         * size until remove is invoked.
         * @return 
         */
        @Override
        public int extractMinimumPt1() {
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
        @Override
        public boolean contains(int i) {
            return bucketMap.containsKey(i);
        }
        
        /**
         * find whether the bucket queue contains the vertex v.
         * @param v
         * @return 
         */
        @Override
        public boolean containsVertex(int v) {
            return reverseBucketMap.containsKey(v);
        }

        @Override
        public TIntSet getBucket(int i) {            
            return bucketMap.get(i);
        }
       
        /**
         * remove value v from the bucket queue for key k.  this method also
         * removes key entries in bucketMap and bucketQueue if the bucket empties.
         * @param k
         * @param v 
         */
        @Override
        public void remove(int k, int v) {
            TIntSet bucket = getBucket(k);
            remove(k, v, bucket);
        }
        
        @Override
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
        @Override
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
                bucket = new TIntHashSet();
                bucketMap.put(toDBucket, bucket);
                bucketQueue.add(toDBucket);
            } else {
                assert(bucketQueue.find(toDBucket) != -1);
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
        @Override
        public void moveItem(int item, int toDBucket) {
            int fromDBucket = reverseBucketMap.get(item);
            moveItem(item, fromDBucket, toDBucket);
        }
        
    }
    
    static class SBucketQueue implements IBucketQueue {
        
        /* 
        This implementation is currently unused, but remains present in case
        a method for multigraphs is ever made.
        
        this datastructure is used in findDegeneracyOrder() in a single iterating
        clause, and so the limiting runtime complexity is the largest of the
        bucket runtime operations.
        
        a SortedSet is the main bucket queue choice when n < maxValue.
        
        Here are the IBucketQueue operation invoked in findDegeneracyOrder():
        <pre>
        while (j < n && !bQ.isEmpty()) {
           extractMinimum();
           bQ.getBucket();
           bucket.iterator().next();
           bQ.remove();
           bQ.moveItem(w, nW);
          ...
        }
        </pre>
        If a TIntObjectMap<TIntSet> bucketMap were used as it is in TBucketQueue,
        its fast O(1) getBucket does not affect the runtime complexity, and so
        the values it would have held can be placed in a SortedMap instead.
        */
        
        // needed in the move operation for neighbors of v
        // key = vertex index,  value = degree of vertex (== # of neighbors)
        protected final TIntIntMap reverseBucketMap = new TIntIntHashMap();
        
        // holds the degrees of vertexes, uniquely, and in a datastructure with
        // operations successor, predecessor, min, max, etc.
        protected final SortedMap<Integer, TIntSet> bucketQueue = new TreeMap<>();
        
        public SBucketQueue(TIntObjectMap<TIntSet> adjMap) {
        
            int n = adjMap.size();
            
            TIntSet adj;
            TIntSet bucket;
            TIntObjectIterator<TIntSet> iter = adjMap.iterator();
            int i;
            int v;
            int nV;
            int maxDegree = 0;
            for (i = 0; i < n; ++i) {
                iter.advance();
                v = iter.key();
                adj = iter.value();
                if (adj == null) {
                    nV = 0;
                } else {
                    nV = adj.size();
                }
                bucket = bucketQueue.get(nV);
                if (bucket == null) {
                    bucket = new TIntHashSet();
                    bucketQueue.put(nV, bucket);
                }
                bucket.add(v);
                if (nV > maxDegree) {
                    maxDegree = nV;
                }
                reverseBucketMap.put(v, nV);
                
                //System.out.printf("bM key=%d, val=%d, max=%d\n", nV, v, maxDegree);
            }
        }
        
        @Override
        public boolean isEmpty() {
            return bucketQueue.isEmpty();
        }

        @Override
        public int minimum() {
            return bucketQueue.firstKey();
        }
        
        @Override
        public int size() {
            return bucketQueue.size();
        }
        
        /**
         * this is purely a minimum() invocation.  The invoker uses remove 
         * afterwards to remove the entry from the bucket values.
         * @return 
         */
        @Override
        public int extractMinimumPt1() {
            return minimum();
        }

        /**
         * find whether the bucket queue contains the degree key i.
         * @param i
         * @return 
         */
        @Override
        public boolean contains(int i) {
            return bucketQueue.containsKey(i);
        }
        
        /**
         * find whether the bucket queue contains the vertex v.
         * @param v
         * @return 
         */
        @Override
        public boolean containsVertex(int v) {
            return reverseBucketMap.containsKey(v);
        }

        @Override
        public TIntSet getBucket(int i) {            
            return bucketQueue.get(i);
        }
       
        /**
         * remove value v from the bucket queue for key k.  this method also
         * removes key entries in bucketMap and bucketQueue if the bucket empties.
         * @param k
         * @param v 
         */
        @Override
        public void remove(int k, int v) {
            TIntSet bucket = getBucket(k);
            remove(k, v, bucket);
        }
        
        @Override
        public void remove(int k, int v, TIntSet bucket) {
            if (bucket == null || bucket.isEmpty()) {
                throw new IllegalStateException("programming error: bucket is null or empty");
            }
            bucket.remove(v);
            if (bucket.isEmpty()) {
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
        @Override
        public void moveItem(int item, int fromDBucket, int toDBucket) {
            
            if (fromDBucket == toDBucket) {
                throw new IllegalArgumentException("fromDBucket cannot equal toDBucket");
            }
            
            TIntSet bucket = bucketQueue.get(fromDBucket);
            
            if (bucket == null || !bucket.contains(item)) {
                throw new IllegalStateException("bucket " + fromDBucket + 
                    " is empty, so cannot remove " + item + " from it");
            }
            
            // remove from current bucket 
            bucket.remove(item);
            if (bucket.isEmpty()) {
                // remove bucket from bucketQueue
                bucketQueue.remove(fromDBucket);
            }
            
            // add to new bucket queue
            bucket = bucketQueue.get(toDBucket);
            if (bucket == null) {
                bucket = new TIntHashSet();
                bucketQueue.put(toDBucket, bucket);
            }
            assert(!bucketQueue.get(toDBucket).contains(item));
            bucket.add(item);
            reverseBucketMap.put(item, toDBucket);
        }
        
        /**
         * remove item from current bucket that its in, and move it to toDBucket.
         * @param item
         * @param toDBucket the bucket degree number to add item to
         */
        @Override
        public void moveItem(int item, int toDBucket) {
            int fromDBucket = reverseBucketMap.get(item);
            moveItem(item, fromDBucket, toDBucket);
        }
        
    }
    
    static interface IBucketQueue {

        /**
         * find whether the bucket queue contains the degree key i.
         *
         * @param i
         * @return
         */
        boolean contains(int i);

        /**
         * find whether the bucket queue contains the vertex v.
         *
         * @param v
         * @return
         */
        boolean containsVertex(int v);

        /**
         * extract the minimum from the internal priority queue called
         * bucketQueue. note that the internal bucketMap size is out of sync w/
         * bucketQueue size until remove is invoked.
         *
         * @return
         */
        int extractMinimumPt1();

        TIntSet getBucket(int i);

        boolean isEmpty();

        int minimum();

        /**
         *
         * @param item
         * @param fromDBucket the bucket degree number to remove item from
         * @param toDBucket the bucket degree number to add item to
         */
        void moveItem(int item, int fromDBucket, int toDBucket);

        /**
         * remove item from current bucket that its in, and move it to
         * toDBucket.
         *
         * @param item
         * @param toDBucket the bucket degree number to add item to
         */
        void moveItem(int item, int toDBucket);

        /**
         * remove value v from the bucket queue for key k. this method also
         * removes key entries in bucketMap and bucketQueue if the bucket
         * empties.
         *
         * @param k
         * @param v
         */
        void remove(int k, int v);

        void remove(int k, int v, TIntSet bucket);

        int size();

    }
}
