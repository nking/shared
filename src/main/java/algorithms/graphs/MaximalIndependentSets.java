package algorithms.graphs;

import algorithms.misc.Misc0;
import algorithms.sort.MiscSorter;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 *
 * @author nichole
 */
public class MaximalIndependentSets {
    
    /**
     * Find a maximal independent set in graph G defined by the given adjacency map
     * using 
     * 
     * 
     * In graph theory, a maximal independent set (MIS) or maximal stable set is 
     * an independent set that is not a subset of any other independent set. 
     * In other words, there is no vertex outside the independent set that may 
     * join it because it is maximal with respect to the independent set property.
     * An independent set, stable set, coclique or anticlique is a set of 
     * vertices in a graph, no two of which are adjacent.
     * 
     * r.t.c. O( log_2(n)^2 )
     * <pre>
      References:
     
      Lecture Notes on a Parallel Algorithm for Generating a Maximal Independent Set
        Eric Vigoda
        Georgia Institute of Technology
        Last updated for 7530 - Randomized Algorithms, Spring 2010.
        https://www.cc.gatech.edu/~vigoda/7530-Spring10/MIS.pdf
        
      https://en.wikipedia.org/wiki/Maximal_independent_set#Listing_all_maximal_independent_sets
      </pre>
     * @param adjMap  adjacency map of undirected graph G with key = vertex u and value = set of
     * neighbors connected to u. 
     * @return returns a maximal independent set of undirected graph G defined by adjMap
     */
    public static TIntSet findOneProbabilistic(TIntObjectMap<TIntSet> adjMap){
        
        /*
        Lecture Notes on a Parallel Algorithm for Generating a Maximal Independent Set
        Eric Vigoda
        Georgia Institute of Technology
        Last updated for 7530 - Randomized Algorithms, Spring 2010.
        
        1. I=∅,V′=V and G′=G. 
        2. While (V′ ̸= ∅) do:
           (a) Choose a random set of vertices S ⊂ V′ 
               by selecting each vertex v independently with 
               probability 1/(2dG′(v)) 
               where dG′(v) is the degree of v in the graph G′.
           (b) For every edge (u, v) ∈ E(G′) if both endpoints are in S 
               then remove the vertex of lower degree from S 
               (Break ties arbitrarily). Call this new set S′.
           (c) I = I union S′. Let V′ = V′\(S′union N_G′(S′))
               where N_G′(S′) are the neighbors of S'.
               Finally, let G′ be the induced subgraph on V′.
        3. Output I
        
        NOTE: choosing S' is where the algorithm can run in parallel
        see 1st paragrpah of page 2 of 
        https://www.cc.gatech.edu/~vigoda/7530-Spring10/MIS.pdf
        "in every round we find a set S which is an independent set. 
        Then we add S to our current independent set I, and we remove 
            S union N(S) from the current graph V′.
        To choose S in parallel, each vertex v independently adds 
        themselves to S with a well chosen probability p(v). 
        We want to avoid adding adjacent vertices to S. 
        Hence, we will prefer to add low degree vertices (which is 1/(2dG′(v))).
        But, if for some edge (u,v), both endpoints were added to S, 
        then we keep the higher degree vertex."
            
        how many vertexes to choose for S'?  
        we want to ensure that by removing S ∪ N(S) from the graph, 
        we remove a constant fraction of the edges
        */
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        System.out.println("seed=" + seed);
        System.out.flush();
        rand.setSeed(seed);
        
        TIntSet indep = new TIntHashSet();
        
        TIntObjectMap<TIntSet> toFromAdjMap = extractInOut(adjMap);
        
        TIntSet vP = extractAllVertices(adjMap);
        int nV = vP.size();
        TIntObjectMap<TIntSet> gP;
        double[] cdf;
        double[] pmf;
        TIntIntMap indexMap;
        int[] chosen;
        int u;
        int v;
        int k;
        TIntSet s = new TIntHashSet();
        
        TIntIterator iterV;
               
        int nhbr;
        TIntSet nhbrsSet;
        TIntIterator iterNhbr;
                
        int d;
        double pV;
        double pVMax;
        double r;
        
        // regarding the probability, can either:
        // (1) create a pmf with the discrete values for each vertex as 1/(2dG′(v)),
        //     then create a normalized cdf from the pmf, then choose k random
        //     numbers between 0 and 1 and use the cdf to determine the variate (=vertex).
        // (2) or use a single bernoulli trial for each vertex: draw a random number 
        //     and compare it to p and accept if less than.  p=1/(2dG′(v))
        //     if all vertexes have same p, would expect to choose |V'|*p vertexes from V'.
        //     the variance in that is p*(1-p).
        //     standard deviation = Math.sqrt( (|V'|-1)*variance ).
        //     can expect to choose  |V'|*p +- Math.sqrt( (|V'|-1)*variance ) for 1 stdev
        //     or for 3 sigma deviation: |V'|*p +- 3*Math.sqrt( (|V'|-1)*variance )
        /*     One of the unit tests is a graph of 8 vertices, each with degree 3,
               hence p = 0.1667 for each of them.
               expected number selected = 8*0.1667 = 1.33
                   st.dev = 0.99,  3*st.dev = 2.96
        
               Two rounds of random doubles in range [0,1] have no numbers < 0.1667
               but the first round has 1 that is > (1-0.1667)
               and the second round has 2 that are > (1-0.1667)
        
               If instead used p=1/(2dG′(v)),
               each p = 0.333
               expected number selected = 8*0.0.333 = 2.66
                   st.dev = 1.25,  3*st.dev = 3.74
               The same 2 rounds of random doubles in range [0,1] 
               have 2 numbers < 0.33 in first round and 3 numbers < 0.33 in second round,
        
               Looks like a better probability function for this would be p=1/(dG′(v))
        
           random         
           number
         r=0.3253 
         r=0.6086 
         r=0.6723 
         r=0.9792 
         r=0.6570 
         r=0.1839 
         r=0.6623 
         r=0.9061 
        --------------------------------------
         r=0.4460 
         r=0.2900 
         r=0.5891 
         r=0.2038 
         r=0.5726 
         r=0.4459 
         r=0.2682 
         r=0.8990 
        */        
       
        // The later method (2) appears to be what would result in removing
        // "a constant fraction of the edges"
        
        int nIter = 0;
        while (!vP.isEmpty()) {
            
            // creates subgraph of adjMap out of vertexes vP.  note that
            //   it also creates edges uv and vu as G is undirected, and
            //   in order to make the degree counting easier below.
            gP = subgraph(vP, adjMap);
            
            /*
            if implementing (1):
            // create discrete probabilities, pdf, 1/(dG′(v))
            pdf = createCDF(gP, indexMap);
            cdf = MiscMath0.cumulativeSum(pdf);
            double norm = cdf[cdf.length - 1];
            for (i = 0; i < cdf.length; ++i) {
                cdf[i] /= norm;
            }
            k = some constant fraction of nV            
            chosen = CDFRandomSelect.chooseKFromBinarySearch(cdf, k, rand);
            s.clear();
            for (i = 0; i < chosen.length; ++i) {                
                v = chosen[i];
                s.add(v);
            }
            */
            
            /*
            if implementing (2):
            */
            s.clear();
            iterV = vP.iterator();
            while (iterV.hasNext()) {
                v = iterV.next();
                d = gP.get(v).size();
                pV = 1./(2.*d);
                
                r = rand.nextDouble();
                System.out.printf("r=%.4f, pV=%.4f\n", r, pV);
                                
                if (r <= pV) {
                    s.add(v);
                }
            }
            
            // ---- look for E(G') in S. calling them edge conflicts ----
            resolveEdgeConflicts(s, gP);
           
            // I = I union S′.
            indep.addAll(s);
            
            //Let V′ = V′ \ (S′ union N_G′(S′)).
            iterV = s.iterator();
            while (iterV.hasNext()) {
                v = iterV.next();
                vP.remove(v);
                if (gP.containsKey(v)) {
                    nhbrsSet = gP.get(v);
                    iterNhbr = nhbrsSet.iterator();
                    while (iterNhbr.hasNext()) {
                        nhbr = iterNhbr.next();
                        vP.remove(nhbr);
                    }
                }
            }
            nIter++;
        }
        
        System.out.printf("|V|=%d, nIter=%d\n", nV, nIter);
        
        return indep;
    }
    
    /**
     * Find a maximal independent set in graph G defined by the given adjacency map
     * using a sequential algorithm.
     * In graph theory, a maximal independent set (MIS) or maximal stable set is 
     * an independent set that is not a subset of any other independent set. 
     * In other words, there is no vertex outside the independent set that may 
     * join it because it is maximal with respect to the independent set property.
     * An independent set, stable set, coclique or anticlique is a set of 
     * vertices in a graph, no two of which are adjacent.
     * 
     * runtime complexity is O(log_2(|E|).
     * <pre>
      References:
     
      The sequential algorithm in this:
          Lecture Notes on a Parallel Algorithm for Generating a Maximal Independent Set
            Eric Vigoda
            Georgia Institute of Technology
            Last updated for 7530 - Randomized Algorithms, Spring 2010.
            https://www.cc.gatech.edu/~vigoda/7530-Spring10/MIS.pdf

      https://en.wikipedia.org/wiki/Maximal_independent_set#Listing_all_maximal_independent_sets
      </pre>
     * @param adjMap and adjacency map with key = vertex u and value = set of
     * neighbors connected to u.
     * @return returns a maximal independent set of graph G defined by adjMap
     */
    public static TIntSet findOne(TIntObjectMap<TIntSet> adjMap){
        /*
        Lecture Notes on a Parallel Algorithm for Generating a Maximal Independent Set
        Eric Vigoda
        Georgia Institute of Technology
        Last updated for 7530 - Randomized Algorithms, Spring 2010.

                1. I = ∅, V′= V.
                2. While (V′ ̸= ∅) do
                    (a) Choose any v∈V′. 
                    (b) Set I = I ∪ v.
                    (c) Set V′ = V′ \ (v ∪ N(v)). where N(v) are the neighbors of v.
                3. Output I
        */
        
        TIntObjectMap<TIntSet> toFromAdjMap = extractInOut(adjMap);
        
        TIntSet indep = new TIntHashSet();
        TIntSet vP = extractAllVertices(adjMap);
        int v;
        TIntSet nhbrs;
        TIntIterator iter;
        while (!vP.isEmpty()) {
            
            v = vP.iterator().next();
            
            indep.add(v);
            
            vP.remove(v);
            
            nhbrs = toFromAdjMap.get(v);
            if (nhbrs != null) {
                iter = nhbrs.iterator();
                while (iter.hasNext()) {
                    vP.remove(iter.next());
                }
            }
        }
        
        return indep;
    }

    static TIntObjectMap<TIntSet> extractInOut(TIntObjectMap<TIntSet> adjMap) {
        
        TIntObjectMap<TIntSet> out = new TIntObjectHashMap<TIntSet>();
        
        TIntObjectIterator<TIntSet> iter = adjMap.iterator();
        TIntSet set, outSetI, outSetJ;
        int i, j;
        TIntIterator iter2;
        
        while (iter.hasNext()) {
            iter.advance();
            i = iter.key();
            set = iter.value();
            
            if (set == null) {
                continue;
            }
            
            outSetI = out.get(i);
            if (outSetI == null) {
                outSetI = new TIntHashSet();
                out.put(i, outSetI);
            }
            
            iter2 = set.iterator();
            while (iter2.hasNext()) {
                j = iter2.next();

                outSetJ = out.get(j);
                if (outSetJ == null) {
                    outSetJ = new TIntHashSet();
                    out.put(j, outSetJ);
                }

                outSetI.add(j);
                outSetJ.add(i);
            }
        }
        return out;
    }

    static TIntSet extractAllVertices(TIntObjectMap<TIntSet> adjMap) {
        TIntSet out = new TIntHashSet();
        
        TIntObjectIterator<TIntSet> iter = adjMap.iterator();
        TIntIterator iter2;
        TIntSet set;
        while (iter.hasNext()) {
            iter.advance();
            out.add(iter.key());
            set = iter.value();
            if (set == null) {
                continue;
            }
            iter2 = set.iterator();
            while (iter2.hasNext()) {
                out.add(iter2.next());
            }
        }
        
        return out;
    }

    static void sort(List<PairInt> edgeConflicts, TIntList edgeDegreesLow) {
        
        int n = edgeConflicts.size();
        
        if (n != edgeDegreesLow.size()){
            throw new IllegalArgumentException("arguments must be the same size");
        }
        
        int[] indexes = new int[n];
        int i;
        for(i = 0; i < n; ++i) {
            indexes[i] = i;
        }
        
        int[] a = Arrays.copyOf(edgeDegreesLow.toArray(), n);
        
        MiscSorter.sortBy1stArg(a, indexes);
        
        PairInt[] s1 = new PairInt[n];
        for (i = 0; i < n; ++i) {
            s1[i] = edgeConflicts.get(indexes[i]);
        }
        edgeConflicts.clear();
        edgeDegreesLow.clear();
        for (i = 0; i < n; ++i) {
            edgeConflicts.add(s1[i]);
            edgeDegreesLow.add(a[i]);
        }
        
    }

    /**
     * creates subgraph of adjMap out of vertexes vP.  note that
       it also creates edges uv and vu as G is undirected, and
       in order to make the degree counting easier where needed.
     * @param vP
     * @param adjMap
     * @return 
     */
    static TIntObjectMap<TIntSet> subgraph(TIntSet vP, TIntObjectMap<TIntSet> adjMap) {
        
        TIntObjectMap<TIntSet> out = new TIntObjectHashMap<TIntSet>();
        
        TIntIterator iter;
        int v, adj;
        TIntIterator iter2;
        TIntSet outSet;
        TIntSet adjSet;
        TIntSet outSet2;
        
        iter = vP.iterator();
        while (iter.hasNext()) {
            v = iter.next();
            outSet = out.get(v);
            if (outSet == null) {
                outSet = new TIntHashSet();
                out.put(v, outSet);
            }
            
            adjSet = adjMap.get(v);
            
            if (adjSet == null || adjSet.isEmpty()) {
                continue;
            }
            
            iter2 = adjSet.iterator();
            while (iter2.hasNext()) {
                adj = iter2.next();
                
                outSet.add(adj);
                
                outSet2 = out.get(adj);
                if (outSet2 == null) {
                    outSet2 = new TIntHashSet();
                    out.put(adj, outSet2);
                }
                
                outSet2.add(v);
            }
        }
        
        return out;
    }

    private static void resolveEdgeConflicts(TIntSet s, TIntObjectMap<TIntSet> gP) {
         
        //For every edge (u, v) ∈ E(G′) if both endpoints are in S then 
        // remove the vertex of lower degree from S (Break ties arbitrarily). 
        // Call this new set S′.

        // avoiding the comparison of all combinations of vertexes in edgeConflicts
        //     to find the max size S, by a greedy pick. 

        // these will be ordered by lowest degree, to remove a vertex in each edge from S
        List<PairInt> edgeConflicts = new ArrayList<PairInt>();
        TIntList edgeDegreesLow = new TIntArrayList();

        int dU;
        TIntObjectIterator<TIntSet> iterGP;
        int v;
        int nhbr;
        TIntSet nhbrsSet;
        TIntIterator iterNhbr;

        iterGP = gP.iterator();
        while (iterGP.hasNext()) {
            iterGP.advance();

            v = iterGP.key();
            nhbrsSet = iterGP.value();

            if (!s.contains(v) || nhbrsSet == null || nhbrsSet.isEmpty()) {
                continue;
            }

            iterNhbr = nhbrsSet.iterator();
            while (iterNhbr.hasNext()) {

                nhbr = iterNhbr.next();

                if (!s.contains(nhbr)) {
                    continue;
                }

                dU = gP.get(nhbr).size();
                if (v < nhbr) {
                    edgeConflicts.add(new PairInt(v, nhbr));
                } else {
                    edgeConflicts.add(new PairInt(nhbr, v));
                }
                edgeDegreesLow.add(Math.min(nhbrsSet.size(), dU));
            }
        }

        sort(edgeConflicts, edgeDegreesLow);

        int u;
        int d;
        
        for (PairInt uv : edgeConflicts) {
            u = uv.getX();
            v = uv.getY();
            //for each conflicting edge, check again for presence in S and then keep only the higher degree vertex
            if (s.contains(u) && s.contains(v)) {
                d = gP.get(v).size();
                dU = gP.get(u).size();
                if (dU > d) {
                    s.remove(v);
                } else {
                    s.remove(u);
                }
            }
        }
                    
    }
}
