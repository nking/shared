package algorithms.graphs;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.List;

/**
 A set cover is the smallest vertex cover of a graph,
 where a vertex cover is is a subset of a graph's vertices which represents at 
 * least one vertex from every edge in the full graph
 * 
 * optimal vertex cover is np-hard (non-deterministic polynomial class
   problems at least as hard as the hardest problems in NP.  No known polynomial
   time algorithm, but one can guess a single solution and verify it.)    

 * @author nichole
 */
public class SetCover {
    
    /**
     * solve or the minimum weighted set cover using an approximation algorithm
     * of 2*log_2(n) where n is the number of vertexes in the final
     * cover (== the unique number in all of sets).
     * @param x the items which must be present in the final cover
     * @param sets
     * @param weights
     * @return the list of indexes of sets which comprise the cover.
     */
    public TIntList weightedApprox2LgN(TIntSet x, List<TIntSet> sets, double[] weights) {
        /*
        material from ecture slides of Principal lecturer: Dr Thomas Sauerwald
        Advanced Algorithms, University of Cambridge.
        VII. Approximation Algorithms: Randomisation and Rounding
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/materials.html
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/rand.pdf
        
        for the linear program:
            minimize: 
                summation_S_in_Cover( c(S) ) = summation_S_in_sets( c(S)*y(S) )
            subject to:
                summation_S_in_sets : x in S ( y(S) ) >= 1
                y(S) <= 1 for each S in sets ----\
            non-negativity constraints:           \ these 2 rules are derived from y(S) ∈ [0,1]
                y(S) >= 0 for each S in sets ----/
        
        for the weighted set cover w/ LP(X, F, c):
            // where each x in X belongs to at least 1 subset in F
            // F is the list of subsets to choose from when building the cover
            // c is the cost for each set in F
            compute y, an optimal solution to the linear program
            C = empty set
            repeat 2*log2(n) times
                for each S in F
                    let C = C union S with probabilty y(S)
            return C
        */
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    /**
     * find a set cover that is O(log-n)-approx, that is no more than O(log-n) times
     * as large as the optimal set cover
     * (e.g., for n=100, this would be up to 6.64 times as large as the optimal solution).
     * This is a greedy approach.
     * The runtime complexity is polynomial.
     * <pre>
     * The algorithm implements pseudocode from
     *  from https://www.ics.uci.edu/~goodrich/teach/graph/notes/Approximation.pdf
     * and Cormen et al "Introduction to Algorithms" chap 35.3.
     * </pre>
     * @param sets
     * @return the list of indexes of sets which comprise the cover.
     */
    public TIntList approxLgN(List<TIntSet> sets) {
        
        // make a copy of the sets to edit it
        TIntObjectMap<TIntSet> setsMap = copy(sets);
        TIntObjectIterator<TIntSet> iter;
        
        TIntList c = new TIntArrayList();
        
        TIntSet si, siMaxN = null;        
        int maxN, i, n, maxNIdx = -1, idx;
        TIntList rm = new TIntArrayList();
        while (!setsMap.isEmpty()) {
            maxN = Integer.MIN_VALUE;
            iter = setsMap.iterator();
            for (i = 0; i < setsMap.size(); ++i) {
                iter.advance();
                si = iter.value();
                idx = iter.key();
                n = si.size();
                //System.out.printf("   idx=%d (%s) n=%d\n", idx, si.toString(), n);
                if (n > 0 && n > maxN) {
                    maxN = n;
                    maxNIdx = idx;
                }
            }
            if (maxN == Integer.MIN_VALUE) {
                break;
            }
            siMaxN = setsMap.remove(maxNIdx);
            //System.out.printf("max idx=%d (%s)\n", maxNIdx, siMaxN.toString());
            c.add(maxNIdx);
            
            // remove siMaxN from each set in s.  and store the empty sets as indexes to remove after use of iterator
            iter = setsMap.iterator();
            for (i = 0; i < setsMap.size(); ++i) {
                iter.advance();
                si = iter.value();
                si.removeAll(siMaxN);
                if (si.isEmpty()) {
                    rm.add(iter.key());
                }
            }
            for (i = 0; i < rm.size(); ++i) {
                setsMap.remove(rm.get(i));
            }
            rm.clear();
        }
        
        return c;
    }

    protected TIntObjectMap<TIntSet> copy(List<TIntSet> s) {
        TIntObjectMap<TIntSet> c = new TIntObjectHashMap<TIntSet>(s.size());
        int i;
        TIntSet cSet, si;
        TIntIterator iter;
        for (i = 0; i < s.size(); ++i) {
            cSet = new TIntHashSet();
            si = s.get(i);
            iter = si.iterator();
            while (iter.hasNext()){
                cSet.add(iter.next());
            }
            c.put(i, cSet);
        }
        return c;
    }

}
