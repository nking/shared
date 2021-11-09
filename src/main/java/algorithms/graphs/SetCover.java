package algorithms.graphs;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

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
     * find a set cover that is O(log-n)-approx, that is no more than O(log-n) times
     * as large as the optimal set cover
     * (e.g., for n=100, this would be up to 6.64 times as large as the optimal solution).
     * This is a greedy approach.
     * <pre>
     * The algorithm implements pseudocode from
     *  from https://www.ics.uci.edu/~goodrich/teach/graph/notes/Approximation.pdf
     * </pre>
     */
    public TIntSet approxLgN(List<TIntSet> sets) {
        
        // make a copy of the sets to edit it
        List<TIntSet> s = copy(sets);
                
        TIntSet c = new TIntHashSet();
        
        TIntSet si, siMaxN = null;        
        int maxN, i, n, maxNIdx = -1;
        while (!s.isEmpty()) {
            maxN = Integer.MIN_VALUE;
            for (i = 0; i < s.size(); ++i) {
                si = s.get(i);
                n = si.size();
                if (n > 0 && n > maxN) {
                    maxN = n;
                    maxNIdx = i;
                }
            }
            if (maxN == Integer.MIN_VALUE) {
                break;
            }
            siMaxN = s.get(maxNIdx);
            c.addAll(siMaxN);
            s.remove(maxNIdx);
            
            // remove siMaxN from each set in s
            for (i = s.size() - 1; i >= 0; i--) {
                si = s.get(i);
                si.removeAll(siMaxN);
                if (si.isEmpty()) {
                    s.remove(i);
                }
            }
        }
        
        return c;
    }

    protected List<TIntSet> copy(List<TIntSet> s) {
        List<TIntSet> c = new ArrayList<TIntSet>(s.size());
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
            c.add(cSet);
        }
        return c;
    }

}
