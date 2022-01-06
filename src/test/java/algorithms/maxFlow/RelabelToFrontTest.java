package algorithms.maxFlow;

import algorithms.util.PairInt;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class RelabelToFrontTest extends TestCase {
    
    public RelabelToFrontTest(String testName) {
        super(testName);
    }
    
    public void testDischarge() {
        
        // Figure 26.10 from Cormen et al. Introduction to Algorithms"
        
        TIntObjectMap<TIntSet> adj = new TIntObjectHashMap<TIntSet>();
        int nV = 5;
        int i;
        for (i = 0; i < nV; ++i) {
            adj.put(i, new TIntHashSet());
        }
        adj.get(0).add(1);
        adj.get(0).add(2);
        adj.get(1).add(2);
        adj.get(1).add(4);
        adj.get(2).add(3);
        adj.get(3).add(1);
        adj.get(3).add(4);
        int srcIdx = 0;
        int sinkIdx = 4;
        
        TObjectDoubleMap<PairInt> cap = new TObjectDoubleHashMap<PairInt>();
        cap.put(new PairInt(0, 1), 12);
        cap.put(new PairInt(0, 2), 14);
        cap.put(new PairInt(1, 2), 5);
        cap.put(new PairInt(1, 4), 16);
        cap.put(new PairInt(2, 3), 8);
        cap.put(new PairInt(3, 1), 7);
        cap.put(new PairInt(3, 4), 10);
        
        //RelabelToFront rtf = new RelabelToFront(adj, cap, srcIdx, sinkIdx);
        
        //add flow to rtf.TObjectDoubleMap<PairInt> to mimic test state
    }
}
