package algorithms.graphs;

import algorithms.graphs.ApproxGraphSearchZeng.Graph;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ApproxGraphSearchZengTest extends TestCase {

    public ApproxGraphSearchZengTest(String testName) {
        super(testName);
    }
    
    /**
     * get the example graph of similar compounds from
     * Fig 3.b of 
     * Fast processing of graph queries on a large database of small and medium-sized data graphs,
       Dipali Pal, Praveen Rao, Vasil Slavov, Anas Katib,
       Journal of Computer and System Sciences, Volume 82, Issue 6, 2016, Pages 1112-1143,
       https://doi.org/10.1016/j.jcss.2016.04.002.
    
     * @param outputDB list of graphs representing the database
     * @return graph representing the query graph
     */
    public static Graph getG0(List<Graph> outputDB) {
        
        TIntObjectMap<TIntSet> adjMapDB;
        TIntIntMap vLabelsDB;
        TObjectIntMap<PairInt> eLabelsDB;
        int nDB = 15, i;
        adjMapDB = new TIntObjectHashMap<TIntSet>();
        for (i = 0; i < nDB; ++i) {
            adjMapDB.put(i, new TIntHashSet());
        }
        adjMapDB.get(0).add(3);
        adjMapDB.get(1).add(3);
        adjMapDB.get(2).add(3);
        adjMapDB.get(3).add(0);
        adjMapDB.get(3).add(1);
        adjMapDB.get(3).add(2);
        adjMapDB.get(3).add(4);
        adjMapDB.get(4).add(3);
        adjMapDB.get(4).add(5);
        adjMapDB.get(4).add(13);
        adjMapDB.get(5).add(4);
        adjMapDB.get(5).add(6);
        adjMapDB.get(5).add(7);
        adjMapDB.get(6).add(5);
        adjMapDB.get(7).add(5);
        adjMapDB.get(7).add(8);
        adjMapDB.get(7).add(10);
        adjMapDB.get(8).add(7);
        adjMapDB.get(9).add(10);
        adjMapDB.get(10).add(7);
        adjMapDB.get(10).add(9);
        adjMapDB.get(10).add(11);
        adjMapDB.get(11).add(10);
        adjMapDB.get(11).add(12);
        adjMapDB.get(11).add(13);
        adjMapDB.get(12).add(11);
        adjMapDB.get(13).add(4);
        adjMapDB.get(13).add(11);
        adjMapDB.get(13).add(14);
        adjMapDB.get(14).add(13);
        
        vLabelsDB = new TIntIntHashMap();
        vLabelsDB.put(0, (int)'h');
        vLabelsDB.put(1, (int)'h');
        vLabelsDB.put(2, (int)'h');
        vLabelsDB.put(3, (int)'c');
        vLabelsDB.put(4, (int)'c');
        vLabelsDB.put(5, (int)'c');
        vLabelsDB.put(6, (int)'o');
        vLabelsDB.put(7, (int)'n');
        vLabelsDB.put(8, (int)'h');
        vLabelsDB.put(9, (int)'o');
        vLabelsDB.put(10, (int)'c');
        vLabelsDB.put(11, (int)'n');
        vLabelsDB.put(12, (int)'h');
        vLabelsDB.put(13, (int)'c');
        vLabelsDB.put(14, (int)'h');
        
        eLabelsDB = new TObjectIntHashMap<PairInt>();
        eLabelsDB.put(new PairInt(0, 3), (int)'s');
        eLabelsDB.put(new PairInt(1, 3), (int)'s');
        eLabelsDB.put(new PairInt(2, 3), (int)'s');
        eLabelsDB.put(new PairInt(3, 0), (int)'s');
        eLabelsDB.put(new PairInt(3, 1), (int)'s');
        eLabelsDB.put(new PairInt(3, 2), (int)'s');
        eLabelsDB.put(new PairInt(3, 4), (int)'s');
        eLabelsDB.put(new PairInt(4, 3), (int)'s');
        eLabelsDB.put(new PairInt(4, 5), (int)'s');
        eLabelsDB.put(new PairInt(4, 13), (int)'d');
        eLabelsDB.put(new PairInt(5, 4), (int)'s');
        eLabelsDB.put(new PairInt(5, 6), (int)'d');
        eLabelsDB.put(new PairInt(5, 7), (int)'s');
        eLabelsDB.put(new PairInt(6, 5), (int)'d');
        eLabelsDB.put(new PairInt(7, 5), (int)'s');
        eLabelsDB.put(new PairInt(7, 8), (int)'s');
        eLabelsDB.put(new PairInt(7, 10), (int)'s');
        eLabelsDB.put(new PairInt(8, 7), (int)'s');
        eLabelsDB.put(new PairInt(9, 10), (int)'d');
        eLabelsDB.put(new PairInt(10, 9), (int)'d');
        eLabelsDB.put(new PairInt(10, 7), (int)'s');
        eLabelsDB.put(new PairInt(10, 11), (int)'s');
        eLabelsDB.put(new PairInt(11, 10), (int)'s');
        eLabelsDB.put(new PairInt(11, 12), (int)'s');
        eLabelsDB.put(new PairInt(11, 13), (int)'s');
        eLabelsDB.put(new PairInt(12, 11), (int)'s');
        eLabelsDB.put(new PairInt(13, 4), (int)'d');
        eLabelsDB.put(new PairInt(13, 11), (int)'s');
        eLabelsDB.put(new PairInt(13, 14), (int)'s');
        eLabelsDB.put(new PairInt(14, 13), (int)'s');
        
        
        Graph db = new Graph(adjMapDB, vLabelsDB, eLabelsDB);
       
        TIntObjectMap<TIntSet> adjMapQ = copy(adjMapDB);
        adjMapQ.get(13).remove(14);
        adjMapQ.remove(14);
        
        TIntIntMap vLabelsQ = copy(vLabelsDB);
        vLabelsQ.remove(14);
        
        TObjectIntMap<PairInt> eLabelsQ = copy(eLabelsDB);
        eLabelsQ.remove(new PairInt(13, 14));
        eLabelsQ.remove(new PairInt(14, 13));
        
        Graph q = new Graph(adjMapDB, vLabelsDB, eLabelsDB);
        
        outputDB.add(db);
        return q;
    }

    /**
     * Test of approxFullSearch method, of class ApproxGraphSearchZeng.
     */
    public void testApproxFullSearch() throws Exception {
        
    }

    public void testApproxSubSearch() {
        
    }

    public void testApproxSubSearchFilter() {
        
    }

    public void testSuboptimalEditDistance() {
        
    }

    public void testSuboptimalEditDistanceV() {
        
    }

    public void testLowerBoundEditDistance() {
        
    }

    public void testMappingDistance() {
        
    }

    public void testRefinedSuboptimalEditDistance() {
        
    }

    public void testOptimalEditDistance() throws Exception {
        
    }

    public void testReverseAssignment() {
        
    }
    
    private static TIntObjectMap<TIntSet> copy(TIntObjectMap<TIntSet> a) {
        TIntObjectMap<TIntSet> c = new TIntObjectHashMap<TIntSet>();
        
        TIntObjectIterator<TIntSet> iter = a.iterator();
        int k;
        TIntSet v;
        while (iter.hasNext()) {
            iter.advance();
            k = iter.key();
            v = iter.value();
            c.put(k, new TIntHashSet(v));
        }
        
        return c;
    }

    private static TIntIntMap copy(TIntIntMap a) {
        TIntIntMap c = new TIntIntHashMap();
        
        TIntIntIterator iter = a.iterator();
        while (iter.hasNext()) {
            iter.advance();
            c.put(iter.key(), iter.value());
        }
        return c;
    }

    private static TObjectIntMap<PairInt> copy(TObjectIntMap<PairInt> a) {
        TObjectIntMap<PairInt> c = new TObjectIntHashMap<PairInt>();
        
        TObjectIntIterator<PairInt> iter = a.iterator();
        PairInt p;
        int v;
        while (iter.hasNext()) {
            iter.advance();
            p = iter.key();
            v = iter.value();
            c.put(p.copy(), v);
        }
        
        return c;
    }
    
}
