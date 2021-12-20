package algorithms.graphs;

import algorithms.graphs.ApproxGraphSearchZeng.Graph;
import algorithms.graphs.ApproxGraphSearchZeng.Norm;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import static junit.framework.Assert.assertEquals;
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
        
        Graph q = new Graph(adjMapQ, vLabelsQ, eLabelsQ);
        
        outputDB.add(db);
        return q;
    }

    public void testReverseAssignment() {
        int[] a = new int[]{3,4,5};
        ApproxGraphSearchZeng ags = new ApproxGraphSearchZeng();
        TIntIntMap r = ags.reverseAssignment(a);
        TIntIntMap expected = new TIntIntHashMap();
        expected.put(3, 0);
        expected.put(4, 1);
        expected.put(5, 2);
        assertEquals(expected.size(), r.size());
        TIntIntIterator iter = r.iterator();
        int k, v;
        while (iter.hasNext()) {
            iter.advance();
            k = iter.key();
            v = iter.value();
            assertEquals(expected.get(k), v);
            expected.remove(k);
        }
        assertTrue(expected.isEmpty());
    }
    
    public void testApproxFullSearch() throws Exception {
        
        //L_M <= lambda <= rho <= tau
        
        
    }

    public void testMappingDistance() {
        List<Graph> dbs = new ArrayList<Graph>();
        Graph q = ApproxGraphSearchZengTest.getG0(dbs);
        int nV = 14; // dB has 15, q has 14
        
        StarStructure[] stars = StarStructure.createStarStructureMultiset(q);
        assertEquals(nV, stars.length);
        
        StarStructure[] starDB = StarStructure.createStarStructureMultiset(dbs.get(0));
        assertEquals(nV+1, starDB.length);
        
        int mDist, mDist2;
        int[] expectedAssignments, assignments;
        int i;
        assignments = new int[nV];
        for (i = 0; i < assignments.length; ++i) {
            assignments[i] = i;
        }

        ApproxGraphSearchZeng ags = new ApproxGraphSearchZeng();
        
        ags.setEdgesAreLabeled(false);
        //4·k1 + 1·k2 + (max{degree(g1),degree(g2)}+1)·k3
        //for edge ins/del cost = 1 instead, we have 2*k1 instead of 4*k1
        //2·k1 + 1·k2 + 2*(max{degree(g1),degree(g2)}+1)·k3; expected total=2*0 + 1*1 + 0*(max degree) = 1
        //k1 edge insertion/deletion operations, 
        //k2 vertex insertion/deletion operations 
        //k3 vertex relabellings
        //i=13 sum=1  :  k1=0,k2=1,k3=0
        mDist = ags.mappingDistance(stars, starDB, assignments);
        /*System.out.printf("not normalized: mapping dist=%d\n", mDist);
        System.out.printf("label(q[13])=%d, label(db[13])=%d\n", 
            stars[13].rootLabel, starDB[13].rootLabel);
        System.out.printf("deg(q[13])=%d, deg(db[13])=%d\n", 
            stars[13].vLabels.length, starDB[13].vLabels.length);
        System.out.printf("label(q[14])=__, label(db[14])=%d\n", 
            starDB[14].rootLabel);
        System.out.printf("deg(q[14])=__, deg(db[14])=%d\n", 
            starDB[14].vLabels.length);
        System.out.printf("nQ=%d  nDB=%d\n", stars.length, starDB.length);*/
        assertEquals(1, mDist);
        
        ags.setEdgesAreLabeled(true);
        //4·k1 + 1·k2 + 1*(max{degree(g1),degree(g2)}+1)·k3
        //for edge ins/del cost = 1 instead, we have 2*k1 instead of 4*k1
        //2·k1 + 1·k2 + 1*(max{degree(g1),degree(g2)}+1)·k3; expected total=2*1 + 0*1 + 0*1(max degree) = 2
        //k1 edge insertion/deletion operations, 
        //k2 vertex insertion/deletion operations 
        //k3 vertex relabellings
        //i=13 ed=2 
        mDist = ags.mappingDistance(stars, starDB, assignments);
        /*System.out.printf("not normalized: edges, mapping dist=%d\n", mDist);
        System.out.printf("label(q[13])=%d, label(db[13])=%d\n", 
            stars[13].rootLabel, starDB[13].rootLabel);
        System.out.printf("deg(q[13])=%d, deg(db[13])=%d\n", 
            stars[13].vLabels.length, starDB[13].vLabels.length);
        System.out.printf("label(q[14])=__, label(db[14])=%d\n", 
            starDB[14].rootLabel);
        System.out.printf("deg(q[14])=__, deg(db[14])=%d\n", 
            starDB[14].vLabels.length);
        System.out.printf("nQ=%d  nDB=%d\n", stars.length, starDB.length);*/
        assertEquals(2, mDist);
        
        Norm norm = ApproxGraphSearchZeng.normalize(stars, starDB);
        stars = norm.sg1;
        starDB = norm.sg2;
        boolean swapped = norm.swapped;
        
        expectedAssignments = new int[nV+1];
        // graphs are identical except for vertex=13
        for (i = 0; i < expectedAssignments.length; ++i) {
            expectedAssignments[i] = i;
        }
        
        double[][] distM;
        //distM = StarStructure.createDistanceMatrixV(stars, starDB);
        distM = StarStructure.createDistanceMatrix(stars, starDB);
        assignments = ApproxGraphSearchZeng.balancedBipartiteAssignment(distM);
        assertTrue(Arrays.equals(expectedAssignments, assignments));
        
        
        ags.setEdgesAreLabeled(false);
        //4·k1 + 1·k2 + (max{degree(g1),degree(g2)}+1)·k3
        //for edge ins/del cost = 1 instead, we have 2*k1 instead of 4*k1
        //2·k1 + 1·k2 + 1*(max{degree(g1),degree(g2)}+1)·k3;  expected total=2*0 + 1*0 + 1*1(max degree=3)=3
        //k1 edge insertion/deletion operations, 
        //k2 vertex insertion/deletion operations 
        //k3 vertex relabellings
        //i=13 ed=1  : nVMismatched=1  =>  1*nVMismatched = 1
        //i=14 ed=2  : nVMismatched=1 + 1 vertex relabelling = 2
        mDist = ags.mappingDistance(stars, starDB, assignments);
        /*System.out.printf("normalized: mapping dist=%d\n", mDist);
        System.out.printf("swapped=%b, mapping dist=%d\n", swapped, mDist);
        System.out.printf("swapped=%b, mapping dist=%d\n", swapped, mDist);
        System.out.printf("label(q[13])=%d, label(db[13])=%d\n", 
            stars[13].rootLabel, starDB[13].rootLabel);
        System.out.printf("deg(q[13])=%d, deg(db[13])=%d\n", 
            stars[13].vLabels.length, starDB[13].vLabels.length);
        System.out.printf("label(q[14])=%d, label(db[14])=%d\n", 
            stars[14].rootLabel, starDB[14].rootLabel);
        System.out.printf("deg(q[14])=%d, deg(db[14])=%d\n", 
            stars[14].vLabels.length, starDB[14].vLabels.length);*/
        assertEquals(3, mDist);
        
        ags.setEdgesAreLabeled(true);
        //for edge ins/del cost = 1 instead, we have 2*k1 instead of 4*k1
        //2·k1 + 1·k2 + 1*(max{degree(g1),degree(g2)}+1)·k3; expected total= 2*1 + 1*0 + 1*1*(max degree=3) = 5
        //i=13 ed=2 
        //i=14 ed=3, sum=5
        //4·k1 + 1·k2 + 2.*(max{degree(g1), degree(g2)}+1)·k3
        //k1 edge insertion/deletion operations, 
        //k2 vertex insertion/deletion operations 
        //k3 vertex relabellings
        mDist = ags.mappingDistance(stars, starDB, assignments);
        /*System.out.printf("swapped=%b,has  edges, mapping dist=%d\n", swapped, mDist);
        System.out.printf("swapped=%b, has edges, mapping dist=%d\n", swapped, mDist);
        System.out.printf("label(q[13])=%d, label(db[13])=%d\n", 
            stars[13].rootLabel, starDB[13].rootLabel);
        System.out.printf("deg(q[13])=%d, deg(db[13])=%d\n", 
            stars[13].vLabels.length, starDB[13].vLabels.length);
        System.out.printf("label(q[14])=%d, label(db[14])=%d\n", 
            stars[14].rootLabel, starDB[14].rootLabel);
        System.out.printf("deg(q[14])=%d, deg(db[14])=%d\n", 
            stars[14].vLabels.length, starDB[14].vLabels.length);*/
        assertEquals(5, mDist);
    }
    
    public void estLowerBoundEditDistance() {
        List<Graph> dbs = new ArrayList<Graph>();
        Graph q = ApproxGraphSearchZengTest.getG0(dbs);
        int nV = 14; // dB has 15, q has 14
        
        StarStructure[] stars = StarStructure.createStarStructureMultiset(q);
        assertEquals(nV, stars.length);
        
        ApproxGraphSearchZeng ags = new ApproxGraphSearchZeng();
        double lM = ags.lowerBoundEditDistance(stars, stars, 0);
        
    }

    public void testRefinedSuboptimalEditDistance() {
        
    }

    public void testOptimalEditDistance() throws Exception {
        
    }

    public void testApproxSubSearch() {
        
    }

    public void testApproxSubSearchFilter() {
        
    }

    public void testSuboptimalEditDistance() {
        
    }

    public void testSuboptimalEditDistanceV() {
        
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
