package algorithms.graphs;

import algorithms.graphs.ApproxGraphSearchZeng.Graph;
import algorithms.graphs.ApproxGraphSearchZeng.Norm;
import algorithms.graphs.ApproxGraphSearchZeng.Result;
import algorithms.matrix.MatrixUtil;
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
                
        List<Graph> dbs = new ArrayList<Graph>();
        Graph q = ApproxGraphSearchZengTest.getG1(dbs);
        
        boolean useAsFilterWithoutOptimal = false;
        ApproxGraphSearchZeng ags = new ApproxGraphSearchZeng();
        ags.setEdgesAreLabeled(true);
        
        
        double w = 1.5;
        
        List<Result> results = ags.approxFullSearch(q, dbs, w, 
            useAsFilterWithoutOptimal);
        
        assertEquals(1, results.size());
        
        assertEquals(2, results.get(0).dbGraphIndex);
        
        int[] assignExpected = new int[]{4, 2, 1, 0, 3};
        
        int[] assignResult = results.get(0).assignment;
        
        System.out.printf("result assign=%s\n", Arrays.toString(assignResult));
        
        assertTrue(Arrays.equals(assignExpected, assignResult));
        
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
    
    public void testLowerBoundEditDistance() {
        
        List<Graph> dbs = new ArrayList<Graph>();
        Graph q = ApproxGraphSearchZengTest.getG0(dbs);
        int nV = 14; // dB has 15, q has 14
        
        StarStructure[] stars = StarStructure.createStarStructureMultiset(q);
        assertEquals(nV, stars.length);
        
        ApproxGraphSearchZeng ags = new ApproxGraphSearchZeng();
        
        StarStructure[] starDB = StarStructure.createStarStructureMultiset(dbs.get(0));
        assertEquals(nV+1, starDB.length);
        
        int[] expectedAssignments, assignments;
        int i, mDist;
        double lM;
        
        assignments = new int[stars.length];
        // graphs are identical except for vertex=13
        for (i = 0; i < assignments.length; ++i) {
            assignments[i] = i;
        }
        mDist = ags.mappingDistance(stars, stars, assignments);
        lM = ags.lowerBoundEditDistance(stars, stars, mDist);
        
        System.out.printf("same graphs: L_M=%.4f\n", lM);
        
        
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
        
        mDist = ags.mappingDistance(stars, starDB, assignments);
        lM = ags.lowerBoundEditDistance(stars, starDB, mDist);
        System.out.printf("V labels only: L_M=%.4f\n", lM);
        
        ags.setEdgesAreLabeled(true);
        mDist = ags.mappingDistance(stars, starDB, assignments);
        lM = ags.lowerBoundEditDistance(stars, starDB, mDist);
        System.out.printf("V and E labeled: L_M=%.4f\n", lM);
    }
    
    public void testSuboptimalCP() {
        //Table 1 of Justice & Hero
        
        int eps = Integer.MIN_VALUE;
        int[] vLabelsA = new int[]{(int)'a', (int)'b', (int)'b', eps, eps};
        int[] vLabelsB = new int[]{(int)'b', eps, (int)'b', (int)'a', eps};
        int[] assign = new int[]{4-1, 1-1, 3-1, 5-1, 2-1};// subtractng 1 to make 0-based indexes
        int[][] P = MatrixUtil.createPermutationMatrix(assign);
        
        int i, j;
        int[][] c = new int[vLabelsA.length][];
        for (i = 0; i < P.length; ++i) {
            c[i] = new int[P[i].length];
            //NLK, changing the fill values from 0 to 1 because the non-matches will cost 1
            Arrays.fill(c[i], 1);
        }
        
        for (i = 0; i < assign.length; ++i) {
            j = assign[i];
            if (vLabelsA[i] == vLabelsB[j]) {
                //NLK, changing to not penalize if vertex labels match
                //c[i][j] = 1;
                c[i][j] = 0;
            }
        }
        
        // this is the optimal permutation and has equal mapped vertex attributes.
        //   if term1 is truly an edit cost it should equal 0.
        double term1 = 0;
        for (i = 0; i < c.length; ++i) {
            for (j = 0; j < c[i].length; ++j) {
                term1 += c[i][j] * P[i][j];
            }
        }
        
        System.out.printf("term1=%.2f\n", term1);
        
    }
    
    public void testSuboptimalAndRefinedEditDistance() throws InterruptedException {
        // tau
        
        List<Graph> dbs = new ArrayList<Graph>();
        Graph q = ApproxGraphSearchZengTest.getG0(dbs);
        int nV = 14; // dB has 15, q has 14
        
        int[] expectedAssignments, assignments, refinedAssign;
        int i, mDist;
        double tau, rho, lambda;
        boolean useEdges = false;
        double[][] distM;
        int[][] a1, a2;
        
        StarStructure[] stars = StarStructure.createStarStructureMultiset(q);
        assertEquals(nV, stars.length);
                
        StarStructure[] starDB = StarStructure.createStarStructureMultiset(dbs.get(0));
        assertEquals(nV+1, starDB.length);
           
        Norm norm = ApproxGraphSearchZeng.normalize(stars, starDB);
        stars = norm.sg1;
        starDB = norm.sg2;
        boolean swapped = norm.swapped;
        
        a1 = ApproxGraphSearchZeng.createAdjacencyMatrix(stars);
        a2 = ApproxGraphSearchZeng.createAdjacencyMatrix(starDB);
                        
        ApproxGraphSearchZeng ags = new ApproxGraphSearchZeng();
        ags.setEdgesAreLabeled(useEdges);
        distM = StarStructure.createDistanceMatrix(stars, stars);
        assignments = ApproxGraphSearchZeng.balancedBipartiteAssignment(distM);
        tau = ags.suboptimalEditDistanceV(stars, stars, a1, a1, assignments);
        System.out.printf("normalized, same graphs, suboptimal: tau=%.4f\n", tau);
        
        refinedAssign = Arrays.copyOf(assignments, assignments.length);
        rho = ags.refinedSuboptimalEditDistance(stars, stars, q.eLabels, q.eLabels, a1, a1, refinedAssign, tau, distM);
        System.out.printf("normalized, same graphs, refined suboptimal: rho=%.4f\n", rho);
        
        //lambda = ags.optimalEditDistance(stars, stars, e1, e1, a1, a1, refinedAssign, tau);
        //System.out.printf("normalized, same graphs, optimal: lambda=%.4f\n", lambda);
        
        
        ags.setEdgesAreLabeled(useEdges);
        distM = StarStructure.createDistanceMatrix(stars, starDB);
        assignments = ApproxGraphSearchZeng.balancedBipartiteAssignment(distM);
        tau = ags.suboptimalEditDistanceV(stars, starDB, a1, a2, assignments);
        System.out.printf("normalized, suboptimal: tau=%.4f\n", tau);
        
        refinedAssign = Arrays.copyOf(assignments, assignments.length);
        if (swapped) {
            rho = ags.refinedSuboptimalEditDistance(stars, starDB, dbs.get(0).eLabels, q.eLabels, a1, a2, refinedAssign, tau, distM);
        } else {
            rho = ags.refinedSuboptimalEditDistance(stars, starDB, q.eLabels, dbs.get(0).eLabels, a1, a2, refinedAssign, tau, distM);
        }
        System.out.printf("normalized, refined suboptimal: rho=%.4f\n", rho);

        //lambda = ags.optimalEditDistance(stars, starDB, e1, e2, a1, a2, refinedAssign, tau);
        //System.out.printf("normalized, optimal: lambda=%.4f\n", lambda);
        
        
        ags.setEdgesAreLabeled(true);
        distM = StarStructure.createDistanceMatrix(stars, starDB);
        assignments = ApproxGraphSearchZeng.balancedBipartiteAssignment(distM);
        if (swapped) {
            tau = ags.suboptimalEditDistance(stars, starDB, dbs.get(0).eLabels, q.eLabels, assignments);
        } else {
            tau = ags.suboptimalEditDistance(stars, starDB, q.eLabels, dbs.get(0).eLabels, assignments);
        }
        
        System.out.printf("normalized, edges, suboptimal: tau=%.4f\n", tau);
        
        refinedAssign = Arrays.copyOf(assignments, assignments.length);
        if (swapped) {
            rho = ags.refinedSuboptimalEditDistance(stars, starDB, dbs.get(0).eLabels, q.eLabels, a1, a2, refinedAssign, tau, distM);
        } else {
            rho = ags.refinedSuboptimalEditDistance(stars, starDB, q.eLabels, dbs.get(0).eLabels, a1, a2, refinedAssign, tau, distM);
        }
        
        System.out.printf("normalized, edges, refined suboptimal: rho=%.4f\n", rho);

        //lambda = ags.optimalEditDistance(stars, starDB, e1, e2, a1, a2, refinedAssign, tau);
        //System.out.printf("normalized, edges, optimal: lambda=%.4f\n", lambda);
        
        System.out.println("expecting L_M <= lambda <= rho <= tau");
        //L_M <= lambda <= rho <= tau
        
    }

    public void testOptimalEditDistance() throws Exception {
        System.out.println("\ntestOptimalEditDistance");
        //example graphs from Figure 1 of
        //A Coding Method for Efficient Subgraph Querying on Vertex- and Edge-Labeled Graphs
        //Zhu et al, May 2014 PLoS ONE 9(5):e97178
        //DOI:10.1371/journal.pone.0097178
        
        // 4!=24, 5!=120, so using a small test graph
        List<Graph> dbs = new ArrayList<Graph>();
        Graph q = ApproxGraphSearchZengTest.getG1(dbs);
        int nQ = 5;
        
        // inserting query at beginning to make sure edit costs are 0
        dbs.add(0, q);
        
        //int nD1 = 5;
        //int nD2 = 4;
        //int nD3 = 5;
        //int nD4 = 5;
        
        ApproxGraphSearchZeng ags = new ApproxGraphSearchZeng();
        
        // q and db.get(3) should be closest.
        
        int i, mappingDist;
        double lM, tau, rho, lambda;
        
        List<Result> results;
        double[][] distM;
        boolean swapped;
        int[][] a1, a2;
        int[] assignments, refinedAssign;
        StarStructure[] stars, starsDB;
        Graph g;
        System.out.println("expecting L_M <= lambda <= rho <= tau");
        //L_M <= lambda <= rho <= tau
        for (boolean useEdges : new boolean[]{false, true}) {
            for (i = 0; i < dbs.size(); ++i) {
                g = dbs.get(i);

                stars = StarStructure.createStarStructureMultiset(q);
                assertEquals(nQ, stars.length);

                starsDB = StarStructure.createStarStructureMultiset(g);
                assertEquals(g.vLabels.size(), starsDB.length);

                Norm norm = ApproxGraphSearchZeng.normalize(stars, starsDB);
                stars = norm.sg1;
                starsDB = norm.sg2;
                swapped = norm.swapped;

                a1 = ApproxGraphSearchZeng.createAdjacencyMatrix(stars);
                a2 = ApproxGraphSearchZeng.createAdjacencyMatrix(starsDB);

                ags.setEdgesAreLabeled(useEdges);
                distM = StarStructure.createDistanceMatrix(stars, starsDB);
                assignments = ApproxGraphSearchZeng.balancedBipartiteAssignment(distM);

                mappingDist = ags.mappingDistance(stars, starsDB, assignments);

                lM = ags.lowerBoundEditDistance(stars, starsDB, mappingDist);

                tau = ags.suboptimalEditDistanceV(stars, starsDB, a1, a2, assignments);

                refinedAssign = Arrays.copyOf(assignments, assignments.length);
                if (swapped) {
                    rho = ags.refinedSuboptimalEditDistance(stars, starsDB, dbs.get(0).eLabels, q.eLabels, a1, a2, refinedAssign, tau, distM);
                    lambda = ags.optimalEditDistance(stars, starsDB, dbs.get(0).eLabels, q.eLabels, a1, a2, refinedAssign, tau);  
                } else {
                    rho = ags.refinedSuboptimalEditDistance(stars, starsDB, q.eLabels, dbs.get(0).eLabels, a1, a2, refinedAssign, tau, distM);
                    lambda = ags.optimalEditDistance(stars, starsDB, q.eLabels, dbs.get(0).eLabels, a1, a2, refinedAssign, tau);  
                }

                System.out.printf("normalized, edges=%b, i=%d:\n   lM(lower)=%.2f, "
                    + "lambda(opt)=%.2f, rho(refSubOpt)=%.2f, tau(subOpt)=%.2f\n",
                    useEdges, i, lM, lambda, rho, tau);  
            }
        }
        System.out.println("expecting L_M <= lambda <= rho <= tau");
        //L_M <= lambda <= rho <= tau
    }

    public void testApproxSubSearch() throws Exception {
        
        System.out.println("testApproxSubSearch()");
                
        List<Graph> dbs = new ArrayList<Graph>();
        Graph q = ApproxGraphSearchZengTest.getG1(dbs);
        
        // adding a copy of D3 Figure 1 but with 2 extra vertexes and labels.
        // 7! = 5040
        Graph d3c = Graph.copy(dbs.get(2));
        d3c.adjMap.get(4).add(5);
        d3c.adjMap.put(5, new TIntHashSet());
        d3c.adjMap.get(5).add(4);
        d3c.adjMap.get(5).add(6);
        d3c.adjMap.put(6, new TIntHashSet());
        d3c.adjMap.get(6).add(5);
        d3c.vLabels.put(5, (int)'E');
        d3c.vLabels.put(6, (int)'F');
        d3c.eLabels.put(new PairInt(4, 5), 'd');
        d3c.eLabels.put(new PairInt(5, 6), 'e');
        dbs.add(d3c);
        
        // adding q to the database graph also to check that a perfect match is also found and has edit dist = 0.
        dbs.add(0, q);
        
        boolean useAsFilterWithoutOptimal = true;
        ApproxGraphSearchZeng ags = new ApproxGraphSearchZeng();
        ags.setEdgesAreLabeled(true);
        
        double w = 0.5;
        
        List<Result> results = ags.approxSubSearch(q, dbs, w, 
            useAsFilterWithoutOptimal);
        
        // expecting to have kept 0, 3, 5
        
        assertEquals(3, results.size());
        
        assertEquals(0, results.get(0).dbGraphIndex);
        assertEquals(3, results.get(1).dbGraphIndex);
        assertEquals(5, results.get(2).dbGraphIndex);
        
        int[] assignExpected0 = new int[]{0, 1, 2, 3, 4};
        int[] assignExpected3 = new int[]{4, 2, 1, 0, 3};
        int[] assignExpected5 = new int[]{4, 2, 1, 0, 3, 5, 6};
        //                                3, 2, 1, 4, 0, 5, 6
        
        for (int i = 0; i < results.size(); ++i) {
            System.out.printf("result(%d=dbi%d) assign=%s\n",
                i, results.get(i).dbGraphIndex, Arrays.toString(results.get(i).assignment));
        }
        // assertTrue(Arrays.equals(assignExpected, assignResult));
        
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
    
    /**
     * get the example graph of similar compounds from
     * Fig 3.b of 
     * Fast processing of graph queries on a large database of small and medium-sized data graphs,
       Dipali Pal, Praveen Rao, Vasil Slavov, Anas Katib,
       Journal of Computer and System Sciences, Volume 82, Issue 6, 2016, Pages 1112-1143,
       https://doi.org/10.1016/j.jcss.2016.04.002.
    
     @param outputDB list of graphs representing the database
     @return graph representing the query graph
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
    
    /**
     * get the example graphs from Figure 1 of
     * A Coding Method for Efficient Subgraph Querying on Vertex- and Edge-Labeled Graphs
       Zhu et al, May 2014 PLoS ONE 9(5):e97178
       DOI:10.1371/journal.pone.0097178
    
     @param outputDB list of graphs representing the database
     @return graph representing the query graph
     */
    public static Graph getG1(List<Graph> outputDB) {
        
        int nQ = 5;
        int nD1 = 5;
        int nD2 = 4;
        int nD3 = 5;
        int nD4 = 5;
        int i;
        
        TIntObjectMap<TIntSet> adjMapQ;
        TIntIntMap vLabelsQ;
        TObjectIntMap<PairInt> eLabelsQ;
        adjMapQ = new TIntObjectHashMap<TIntSet>();
        for (i = 0; i < nQ; ++i) {
            adjMapQ.put(i, new TIntHashSet());
        }
        adjMapQ.get(0).add(1);
        adjMapQ.get(1).add(0);
        adjMapQ.get(1).add(2);
        adjMapQ.get(2).add(1);
        adjMapQ.get(2).add(3);
        adjMapQ.get(2).add(4);
        adjMapQ.get(3).add(2);
        adjMapQ.get(4).add(2);
        
        vLabelsQ = new TIntIntHashMap();
        vLabelsQ.put(0, (int)'D');
        vLabelsQ.put(1, (int)'B');
        vLabelsQ.put(2, (int)'A');
        vLabelsQ.put(3, (int)'C');
        vLabelsQ.put(4, (int)'C');
        
        eLabelsQ = new TObjectIntHashMap<PairInt>();
        eLabelsQ.put(new PairInt(0, 1), (int)'b');
        eLabelsQ.put(new PairInt(1, 2), (int)'a');
        eLabelsQ.put(new PairInt(2, 3), (int)'a');
        eLabelsQ.put(new PairInt(2, 4), (int)'c');
        
        Graph q = new Graph(adjMapQ, vLabelsQ, eLabelsQ);
       
        // ---D4--------
        TIntObjectMap<TIntSet> adjMapD4 = new TIntObjectHashMap<TIntSet>();
        TIntIntMap vLabelsD4;
        TObjectIntMap<PairInt> eLabelsD4;
        for (i = 0; i < nD4; ++i) {
            adjMapD4.put(i, new TIntHashSet());
        }
        adjMapD4.get(0).add(1);
        adjMapD4.get(0).add(2);
        adjMapD4.get(0).add(3);
        adjMapD4.get(0).add(4);
        adjMapD4.get(1).add(0);
        adjMapD4.get(2).add(0);
        adjMapD4.get(3).add(0);
        adjMapD4.get(4).add(0);
        
        vLabelsD4 = new TIntIntHashMap();
        vLabelsD4.put(0, (int)'B');
        vLabelsD4.put(1, (int)'A');
        vLabelsD4.put(2, (int)'C');
        vLabelsD4.put(3, (int)'D');
        vLabelsD4.put(4, (int)'C');
        
        eLabelsD4 = new TObjectIntHashMap<PairInt>();
        eLabelsD4.put(new PairInt(0, 1), (int)'a');
        eLabelsD4.put(new PairInt(0, 2), (int)'c');
        eLabelsD4.put(new PairInt(0, 3), (int)'b');
        eLabelsD4.put(new PairInt(0, 4), (int)'a');
        
        Graph d4 = new Graph(adjMapD4, vLabelsD4, eLabelsD4);
        
        // ---D3--------
        TIntObjectMap<TIntSet> adjMapD3 = new TIntObjectHashMap<TIntSet>();
        TIntIntMap vLabelsD3;
        TObjectIntMap<PairInt> eLabelsD3;
        for (i = 0; i < nD3; ++i) {
            adjMapD3.put(i, new TIntHashSet());
        }
        adjMapD3.get(0).add(1);
        adjMapD3.get(1).add(0);
        adjMapD3.get(1).add(2);
        adjMapD3.get(1).add(3);
        adjMapD3.get(2).add(1);
        adjMapD3.get(2).add(3);
        adjMapD3.get(2).add(4);
        adjMapD3.get(3).add(1);
        adjMapD3.get(3).add(2);
        adjMapD3.get(4).add(2);
        
        vLabelsD3 = new TIntIntHashMap();
        vLabelsD3.put(0, (int)'C');
        vLabelsD3.put(1, (int)'A');
        vLabelsD3.put(2, (int)'B');
        vLabelsD3.put(3, (int)'C');
        vLabelsD3.put(4, (int)'D');
        
        eLabelsD3 = new TObjectIntHashMap<PairInt>();
        eLabelsD3.put(new PairInt(0, 1), (int)'a');
        eLabelsD3.put(new PairInt(1, 2), (int)'a');
        eLabelsD3.put(new PairInt(1, 3), (int)'c');
        eLabelsD3.put(new PairInt(2, 3), (int)'c');
        eLabelsD3.put(new PairInt(2, 4), (int)'b');
        
        Graph d3 = new Graph(adjMapD3, vLabelsD3, eLabelsD3);
        
        // ---D2--------
        TIntObjectMap<TIntSet> adjMapD2 = new TIntObjectHashMap<TIntSet>();
        TIntIntMap vLabelsD2;
        TObjectIntMap<PairInt> eLabelsD2;
        for (i = 0; i < nD2; ++i) {
            adjMapD2.put(i, new TIntHashSet());
        }
        adjMapD2.get(0).add(1);
        adjMapD2.get(0).add(2);
        adjMapD2.get(0).add(3);
        adjMapD2.get(1).add(0);
        adjMapD2.get(1).add(3);
        adjMapD2.get(2).add(0);
        adjMapD2.get(2).add(3);
        adjMapD2.get(3).add(0);
        adjMapD2.get(3).add(1);
        adjMapD2.get(3).add(2);
        
        vLabelsD2 = new TIntIntHashMap();
        vLabelsD2.put(0, (int)'C');
        vLabelsD2.put(1, (int)'D');
        vLabelsD2.put(2, (int)'A');
        vLabelsD2.put(3, (int)'B');
        
        eLabelsD2 = new TObjectIntHashMap<PairInt>();
        eLabelsD2.put(new PairInt(0, 1), (int)'b');
        eLabelsD2.put(new PairInt(0, 2), (int)'a');
        eLabelsD2.put(new PairInt(0, 3), (int)'a');
        eLabelsD2.put(new PairInt(1, 3), (int)'b');
        eLabelsD2.put(new PairInt(2, 3), (int)'c');
        
        Graph d2 = new Graph(adjMapD2, vLabelsD2, eLabelsD2);
        
        // ---D1--------
        TIntObjectMap<TIntSet> adjMapD1 = new TIntObjectHashMap<TIntSet>();
        TIntIntMap vLabelsD1;
        TObjectIntMap<PairInt> eLabelsD1;
        for (i = 0; i < nD1; ++i) {
            adjMapD1.put(i, new TIntHashSet());
        }
        adjMapD1.get(0).add(1);
        adjMapD1.get(0).add(2);
        adjMapD1.get(1).add(0);
        adjMapD1.get(1).add(2);
        adjMapD1.get(1).add(3);
        adjMapD1.get(2).add(0);
        adjMapD1.get(2).add(1);
        adjMapD1.get(2).add(4);
        adjMapD1.get(3).add(1);
        adjMapD1.get(4).add(2);
        
        vLabelsD1 = new TIntIntHashMap();
        vLabelsD1.put(0, (int)'A');
        vLabelsD1.put(1, (int)'B');
        vLabelsD1.put(2, (int)'C');
        vLabelsD1.put(3, (int)'C');
        vLabelsD1.put(4, (int)'D');
        
        eLabelsD1 = new TObjectIntHashMap<PairInt>();
        eLabelsD1.put(new PairInt(0, 1), (int)'a');
        eLabelsD1.put(new PairInt(0, 2), (int)'a');
        eLabelsD1.put(new PairInt(1, 2), (int)'b');
        eLabelsD1.put(new PairInt(1, 3), (int)'b');
        eLabelsD1.put(new PairInt(2, 4), (int)'b');
        
        Graph d1 = new Graph(adjMapD1, vLabelsD1, eLabelsD1);
        
        outputDB.add(d1);
        outputDB.add(d2);
        outputDB.add(d3);
        outputDB.add(d4);
        
        return q;
    }
}
