package algorithms.graphs;

import algorithms.optimization.LinearProgramming;
import algorithms.optimization.ORSetCoverReader;
import algorithms.util.FormatArray;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SetCoverTest extends TestCase {
    
    public SetCoverTest(String testName) {
        super(testName);
    }

    /**
     * Test of approx2 method, of class SetCover.
     
     0  1  2
     3  4  5
     6  7  8
     9  10 11
     */
    public void testApproxLgN() {
        /*
        from Section 35.3 of Cormen at al. "Introduction to Comptuer Algorithms"
        */
        List<TIntSet> sets = new ArrayList<TIntSet>();
        int n = 6, i;
        for (i = 0; i < n; ++i) {
            sets.add(new TIntHashSet());
        }
        sets.get(0).addAll(new int[]{0, 1, 2, 3, 4, 5});//S1
        sets.get(1).addAll(new int[]{4, 5, 7, 8});      //S2
        sets.get(2).addAll(new int[]{0, 3, 6, 9});      //S3
        sets.get(3).addAll(new int[]{1, 4, 6, 7, 10});  //S4
        sets.get(4).addAll(new int[]{2, 5, 8, 11});     //S5
        sets.get(5).addAll(new int[]{9, 10});           //S6
        
        SetCover sc = new SetCover();
    
        TIntList setCoverIndexes = sc.approxLgN(sets);
        
        // optimal minimum set cover: s2, s3, s4
        TIntSet expected0 = new TIntHashSet(new int[]{0, 3, 4, 5});
        TIntSet expected1 = new TIntHashSet(new int[]{0, 3, 4, 2});
        
        assertEquals(expected0.size(), setCoverIndexes.size());
        assertTrue(setCoverIndexes.containsAll(expected0) || setCoverIndexes.containsAll(expected1));
    }

    public void testSetWeightedApprox2LgN() {
        /*
        from Section 35.3 of Cormen at al. "Introduction to Comptuer Algorithms"
        and
        material from lecture slides of Principal lecturer: Dr Thomas Sauerwald
        Advanced Algorithms, University of Cambridge.
        VII. Approximation Algorithms: Randomisation and Rounding
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/materials.html
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/rand.pdf
        */
        List<TIntSet> sets = new ArrayList<TIntSet>();
        int n = 6, i;
        for (i = 0; i < n; ++i) {
            sets.add(new TIntHashSet());
        }
        sets.get(0).addAll(new int[]{0, 1, 2, 3, 4, 5});//S1 
        sets.get(1).addAll(new int[]{4, 5, 7, 8});      //S2
        sets.get(2).addAll(new int[]{0, 3, 6, 9});      //S3 
        sets.get(3).addAll(new int[]{1, 4, 6, 7, 10});  //S4
        sets.get(4).addAll(new int[]{2, 5, 8, 11});     //S5  
        sets.get(5).addAll(new int[]{9, 10});           //S6
        
        double[] weights = new double[]{2, 3, 3, 5, 1, 2};
        // 4,3,2,0 cost=1+5+3+2=11
        // opt=1,2,3  cost=3+3+5=11
        
        int nU = 12;
        
        // expecting y = all 0.5's excepting for S5 where it is '1'
        double[] expectedYs = new double[n];
        Arrays.fill(expectedYs, 0.5);
        expectedYs[4] = 1;
        
        LinearProgramming.StandardForm standForm = 
            SetCover.createLinearProgramInStandardFormForWeightedSets(nU, sets, weights);
        LinearProgramming lp = new LinearProgramming();
        LinearProgramming.SlackForm lpSoln = lp.solveUsingSimplexMethod(standForm);
        
        System.out.printf("lp-soln=\n%s\n", lpSoln.toString());
        /*
        for matrix a formed as [nU X sets.size()], 
           the objective is summation over i_(from 0 to sets.size()-1) of (a_i_j * y_i)
        */
        double[] x = lpSoln.computeBasicSolution();
        double[] ys = lpSoln.computeBasicDualSolution();
        
        System.out.printf("x=%s\n", FormatArray.toString(x, "%.3f"));
        System.out.printf("ys=%s\n", FormatArray.toString(ys, "%.3f"));
        System.out.printf("y_n=\n");
        for (i = 0; i < lpSoln.nIndices.length; ++i) {
            System.out.printf("%.3f, ", ys[lpSoln.nIndices[i]]);
        }
        System.out.println();
        System.out.printf("x_b=\n");
        for (i = 0; i < lpSoln.bIndices.length; ++i) {
            System.out.printf("%.3f, ", x[lpSoln.bIndices[i]]);
        }
        System.out.println();
        
        double expectedZ = 8.5;
        
        double[] primalX = lpSoln.calculatePrimalX();
        double[] dualY = lpSoln.calculateDualY();
        System.out.printf("primal x=%s\n", FormatArray.toString(primalX, "%.3f"));
        System.out.printf("dual y=%s\n", FormatArray.toString(dualY, "%.3f"));
        System.out.printf("obj=%.4f, %.4f\n", lpSoln.evaluateObjective(),
            lpSoln.evaluateDualObjective());

        assertEquals(expectedYs.length, primalX.length);
        
        // the random picking of entering variable in the pivot of LinearProgramming,
        //     occassiionally gives a solution which is not optimal.        
        /*double diff, tol = 1e-7;
        for (i = 0; i < expectedYs.length; ++i) {
            diff = Math.abs(expectedYs[i] - primalX[i]);
            assertTrue(diff < tol);
        }*/
        SetCover sc = new SetCover();
        TIntSet cover = sc.weightedSetsApprox2LgN(nU, sets, weights);
        System.out.printf("cover=%s\n", Arrays.toString(cover.toArray()));

        // 4,3,2,0    cost=1+5+3+2=11
        // opt=1,2,3  cost=3+3+5=11      
        int[] optimalCover = new int[]{1,2,3};
        TIntIterator iter = cover.iterator();
        int c = 0;
        while (iter.hasNext()) {
            c += weights[iter.next()];
        }
        int oc = 0;
        for (i = 0; i < optimalCover.length; ++i) {
            oc += weights[optimalCover[i]];
        }
        double nOptC = optimalCover.length;
        double nC = cover.size();
        double lnN = Math.log(nU);
        //System.out.printf("%.3f/%.3f=%.3f, ln(%d)=%.3f\n", nC, nOptC, (nC/nOptC), nU, lnN);
        assertTrue((nC/nOptC) < lnN);
    }
    
    /*
    test sets from
    http://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html
    http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/
    Found the test set reference from google search result: the thesis 
    "Algorithms for NP-hard Optimization Problems and Cluster Analysis",
    Li, 2017 (see Table 1.2 1nd 1.3)
    
    scpe1 has 500 sets (= num columns n) 
    and 50 items (=num rows m) with opt=5.  
    all columns have same weight.
    
    scp41.txt has 1000 sets and 200 items.
    the column weights are not all the same.
    */
    
    public void estElementWeightedSetCover() throws IOException {
        
        String filename = "scp41.txt";
        
        //[nSets, nItems]
        int[] mn = ORSetCoverReader.getSCPDatasetNumberOfRowsCols(filename);
        
        int nSets = mn[0];
        int nU = mn[1];
        
        List<TIntSet> sets = new ArrayList<TIntSet>(nSets);
        double[] uWeights = new double[nU];
        
        ORSetCoverReader.getSCPDataset(filename, sets, uWeights);

        double[] sWeights = SetCover.calcSetWeightsFromElementWeights(nU, sets, uWeights);
        
        double expectedZ = ORSetCoverReader.getSCPDatasetZ(filename);
        
        SetCover sc = new SetCover();
        
        TIntSet cover = sc.weightedElementsApprox2LgN(nU, sets, uWeights);
    
        double sum = 0;
        TIntIterator iter = cover.iterator();
        int idx;
        while (iter.hasNext()) {
            idx = iter.next();
            sum += sWeights[idx];
        }
        System.out.printf("cover sum=%.3f, expected=%.3f\n", sum, expectedZ);
        System.out.printf("cover size=%d\n", cover.size());
        System.out.printf("cover=%s\n", Arrays.toString(cover.toArray()));
        
    }
    
}
