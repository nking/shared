package thirdparty.graphMatchingToolkit.algorithms;

import algorithms.bipartite.Graph;
import algorithms.bipartite.GraphWithoutWeights;
import algorithms.bipartite.HopcroftKarp;
import algorithms.bipartite.MinCostUnbalancedAssignment3Test;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.logging.Logger;
import junit.framework.TestCase;
import thirdparty.HungarianAlgorithm;

/**
 *
 * @author nichole
 */
public class VolgenantJonkerTest extends TestCase {
    
    private Logger log = 
        Logger.getLogger(this.getClass().getName());
    
    public VolgenantJonkerTest(String testName) {
        super(testName);
    }
    
    public void test0() throws NoSuchAlgorithmException {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1466321947621L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        int size = 100;
        int maxCost = 10000;
        log.info("size=" + size + " maxCost=" + maxCost);
                
        Graph g = MinCostUnbalancedAssignment3Test.getTestGraph1(size, maxCost, sr);
        
        float[][] fMatrix = MinCostUnbalancedAssignment3Test.convert(g);
        double[][] dMatrix = MinCostUnbalancedAssignment3Test.convert(fMatrix);
                
        HungarianAlgorithm ha = new HungarianAlgorithm();
                
        int[][] haMatching = ha.computeAssignments(fMatrix);
        
        VolgenantJonker vj = new VolgenantJonker();
        double vJCost = vj.computeAssignment(dMatrix);
        int[] vjRowAssignments = vj.getAssignment();
        
        int i, j;
        /*for (i = 0; i < haMatching.length; ++i) {
            System.out.printf("HA: (%d, %d)\n", haMatching[i][0], haMatching[i][1]);
        }
        for (i = 0; i < vjRowAssignments.length; ++i) {
            System.out.printf("VJ: (%d, %d)\n", i, vjRowAssignments[i]);
        }*/
        int r, c;
        for (i = 0; i < haMatching.length; ++i) {
            r = haMatching[i][0];
            c = haMatching[i][1];
            assertEquals(vjRowAssignments[r], c);
        }
    }
    
    public void testUnweighted() throws NoSuchAlgorithmException {
        
        int size = 10;
        int maxCost = 10;
    
        /*
        - graph of size n for both sets
        - random number of edges for each
        - all edges have same cost
        --> expecting same results as hopcroft-karp, that
            is a maximal matching, but the random links
            may create less than maximal possibilities.
            the method is mostly to exercise the
            code to explore it further.
        */
                
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1464995162443L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        Graph g = MinCostUnbalancedAssignment3Test.getTestGraph3(sr, size, maxCost);
        
        HopcroftKarp hk = new HopcroftKarp();
        int[] hkMatched = hk.hopcroftKarpV0(new GraphWithoutWeights(g));
        
        float[][] fMatrix = MinCostUnbalancedAssignment3Test.convert(g);
        double[][] dMatrix = MinCostUnbalancedAssignment3Test.convert(fMatrix);
                
        HungarianAlgorithm ha = new HungarianAlgorithm();
                
        int[][] haMatching = ha.computeAssignments(fMatrix);
        
        VolgenantJonker vj = new VolgenantJonker();
        double vJCost = vj.computeAssignment(dMatrix);
        int[] vjRowAssignments = vj.getAssignment();
        
        int i, j;
        /*for (i = 0; i < haMatching.length; ++i) {
            System.out.printf("HA: (%d, %d)\n", haMatching[i][0], haMatching[i][1]);
        }
        for (i = 0; i < vjRowAssignments.length; ++i) {
            System.out.printf("VJ: (%d, %d)\n", i, vjRowAssignments[i]);
        }
        for (i = 0; i < vjRowAssignments.length; ++i) {
            System.out.printf("HK: (%d, %d)\n", i, hkMatched[i]);
        }*/
        assertEquals(hkMatched.length, fMatrix.length);
        assertEquals(haMatching.length, fMatrix.length);
        assertEquals(vjRowAssignments.length, fMatrix.length);
    }
    
}
