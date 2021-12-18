package algorithms.graphs;

import algorithms.graphs.ApproxGraphSearchZeng.Graph;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class StarStructureTest extends TestCase {
    
    public StarStructureTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
        List<Graph> dbs = new ArrayList<Graph>();
        Graph q = ApproxGraphSearchZengTest.getG0(dbs);
        int nV = 14; // dB has 15, q has 14
        
        StarStructure[] stars = StarStructure.createStarStructureMultiset(q);

        assertEquals(nV, stars.length);
        
        int ged;
        
        ged = StarStructure.calculateEditDistanceV(stars[0], stars[1]);
        assertEquals(0, ged);
        ged = StarStructure.calculateEditDistance(stars[0], stars[1]);
        assertEquals(0, ged);
        
        ged = StarStructure.calculateEditDistanceV(stars[0], stars[8]);
        assertEquals(1, ged);
        ged = StarStructure.calculateEditDistance(stars[0], stars[8]);
        assertEquals(1, ged);
        
        int s3 = StarStructure.calculateSupport(stars, 3);
        assertEquals(6, s3);
        int s4 = StarStructure.calculateSupport(stars, 4);
        assertEquals(3+3+3, s4);
        
        s4 = StarStructure.calculateSupport(
            StarStructure.createStarStructureMultiset(dbs.get(0)), 4);
        assertEquals(4+3+3, s4);
        
        // ==== test a few properties of copy
        StarStructure[] stars2 = StarStructure.copy(stars);
        assertEquals(stars.length, stars2.length);
        
        ged = StarStructure.calculateEditDistanceV(stars2[0], stars2[1]);
        assertEquals(0, ged);
        ged = StarStructure.calculateEditDistance(stars2[0], stars2[1]);
        assertEquals(0, ged);
        
        ged = StarStructure.calculateEditDistanceV(stars2[0], stars2[8]);
        assertEquals(1, ged);
        ged = StarStructure.calculateEditDistance(stars2[0], stars2[8]);
        assertEquals(1, ged);
        
        s3 = StarStructure.calculateSupport(stars2, 3);
        assertEquals(6, s3);
        s4 = StarStructure.calculateSupport(stars2, 4);
        assertEquals(3+3+3, s4);
        
    }

    /**
     * Test of calculateEditDistanceNoRelabeling method, of class StarStructure.
     */
    public void testCalculateEditDistanceNoRelabeling() {
        
    }

    /**
     * Test of calculateEditDistanceNoRelabelingV method, of class StarStructure.
     */
    public void testCalculateEditDistanceNoRelabelingV() {
        
    }

    /**
     * Test of createDistanceMatrix method, of class StarStructure.
     */
    public void testCreateDistanceMatrix() {
        
    }

    /**
     * Test of createDistanceMatrixV method, of class StarStructure.
     */
    public void testCreateDistanceMatrixV() {
        
    }

    /**
     * Test of createDistanceMatrixNoRelabeling method, of class StarStructure.
     */
    public void testCreateDistanceMatrixNoRelabeling() {
        
    }

    /**
     * Test of createDistanceMatrixNoRelabelingV method, of class StarStructure.
     */
    public void testCreateDistanceMatrixNoRelabelingV() {
        
    }

}
