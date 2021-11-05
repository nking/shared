package algorithms.optimization;

import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import junit.framework.TestCase;
import thirdparty.HungarianAlgorithm;

/**
 *
 * @author nichole
 */
public class TSPChristofidesSerdyukovTest extends TestCase {
    
    public TSPChristofidesSerdyukovTest(String testName) {
        super(testName);
    }

    /**
     * Test of approxTSPTour method, of class TSPChristofidesSerdyukov.
     */
    public void estApproxTSPTour() {
        System.out.println("approxTSPTour");
        int nVertexes = 0;
        TIntObjectMap<TIntIntMap> adjCostMap = null;
        TSPChristofidesSerdyukov instance = new TSPChristofidesSerdyukov();
        int[] expResult = null;
        int[] result = instance.approxTSPTour(nVertexes, adjCostMap);
        //assertEquals(expResult, result);
    }
    
    public void testCalculateDegrees() {
        
        int nVertexes = 8;
        Map<Integer, LinkedList<Integer>> mstTree = new HashMap<Integer, LinkedList<Integer>>();
        int i, j;
        LinkedList<Integer> list;
        TIntObjectMap<TIntIntMap> adjCostMap = new TIntObjectHashMap<TIntIntMap>();        
        for (i = 0; i < nVertexes; ++i) {
            list = new LinkedList<Integer>();
            mstTree.put(i, list);
            adjCostMap.put(i, new TIntIntHashMap());
        }
        /* make adjacency list 
                    0
              *1           *2
            3   *4       *5   *6
          *7  
        */
        mstTree.get(0).add(1); mstTree.get(0).add(2);
        mstTree.get(1).add(0); mstTree.get(1).add(3); mstTree.get(1).add(4); //*
        mstTree.get(2).add(0); mstTree.get(2).add(5); mstTree.get(2).add(6); //*
        mstTree.get(3).add(1); mstTree.get(3).add(7);
        mstTree.get(4).add(1); //*
        mstTree.get(5).add(2); //*
        mstTree.get(6).add(2); //*
        mstTree.get(7).add(3); //*
        
        int costHigh = 100; int costLow = 2;
        adjCostMap.get(0).put(1, costLow+2); adjCostMap.get(0).put(2, costLow+2);
        adjCostMap.get(1).put(0, costLow+2); adjCostMap.get(1).put(3, costLow); adjCostMap.get(1).put(4, costLow);
        adjCostMap.get(2).put(0, costLow+2); adjCostMap.get(2).put(5, costLow); adjCostMap.get(2).put(6, costLow+1);
        adjCostMap.get(3).put(1, costLow); adjCostMap.get(3).put(7, costLow);
        adjCostMap.get(4).put(1, costLow); 
        adjCostMap.get(5).put(2, costLow);
        adjCostMap.get(6).put(2, costLow+1);
        adjCostMap.get(7).put(3, costLow);
        
        
        int[] expecctedD = new int[]{2, 3, 3, 2, 1, 1, 1, 1};
        
        TSPChristofidesSerdyukov tsp = new TSPChristofidesSerdyukov();
        
        int[] d = tsp.calculateDegrees(mstTree, nVertexes);

        assertTrue(Arrays.equals(expecctedD, d));
        
        int[] oddDExpected = new int[]{1, 2, 4, 5, 6, 7};
        
        int[] oddD = tsp.oddPassFilter(d);
        
        assertEquals(oddDExpected.length, oddD.length);

        assertTrue(Arrays.equals(oddDExpected, oddD));
        
        
        float[][] expectedCostMatrix = new float[oddDExpected.length][];
        for (i = 0; i < expectedCostMatrix.length; ++i) {
            expectedCostMatrix[i] = new float[oddDExpected.length];
            Arrays.fill(expectedCostMatrix[i], Float.MAX_VALUE);
        }
        //int[] oddDExpected = new int[]{1, 2, 4, 5, 6, 7};
        //                               0  1  2  3  4  5
        //1:4 is cost 2 => 0:2
        //2:5 is cost 2 => 1:3
        //2:6 is cost 3 => 1:4
        expectedCostMatrix[0][2] = costLow;
        expectedCostMatrix[1][3] = costLow;
        expectedCostMatrix[1][4] = costLow + 1;
        
        expectedCostMatrix[2][0] = costLow;
        expectedCostMatrix[3][1] = costLow;
        expectedCostMatrix[4][1] = costLow + 1;
        
        float[][] costMatrix = tsp.buildCostMatrix(oddD, adjCostMap);
        assertEquals(expectedCostMatrix.length, costMatrix.length);
        assertEquals(expectedCostMatrix[0].length, costMatrix[0].length);
        double diff;
        double tol = 1e-7;
        for (i = 0; i < expectedCostMatrix.length; ++i) {
            for (j = 0; j < expectedCostMatrix[i].length; ++j) {
                diff = Math.abs(expectedCostMatrix[i][j] - costMatrix[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        // expecting 1:4 and 2:5 to be matched which are 0:2 and 1:3 before 
        //    transformed back to original indexes
        HungarianAlgorithm ha = new HungarianAlgorithm();
        int[][] assignments = ha.computeAssignments(costMatrix);
        
        // use min cost unbalanced here too

        int z = 1;

    }
    
}
