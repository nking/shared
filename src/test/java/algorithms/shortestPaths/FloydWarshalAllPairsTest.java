package algorithms.shortestPaths;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FloydWarshalAllPairsTest extends TestCase {

    public FloydWarshalAllPairsTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void test0() {
        
        int[][] edgeWeights = new int[5][5];
        for (int i = 0; i < edgeWeights.length; i++) {
            edgeWeights[i] = new int[5];
            Arrays.fill(edgeWeights[i], Integer.MAX_VALUE);
        }
        
        /*           < 1 >
         *          * | | *    
         *        3   | |   4
         *      /    |  |     \
         *  < 0 > --/----\-8--*< 2 >
         *    \ *   |     |      *
         *     \  2 |     1|    /
         *     -4   | \ _  |  -5
         *       * *      \ * /
         *      < 4 >--6--*< 3 >
         * 
         */
        
        edgeWeights[0][1] = 3;
        edgeWeights[0][2] = 8;
        edgeWeights[0][4] = -4;
        
        edgeWeights[1][3] = 1;
        edgeWeights[1][4] = 7;
        
        edgeWeights[2][1] = 4;
        
        edgeWeights[3][0] = 2;
        edgeWeights[3][2] = -5;
        
        edgeWeights[4][3] = 6;
        
        for (int i = 0; i < edgeWeights.length; i++) {
            edgeWeights[i][i] = 0;
        }
        
        FloydWarshalAllPairs apsp = new FloydWarshalAllPairs();
        
        apsp.findShortestPaths(edgeWeights);
        
        int[][] sps = apsp.dist;
        
        assertTrue(sps[0][0] == 0);
        assertTrue(sps[0][1] == 1);
        assertTrue(sps[0][2] == -3);
        assertTrue(sps[0][3] == 2);
        assertTrue(sps[0][4] == -4);
        
        assertTrue(sps[1][0] == 3);
        assertTrue(sps[1][1] == 0);
        assertTrue(sps[1][2] == -4);
        assertTrue(sps[1][3] == 1);
        assertTrue(sps[1][4] == -1);
       
        assertTrue(sps[2][0] == 7);
        assertTrue(sps[2][1] == 4);
        assertTrue(sps[2][2] == 0);
        assertTrue(sps[2][3] == 5);
        assertTrue(sps[2][4] == 3);
        
        assertTrue(sps[3][0] == 2);
        assertTrue(sps[3][1] == -1);
        assertTrue(sps[3][2] == -5);
        assertTrue(sps[3][3] == 0);
        assertTrue(sps[3][4] == -2);
        
        assertTrue(sps[4][0] == 8);
        assertTrue(sps[4][1] == 5);
        assertTrue(sps[4][2] == 1);
        assertTrue(sps[4][3] == 6);
        assertTrue(sps[4][4] == 0);
    
        int[][] p = apsp.prev;
        int[][] expectedP = new int[5][5];
        int nil = Integer.MIN_VALUE;
        expectedP[0] = new int[]{nil, 3-1, 4-1, 5-1, 1-1};
        expectedP[1] = new int[]{4-1, nil, 4-1, 2-1, 1-1};
        expectedP[2] = new int[]{4-1, 3-1, nil, 2-1, 1-1};
        expectedP[3] = new int[]{4-1, 3-1, 4-1, nil, 1-1};
        expectedP[4] = new int[]{4-1, 3-1, 4-1, 5-1, nil};
        assertEquals(expectedP.length, p.length);
        
        for (int i = 0; i < expectedP.length; ++i) {
            assertEquals(expectedP[i].length, p[i].length);
            assertTrue(Arrays.equals(expectedP[i], p[i]));
        }
    }
    
    public void est1() {
        
        int[][] edgeWeights = new int[5][5];
        for (int i = 0; i < edgeWeights.length; i++) {
            edgeWeights[i] = new int[5];
            Arrays.fill(edgeWeights[i], Integer.MAX_VALUE);
        }
        
        /*    < 0 >*<-- 3 -- < 1 >
         *      |  \            *  *
         *      |    \         /\    \ 4
         *    6 |      3       |       \
         *      |        \     1      < 4 >
         *      |          \   |       / 
         *      *      1     * |    * 2
         *    < 2 >*-------- < 3 >
         *          --------*
         *             2
         */
        
        edgeWeights[0][2] = 6;
        edgeWeights[0][3] = 3;
        edgeWeights[1][0] = 3;
        edgeWeights[2][3] = 2;
        edgeWeights[3][1] = 1;
        edgeWeights[3][2] = 1;        
        edgeWeights[4][1] = 4;
        edgeWeights[4][3] = 2;
        
        for (int i = 0; i < edgeWeights.length; i++) {
            edgeWeights[i][i] = 0;
        }
        
        FloydWarshalAllPairs apsp = new FloydWarshalAllPairs();
        apsp.setDebug(true);
        
        apsp.findShortestPaths(edgeWeights);
        
        int[][] sps = apsp.dist;
        
        assertTrue(sps.length == 5);
        assertTrue( Arrays.equals(sps[0], new int[]{0, 4, 4, 3, Integer.MAX_VALUE}));
        assertTrue( Arrays.equals(sps[1], new int[]{3, 0, 7, 6, Integer.MAX_VALUE}));
        assertTrue( Arrays.equals(sps[2], new int[]{6, 3, 0, 2, Integer.MAX_VALUE}));
        assertTrue( Arrays.equals(sps[3], new int[]{4, 1, 1, 0, Integer.MAX_VALUE}));
        assertTrue( Arrays.equals(sps[4], new int[]{6, 3, 3, 2, 0}));
    
        int[][] prev = apsp.prev;
        assertTrue(prev.length == 5);
        assertTrue( Arrays.equals(prev[0], new int[]{Integer.MIN_VALUE, 3, 3, 0, Integer.MIN_VALUE}));
        assertTrue( Arrays.equals(prev[1], new int[]{1, Integer.MIN_VALUE, 3, 0, Integer.MIN_VALUE}));
        assertTrue( Arrays.equals(prev[2], new int[]{1, 3, Integer.MIN_VALUE, 2, Integer.MIN_VALUE}));
        assertTrue( Arrays.equals(prev[3], new int[]{1, 3, 3, Integer.MIN_VALUE, Integer.MIN_VALUE}));
        assertTrue( Arrays.equals(prev[4], new int[]{1, 3, 3, 4, Integer.MIN_VALUE}));
    }

}
