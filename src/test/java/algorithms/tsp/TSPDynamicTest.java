package algorithms.tsp;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TSPDynamicTest extends TestCase {
    
    private int n = 6;
    private TSPDynamic tsp;
    
    public TSPDynamicTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        int[][] dist = new int[n][];
        for (int i = 0; i < n; ++i) {
            dist[i] = new int[n];
        }
        // set 0:4 and 3:0 in dist to 412 and 53:
        dist[0][4] = 412;
        dist[3][0] = 53;
        
        tsp = new TSPDynamic(dist);
    }
    
    public void testSetBits() {
        // for n=6, there are 5 nodes that aren't the startnode 0: 1-5
        // w = 4
        //4 = 100
        //1 = 001
        //2 = 010
        //5 = 101
        //3 = 011
        // path of 4,1,2 = 010 001 100
        
        //path of 4,1,2,5,3 =  011  101 010 001 100
        
        int nSet, nUnset;
        int base10Node = 4;
        long path = 0;
        int pathNodeNumber = 0;
        long path4 = tsp.setBits(base10Node, path, pathNodeNumber);
        assertEquals(path4, 4);
        
        int b10 = tsp.getBase10NodeIndex(pathNodeNumber, path4);
        assertEquals(base10Node, b10);
        
        assertEquals(n-1, tsp.numberOfUnsetNodes(path));
        assertEquals(0, tsp.numberOfSetNodes(path));
        
        assertEquals(n-2, tsp.numberOfUnsetNodes(path4));
        assertEquals(1, tsp.numberOfSetNodes(path4));
        
        
        base10Node = 1;
        pathNodeNumber++;
        long path41 = tsp.setBits(base10Node, path4, pathNodeNumber);
        assertEquals(path41, 12);
        
        b10 = tsp.getBase10NodeIndex(pathNodeNumber, path41);
        assertEquals(base10Node, b10);
        
        base10Node = 2;
        pathNodeNumber++;
        long path412 = tsp.setBits(base10Node, path41, pathNodeNumber);
        assertEquals(path412, 140);
        
        b10 = tsp.getBase10NodeIndex(pathNodeNumber, path412);
        assertEquals(base10Node, b10);
        
        int[] s412 = new int[]{4, 1, 2};
        long path412a = tsp.createThe3NodeBitstring(s412);
        assertEquals(path412, path412a);
        
        int nSet412 = tsp.numberOfSetNodes(path412);
        assertEquals(3, nSet412);
        
        long path412b = tsp.concatenate(0, 0, s412);
        assertEquals(path412, path412b);
        
        //path of 4,1,2,5,3 =  011 101 010 001 100
        int[] s53 = new int[]{5, 3};
        long path41253 = tsp.concatenate(path412, nSet412, s53);
        assertEquals(14988, path41253);
        
        // set 0:4 and 3:0 in dist to 412 and 53 in setUp:
        tsp.compareToMin(path41253, 0);
        long expected = 412 + 53;
        assertEquals(expected, tsp.getMinCost());
    }
    
}
