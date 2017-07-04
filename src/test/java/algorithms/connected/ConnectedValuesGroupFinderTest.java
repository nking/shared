package algorithms.connected;

import gnu.trove.set.TIntSet;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ConnectedValuesGroupFinderTest extends TestCase {
    
    public ConnectedValuesGroupFinderTest(String testName) {
        super(testName);
    }
    
    public void test() {
        
        int w = 5;
        int h = 5;
        
        int[][] data = new int[w][];
        data[0] = new int[]{0, 0, 0, 0, 0};
        data[1] = new int[]{0, 1, 1, 3, 0};
        data[2] = new int[]{0, 1, 2, 3, 0};
        data[3] = new int[]{0, 2, 3, 1, 0};
        data[4] = new int[]{0, 0, 0, 0, 0};
        
        ConnectedValuesGroupFinder finder = new ConnectedValuesGroupFinder();
        finder.setMinimumNumberInCluster(2);
        finder.setToUse8Neighbors();
        
        List<TIntSet> groups = finder.findGroups(data);
        
        assertEquals(4, groups.size());
        
        
        //-------
        /*
        data[0] = new int[]{0, 0, 0, 0, 0};
        data[1] = new int[]{0, 1, 1, 3, 0};
        data[2] = new int[]{0, 1, 2, 3, 0};
        data[3] = new int[]{0, 2, 3, 1, 0};
        data[4] = new int[]{0, 0, 0, 0, 0};
        */
        finder = new ConnectedValuesGroupFinder();
        finder.setMinimumNumberInCluster(2);
        assertTrue(finder.use4Neighbors);
        
        groups = finder.findGroups(data);
        
        assertEquals(3, groups.size());
    }
}
