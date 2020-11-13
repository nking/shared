package algorithms.sort;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class InsertionSortTest extends TestCase {
    
    public InsertionSortTest(String testName) {
        super(testName);
    }
    
    public void testSort(){
        /*
        from Cormen et al. "Intorudction toAlgorithms" Fig 2.2
        */
        
        int[] a = new int[]{5, 2, 4, 6, 1, 3};
        
        int[] expected = new int[]{1, 2, 3, 4, 5, 6};
        
        InsertionSort.sort(a);
        
        for (int i = 0; i < expected.length; ++i) {
            assertEquals(expected[i], a[i]);
        }
    }
    
}
