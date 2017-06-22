package algorithms.util;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MiscSorterTest extends TestCase {

    public void testSortByFirstArgument_1() {

        int[] a = new int[]{3, 7, 1, 5};
        int[] b = new int[]{0, 1, 2, 3};
        
        MiscSorter.sortBy1stArg(a, b);
        
        assertTrue(a[0] == 1);
        assertTrue(b[0] == 2);
        
        assertTrue(a[1] == 3);
        assertTrue(b[1] == 0);
        
        assertTrue(a[2] == 5);
        assertTrue(b[2] == 3);
        
        assertTrue(a[3] == 7);
        assertTrue(b[3] == 1); 
    }
 
}
