package algorithms.util;

import junit.framework.TestCase;

import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class PolyIntTest extends TestCase {

    public PolyIntTest() {
    }
    
    public void test0() {
        int[] a = new int[]{0, 10, 20, 44, 12};
        int[] a2 = Arrays.copyOf(a, a.length);
        int[] b = new int[]{0, 10, 20, 44};
        int[] c = new int[]{10, 20, 0, 12, 44};
        int[] d = new int[a.length];

        assertEquals( (new PolyInt(a)), (new PolyInt(a2)));
        assertTrue( (new PolyInt(a)).equals(new PolyInt(a2)));
        assertFalse( (new PolyInt(a)).equals(new PolyInt(b)));
        assertFalse( (new PolyInt(a)).equals(new PolyInt(c)));
        assertFalse( (new PolyInt(a)).equals(new PolyInt(d)));
    }
    
}
