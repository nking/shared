package algorithms.graphs;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TransitiveClosureTest extends TestCase {
    
    public TransitiveClosureTest(String testName) {
        super(testName);
    }
    
    public void testCalc() {
        
        boolean[][] c = new boolean[4][4];
        c[0] = new boolean[]{true, false, false, false};
        c[1] = new boolean[]{false, true, true, true};
        c[2] = new boolean[]{false, true, true, false};
        c[3] = new boolean[]{true, false, true, true};
        
        boolean[][] e = new boolean[4][4];
        e[0] = new boolean[]{true, false, false, false};
        e[1] = new boolean[]{true, true, true, true};
        e[2] = new boolean[]{true, true, true, true};
        e[3] = new boolean[]{true, true, true, true};
        
        TransitiveClosure tc = new TransitiveClosure();
        //tc.setDebug(true);
        tc.calc(c);
        
        boolean[][] t = tc.t;
        assertEquals(e.length, t.length);
        for (int i = 0; i < e.length; ++i) {
            assertTrue(Arrays.equals(e[i], t[i]));
        }
        
    }
}
