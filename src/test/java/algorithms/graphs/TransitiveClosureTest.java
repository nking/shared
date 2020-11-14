package algorithms.graphs;

import algorithms.VeryLongBitString;
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
        
        int n = c.length;
        
        /*System.out.println("input:");
        for (int i = 0; i < n; i++) {
            System.out.println("c i=" + i + " : " + Arrays.toString(c[i]));
        }
        System.out.println("output:");
        for (int i = 0; i < n; i++) {
            System.out.println("e i=" + i + " : " + Arrays.toString(e[i]));
        }*/
        
        TransitiveClosure tc = new TransitiveClosure();
        //tc.setDebug(true);
        boolean[][] t = tc.calc(c);
        
        assertEquals(e.length, t.length);
        for (int i = 0; i < e.length; ++i) {
            assertTrue(Arrays.equals(e[i], t[i]));
        }
        
        
        //tc.setDebug(true);
        VeryLongBitString[] cBS = TransitiveClosure.convert(c);
        VeryLongBitString[] eBS = TransitiveClosure.convert(e);
        /*System.out.println("input:");
        for (int i = 0; i < n; i++) {
            System.out.println("t i=" + i + " : " + cBS[i].toString());
        }
        System.out.println("output:");
        for (int i = 0; i < n; i++) {
            System.out.println("t i=" + i + " : " + eBS[i].toString());
        }*/
        //tc.setDebug(true);
        VeryLongBitString[] tBS = tc.calc(cBS);
        assertEquals(e.length, t.length);
        for (int i = 0; i < e.length; ++i) {
            assertTrue(Arrays.equals(eBS[i].getSetBits(), tBS[i].getSetBits()));
        }
    }
}
