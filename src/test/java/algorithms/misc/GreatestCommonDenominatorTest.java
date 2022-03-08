package algorithms.misc;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GreatestCommonDenominatorTest extends TestCase {

    public GreatestCommonDenominatorTest(String testName) {
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
    
    public void testEuclid0() {
        System.out.println("testEuclid0");
        int a = 11;
        for (int b = (a-1); b > 0; b--) {
            int result = GreatestCommonDenominator.euclid(a, b);
            System.out.println(" result=" + result);
            if (result == 1) {
                // a is possibly prime
            } else {
                // a is not prime
            }
        }
    }

    /**
     * Test of euclid method, of class GreatestCommonDenominator.
     */
    public void testEuclid() {
        System.out.println("testEuclid");
        int a = 30;
        int b = 21;
        int expResult = 3;
        int result = GreatestCommonDenominator.euclid(a, b);
        assertTrue(expResult == result);
        
        result = GreatestCommonDenominator.euclid(b, a);
        assertTrue(expResult == result);
        
        a = 561;
        b = 21;
        expResult = 3;
        //561 = 3*11*17:
        result = GreatestCommonDenominator.euclid(a, b);
        assertTrue(expResult == result);
    }

    public void testExtendedEuclid() {
        //from CLRS Fig 31.1 (Cormen et al. Introduction to Algorithms)
        System.out.println("testExtendedEuclid");
        long a = 99;
        long b = 78;
        long[] expResult = new long[]{3, -11, 14};
        long[] result = GreatestCommonDenominator.extendedEuclid(a, b);
        assertTrue( Arrays.equals(expResult, result));
        
        result = GreatestCommonDenominator.extendedEuclid(78, 21);
        expResult = new long[]{3, 3, -11};
        assertTrue( Arrays.equals(expResult, result));
        
        result = GreatestCommonDenominator.extendedEuclid(21, 15);
        expResult = new long[]{3, -2, 3};
        assertTrue( Arrays.equals(expResult, result));
        
        result = GreatestCommonDenominator.extendedEuclid(15, 6);
        expResult = new long[]{3, 1, -2};
        assertTrue( Arrays.equals(expResult, result));
        
        result = GreatestCommonDenominator.extendedEuclid(6, 3);
        expResult = new long[]{3, 0, 1};
        assertTrue( Arrays.equals(expResult, result));
        
        result = GreatestCommonDenominator.extendedEuclid(3, 0);
        expResult = new long[]{3, 1, 0};
        assertTrue( Arrays.equals(expResult, result));
        
        
    }

    public void testZStarGenerator() {
        int z0;
        int a;
        int k = 1;
        int n = 7;//15;
        int euclid;
        int euclid1;
        long[] dxy;
        long axny;
        for (a = 0; a < n; ++a) {
            z0 = ((a + k*n) % n);
            euclid = GreatestCommonDenominator.euclid(a, n);
            euclid1 = GreatestCommonDenominator.euclid(z0, n);
            dxy = GreatestCommonDenominator.extendedEuclid(a, n);
            axny = a*dxy[1] + n*dxy[2];
            System.out.printf("a=%d, z=%d, gcd=(%d,%d) (gcd, x, y) = (%s) a*x+n*y=%d\n", 
                a, z0, euclid, euclid1, 
                Arrays.toString(dxy), axny);
        }
    }
    
    public void testGcdModularLinearEqnSolver() {
        // example form Cormen et al. Introuduction to Algorithms, Sect 31.4
        long a = 14;
        long b = 30;
        long n = 100;
        
        long[] dXY = GreatestCommonDenominator.extendedEuclid(a, n);
        assertEquals(2L, dXY[0]);
        assertEquals(-7L, dXY[1]);
        assertEquals(1L, dXY[2]);
                
        long[] s = GreatestCommonDenominator.gcdModularLinearEqnSolver(a, b, n);
        assertEquals(95L, s[0]);
        assertEquals(45L, s[1]);
    }
}
