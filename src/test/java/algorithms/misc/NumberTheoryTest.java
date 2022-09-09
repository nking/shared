package algorithms.misc;

import java.math.BigInteger;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class NumberTheoryTest extends TestCase {

    public NumberTheoryTest(String testName) {
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
            int result = NumberTheory.euclid(a, b);
            System.out.println(" result=" + result);
            if (result == 1) {
                // a is possibly prime
            } else {
                // a is not prime
            }
        }
    }

    /**
     * Test of euclid method, of class NumberTheory.
     */
    public void testEuclid() {
        System.out.println("testEuclid");
        int a = 30;
        int b = 21;
        int expResult = 3;
        int result = NumberTheory.euclid(a, b);
        assertTrue(expResult == result);
        
        result = NumberTheory.euclid(b, a);
        assertTrue(expResult == result);
        
        a = 561;
        b = 21;
        expResult = 3;
        //561 = 3*11*17:
        result = NumberTheory.euclid(a, b);
        assertTrue(expResult == result);
    }

    public void testExtendedEuclid() {
        //from CLRS Fig 31.1 (Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms)
        System.out.println("testExtendedEuclid");
        long a = 99;
        long b = 78;
        long[] expResult = new long[]{3, -11, 14};
        long[] result = NumberTheory.extendedEuclid(a, b);
        assertTrue( Arrays.equals(expResult, result));
        
        result = NumberTheory.extendedEuclid(78, 21);
        expResult = new long[]{3, 3, -11};
        assertTrue( Arrays.equals(expResult, result));
        
        result = NumberTheory.extendedEuclid(21, 15);
        expResult = new long[]{3, -2, 3};
        assertTrue( Arrays.equals(expResult, result));
        
        result = NumberTheory.extendedEuclid(15, 6);
        expResult = new long[]{3, 1, -2};
        assertTrue( Arrays.equals(expResult, result));
        
        result = NumberTheory.extendedEuclid(6, 3);
        expResult = new long[]{3, 0, 1};
        assertTrue( Arrays.equals(expResult, result));
        
        result = NumberTheory.extendedEuclid(3, 0);
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
            euclid = NumberTheory.euclid(a, n);
            euclid1 = NumberTheory.euclid(z0, n);
            dxy = NumberTheory.extendedEuclid(a, n);
            axny = a*dxy[1] + n*dxy[2];
            System.out.printf("a=%d, z=%d, gcd=(%d,%d) (gcd, x, y) = (%s) a*x+n*y=%d\n", 
                a, z0, euclid, euclid1, 
                Arrays.toString(dxy), axny);
        }
        
        System.out.printf("ee(%d,%d)=%s\n",
            2,5, Arrays.toString(NumberTheory.extendedEuclid(2, 5)));
        System.out.printf("ee(%d,%d)=%s\n",
            3,13, Arrays.toString(NumberTheory.extendedEuclid(3, 13)));
    }
    
    public void testGcdModularLinearEqnSolver() {
        // example form Cormen, Leiserson, Rivest, and Stein Introuduction to Algorithms, Sect 31.4
        long a = 14;
        long b = 30;
        long n = 100;
        
        long[] dXY = NumberTheory.extendedEuclid(a, n);
        assertEquals(2L, dXY[0]);
        assertEquals(-7L, dXY[1]);
        assertEquals(1L, dXY[2]);
                
        long[] s = NumberTheory.gcdModularLinearEqnSolver(a, b, n);
        assertEquals(95L, s[0]);
        assertEquals(45L, s[1]);
    }
    
    public void testModularExponentiation() {
        // test from CLRS (Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms)
        // Chap 31
        int a = 7;
        int b = 560; 
        int n = 561;
        /*
        i   9   8   7   6   5   4   3    2   1   0 
        bi  1   0   0   0   1   1   0    0   0   0
         c  1   2   4   8  17  35  70  140 280 560 
         d  7  49 157 526 160 241 298  166  67   1
        */
        int m = NumberTheory.modularExponentiation(a, b, n);
        assertEquals(1, m);
    }
    
    public void testRSA() {
        // following Section 31.7 of Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms (CLRS)
        // exercise 31.7-1
        
        // NOTE: see also https://en.wikipedia.org/wiki/Paillier_cryptosystem
        
        int p = 11;
        int q = 29;
        int n = p*q;//319
        int phi = (p-1)*(q-1);//p*q*(p-1)*(q-1); 280
        long e = 3;
        long d = NumberTheory.euclid(e, phi);
        System.out.printf("phi=%d, e=%d gcd=%d\n", phi, e, d);
        assertEquals(1, d);
        
        // ** step 3 in RSA **
        // solve for e:
        int a = 1;
        long[] x = null;
        
        long eMultInv = 0;
        
        for (int b = phi+1; b < phi+10; ++b) {
            if ((b & 1) == 0) {
                continue;
            }
        
            //solves for x in the equation a * x ≡ b (mod n) as 'e'
            x = NumberTheory.gcdModularLinearEqnSolver(a, b, phi);
            System.out.printf("phi=%d, a=%d, b=%d, x=%s\n",
                phi, a, b, Arrays.toString(x));
            
            if (x.length == 1 && x[0] > 1) {
                // calculate the multiplicative inverse of e, modulo phi(n)                

                long[] dxy = NumberTheory.extendedEuclid(e, phi);
                System.out.printf("extEuc(e,phi) dXY=%s\n", Arrays.toString(dxy));
                assertEquals(1, dxy[0]);
                // finds [1,-93,1] so x=-93
                //        noticing -93 is equal to -(phi/e)
                // ɸ(n)=280
                // e*x mod phi = 1;  yes, (3*-93) % 280 = 1
                // the multiplicative inverse of e is -93 % 280 = 187 (noticing it's = -(phi/e) mod phi)
                eMultInv = Math.floorMod(dxy[1], phi);
                assertEquals(187, eMultInv);
                break;
            }  
        } 
        assertNotNull(x);
        assertEquals(e, x[0]);
        
        // public pair = (e, n)
        // private pair = (eMultInv, n)
        System.out.printf("public pair=(%d,%d)\n", e, n);
        System.out.printf("private pair=(%d,%d)\n", eMultInv, n);
        
        int m = 100;
        
        //P(M) = M^e mod n = 100^3 mod 319 = 254
        long mPublicEncryp = Math.floorMod((long)Math.pow(m, e), n);
        assertEquals(254, mPublicEncryp);
        
        //apply private key:
        // //(M^e mod n)^d mod n = M
        // this is incorrect due to overflow.  need to use big integer
        //long mDecryp = Math.floorMod((long)Math.pow(mPublicEncryp, eMultInv), n); //40
        //System.out.printf("m=%d, P(m)=%d, S(P(m))=%d\n", m, mPublicEncryp, mDecryp);
        BigInteger c = new BigInteger(Integer.toString((int)mPublicEncryp));
        BigInteger cPowD = c.pow((int)eMultInv);
        BigInteger mDecryp2 = cPowD.mod(new BigInteger(Integer.toString(n)));
        System.out.printf("biginteger: m=%d, P(m)=%s, S(P(m))=%s\n", m, c.toString(), 
            mDecryp2.toString());
        assertEquals(Integer.toString(m), mDecryp2.toString());
    }
    
    public void testLCM() {
        //https://en.m.wikipedia.org/wiki/Lowest_common_denominator
        
        assertEquals(6, NumberTheory.leastCommonMultiple(2, 3));
        
        assertEquals(36, NumberTheory.leastCommonMultiple(12, 18));
    }
}
