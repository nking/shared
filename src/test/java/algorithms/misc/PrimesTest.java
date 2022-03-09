
package algorithms.misc;

import gnu.trove.iterator.TLongIterator;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PrimesTest extends TestCase {
    
    public PrimesTest(String testName) {
        super(testName);
    }
    
    public void test0() throws NoSuchAlgorithmException {
        
        // Figure 31.7 of Cormen et al. Introduction to Alforithms
        
        long n;
        long[] expected;
        TLongSet result;
        long[] sorted;
        
        
        n = 1387;
        expected = new long[]{19, 73};
        result = Primes.pollardRhoFactorization(n);
        sorted = result.toArray();
        Arrays.sort(sorted);
        System.out.println("\n1387: " + Arrays.toString(sorted));
       // assertTrue(Arrays.equals(expected, sorted));
       
        //----
        n = 825;
        //  prime answer using LCDs = 11, 5, 3
        expected = new long[]{3, 5, 11};
        result = Primes.pollardRhoFactorization(n);
        sorted = result.toArray();
        Arrays.sort(sorted);
        System.out.println("\n825: " + Arrays.toString(sorted));
       // assertTrue(Arrays.equals(expected, sorted));
    }
    
    public void estEE() throws NoSuchAlgorithmException {
        
        long n;
        TLongSet expected;
        
        //----
        n = 825;
        expected = new TLongHashSet();
        // 11 and 75  OR 15 and 55
        //  ideally prime answer: 11, 5, 3
        expected.add(Long.valueOf(11));
        expected.add(Long.valueOf(75));
        //expected.add(Long.valueOf(11));
        //expected.add(Long.valueOf(5));
        //expected.add(Long.valueOf(3));
        long[] r;
        long a, b;
        a = n; b = 75;
        r = NumberTheory.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        b = n; a = 75;
        r = NumberTheory.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        a = n; b = 55;
        r = NumberTheory.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        b = n; a = 55;
        r = NumberTheory.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        a = n; b = 15;
        r = NumberTheory.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        b = n; a = 15;
        r = NumberTheory.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        a = n; b = 5;
        r = NumberTheory.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        b = n; a = 5;
        r = NumberTheory.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        a = n; b = 3;
        r = NumberTheory.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        b = n; a = 3;
        r = NumberTheory.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
    }
    
    public void estMillerRabin() throws Exception {
        long n = 561;
        int t = 4;
        int u = 35;
        long a = 7;
        
        //x0 ≡ a^(35) ≡ 241 (mod 561)
        // which computes the sequence X = <241, 298, 166, 67, 1>
        // so finds nontrivial square root of 1 in the last squaring step, 
        // since a^(280) ≡ 67 (mod n) and a^(560) ≡ 1 (mod n)
        long[] x = new long[t + 1];
        int i;
        
        // compute x0 = a^u mod n
        x[0] = NumberTheory.modularExponentiation(a, u, n);
        
        System.out.printf("x[%d]=%d\n", 0, x[0]);
        
        boolean c1 = false;
        boolean c2 = false;
        for (i = 1; i <= t; ++i) {
            x[i] = Math.floorMod(x[i - 1]*x[i - 1], n);
            System.out.printf("x[%d]=%d\n", i, x[i]);
            if (x[i] == 1 && x[i - 1] != 1 && x[i - 1] != (n-1)) {
                // x[i-1] is a nontrivial square root of 1, modulo n.
                // n is composite.
                c1 = true;
                break;
            }
        }
        if (x[t] != 1){
            //x_t ≢ (a^(n-1)) (mod n) != 1
            c2 = true;
        }
        assertTrue(Arrays.equals(new long[]{241, 298, 166, 67, 1}, x));
        assertTrue(c1);
        assertFalse(c2);
    }
    
    public void estWitnessAndMillerRabin() throws Exception {
        
        ThreadLocalRandom rand = ThreadLocalRandom.current();
                
        // carmichael number 561 = 3*11*17
        // carmichael number 41041 = 7*11*13*41
        // carmichael number 62745 = 3*5*47*89
        // carmichael number 825265 = 5*7*17*19*73
        
        long a = 7;  // [1, n-1]
        int s = 10;
        
        // test carmichael number as they are composite, not prime
        long[] car = new long[]{561, 41041, 62745, 825265};
        for (long n : car) {
            
            a = rand.nextLong(1, n-1);
            
            // witness might fail invoked singly, but has low probability to
            // fail when invoked several times in the Miller-Rabin probablyPrime
            assertTrue(Primes.witness(a, n, rand) || !Primes.probablyPrime(n, s));
          
            System.out.printf("pollardRhoFactorization(%d)=%s\n",
               n, Arrays.toString(Primes.pollardRhoFactorization(n).toArray()));
        } 
        
    }
    
    public void estNaivePrimeGenerator() {
        
        long p;
        for (int bitLength = 3; bitLength < 5; ++bitLength) {
            p = Primes.naivePrimeGenerator(bitLength);
            System.out.flush();
        }
    }
}
