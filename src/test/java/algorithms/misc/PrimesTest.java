
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
        
        long n;
        TLongSet expected;
        
        n = 1387;
        expected = new TLongHashSet();
        expected.add(Long.valueOf(19));
        expected.add(Long.valueOf(73));
        
        assertTrue(checkPollardRhoFactorization(n, expected));
       
        //----
        n = 825;
        expected = new TLongHashSet();
        //  prime answer using LCDs = 11, 5, 3
        expected.add(Long.valueOf(11));
        expected.add(Long.valueOf(75));
        expected.add(Long.valueOf(5));
        expected.add(Long.valueOf(3));
        
        assertTrue(checkPollardRhoFactorization2(n, expected));
        
        //---
        n = 197;
        expected = new TLongHashSet();
        //expected.add(Long.valueOf(197));
        
        assertTrue(checkPollardRhoFactorization(n, expected));
       
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

    /**
     * check that all in expected are found
     * @param n
     * @param expected
     * @return
     * @throws NoSuchAlgorithmException 
     */
    private boolean checkPollardRhoFactorization(long n, TLongSet expected) throws NoSuchAlgorithmException {
        TLongSet result;
        TLongIterator iter;
        long r;
        
    //TODO: revisit the code and fix this test    
        
        int count = 0;
        int nTries = 10;
        
        while (count < nTries) {
            //System.out.printf("===> try iteration %d <===\n", count);
            count++;
            result = Primes.pollardRhoFactorization(n);
            if (expected.size() != result.size()) {
                continue;
            }
            TLongSet copy = new TLongHashSet();
            copy.addAll(expected);
            iter = expected.iterator();
            while (iter.hasNext()) {
                r = iter.next();
                if (result.contains(r)) {
                    copy.remove(r);
                }
            }
            if (copy.isEmpty()) {
                //System.out.println("nTries=" +count);
                return true;
            }
        }
        return false;
    }
    
    /**
     * check that any in expected are present
     * @param n
     * @param expected
     * @return
     * @throws NoSuchAlgorithmException 
     */
    private boolean checkPollardRhoFactorization2(long n, TLongSet expected) throws NoSuchAlgorithmException {
        TLongSet result;
        TLongIterator iter;
        long r;
        
        int count = 0;
        int nTries = 10;
        boolean found;
        
        while (count < nTries) {
            //System.out.printf("===> try iteration %d <===\n", count);
            count++;
            result = Primes.pollardRhoFactorization(n);
            if (result.isEmpty() && !expected.isEmpty()) {
                continue;
            }
            iter = result.iterator();
            while (iter.hasNext()) {
                r = iter.next();
                if (expected.contains(r)) {
                    //System.out.println("nTries=" +count);
                    return true;
                }
            }
        }
        return false;
    }
    
    public void testWitnessAndMillerRabin() throws Exception {
        
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
