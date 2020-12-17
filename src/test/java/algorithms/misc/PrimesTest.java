
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
        
        assertTrue(check(n, expected));
       
        //----
        n = 825;
        expected = new TLongHashSet();
        //  prime answer using LCDs = 11, 5, 3
        expected.add(Long.valueOf(11));
        expected.add(Long.valueOf(75));
        expected.add(Long.valueOf(5));
        expected.add(Long.valueOf(3));
        
        assertTrue(check2(n, expected));
        
        //---
        n = 197;
        expected = new TLongHashSet();
        //expected.add(Long.valueOf(197));
        
        assertTrue(check(n, expected));
       
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
        r = GreatestCommonDenominator.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        b = n; a = 75;
        r = GreatestCommonDenominator.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        a = n; b = 55;
        r = GreatestCommonDenominator.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        b = n; a = 55;
        r = GreatestCommonDenominator.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        a = n; b = 15;
        r = GreatestCommonDenominator.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        b = n; a = 15;
        r = GreatestCommonDenominator.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        a = n; b = 5;
        r = GreatestCommonDenominator.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        b = n; a = 5;
        r = GreatestCommonDenominator.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        a = n; b = 3;
        r = GreatestCommonDenominator.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
        
        b = n; a = 3;
        r = GreatestCommonDenominator.extendedEuclid(a, b);
        System.out.printf("EE(%d,%d)=%s\n", a, b, Arrays.toString(r));
    }

    /**
     * check that all in expected are found
     * @param n
     * @param expected
     * @return
     * @throws NoSuchAlgorithmException 
     */
    private boolean check(long n, TLongSet expected) throws NoSuchAlgorithmException {
        TLongSet result;
        TLongIterator iter;
        long r;
        
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
    private boolean check2(long n, TLongSet expected) throws NoSuchAlgorithmException {
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
    
    public void testRabinMiller() throws Exception {
        
        ThreadLocalRandom rand = ThreadLocalRandom.current();
        int number = 3;
        int nTries = 10;
        boolean ans = Primes.witness(2, number, rand);
        assertTrue(ans);
        
        boolean ans0 = Primes.probablyPrime(number, nTries);
        assertTrue(ans0);
    }
}
