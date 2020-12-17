package algorithms.misc;

import gnu.trove.iterator.TLongIterator;
import gnu.trove.map.TLongObjectMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;

/**
 *
 * @author nichole
 */
public class Primes {
    
    /**
     * integer factorization following Pollard-Rho algorithm in Cormen et al. 
     * "Introduction to Algorithms".  Usually can find at least one small integer
     * that divides the number n.  The runtime is usually O(n^1/4).
     * 
     * NOTE: for factoring large numbers, may want to implement:
     * "Factoring integers with the number field sieve"
        J. P. Buhler,  H. W. Lenstra, Jr., Carl Pomerance
     * http://www.math.leidenuniv.nl/~hwl/PUBLICATIONS/1993e/art.pdf
     * 
     * @param n
     * @return 
     * @throws java.security.NoSuchAlgorithmException 
     */
    public static TLongSet pollardRhoFactorization(final long n) throws NoSuchAlgorithmException {
        
        ThreadLocalRandom rand = ThreadLocalRandom.current();
        
        boolean useEdits = true;
        
        TLongSet factors = new TLongHashSet();
        
        /* only the most recent value of x_i needs to be retained, so to reduce
        space, will only keep x_latest and x_latest_index. */
        
        long i = 1;
        
        long x_latest = rand.nextLong(n - 1);

        long y = x_latest;
        long k = 2;
        long[] dxy;
        
        long maxIter = (long)Math.ceil(Math.pow(n, 0.25));
        long m, r, dabs;
        
        // store y and x_latest pairs and break when cycle returns
        TLongObjectMap<TLongSet> pairs = new TLongObjectHashMap<TLongSet>();
        TLongSet pairsV;
        
        while (true) {
            
            //System.out.printf("  n=%d i=%d) x_latest=%d\n", n, i, x_latest);  System.out.flush();
            
            i++;
            x_latest = ((x_latest * x_latest) - 1) % n;
            
            long tmp = y - x_latest;
            
            // store for cycle checks:
            if (pairs.containsKey(y)) {
                if (pairs.get(y).contains(x_latest)) {
                    break;
                } else {
                    pairs.get(y).add(x_latest);
                }
            } else {
                pairsV = new TLongHashSet();
                pairsV.add(x_latest);
                pairs.put(y, pairsV);
            }
                        
            dxy = GreatestCommonDenominator.extendedEuclid(tmp, n);
            //long d = GreatestCommonDenominator.euclid(tmp, n);
            dabs = Math.abs(dxy[0]);
            
            //System.out.printf("    y=%d, x_latest=%d, EE(%d,%d)=%s\n", 
            //    y, x_latest, tmp, n, Arrays.toString(dxy));  System.out.flush();
            
            if ((dabs != 1) && (dabs != Math.abs(n))) {
                factors.add(dabs);
                //System.out.printf(" * store %d (size=%d)\n", dabs, factors.size());
            }
            
            // check for stopping conditions
            m = multiply(factors);
            
            /*if(!factors.isEmpty()) {
                System.out.printf("   %s=> m=%d (n=%d x_latest=%d)\n", 
                    Arrays.toString(factors.toArray()), m, n, x_latest);
                System.out.flush();
            }*/
            
            if (m == Math.abs(n)) {
                break;
            }
            if (useEdits && m > Math.abs(n)) {
                return new TLongHashSet();
                /*while (m > n && !factors.isEmpty()) {
                    long rm = max(factors);
                    factors.remove(rm);
                    m = multiply(factors);
                }*/
            }
            
            if (i == k) {
                y = x_latest;
                k *= 2;
            }            
        }
        //System.out.printf("  i=%d  n^(1/4)=%d\n", i, maxIter);
        //System.out.flush();

        return factors;
    }
    
    /*
    from Aho & Ullman:
    
    Fermat’s theorem states that if p is a prime, and a is any integer between 1 and p − 1, 
    then a^(p−1) leaves a remainder of 1 when divided by p.
    
    usually, if a is chosen at random between 1 and p−1, then the probability
    is at least 1/2 that a^(p−1) will have a remainder other than 1 when divided by p.
    
    e.g., let p = 7. Then 1^6, 2^6,...6^6 are respectively 1, 64, 729, 4096, 15625, and 46656. 
    Their remainders when divided by 7 are all 1. 
    However, if p = 6, a composite number, 
    then 1^5, 2^5,...5^5 are respectively 1, 32, 243, 1024, and 3125,
    Their remainders when divided by 6 are 1, 2, 3, 4, and 5. 
    ==> Only 20% are 1.
    
    Thus, the “algorithm” for testing whether a number p is a prime is to select
    k integers from 1 to p − 1, independently and at random. 
    
    If for any selected a we find the remainder of ( a^(p−1) )/p to be other 
    than 1, we say p is composite; otherwise, we say it is prime. 
    If it weren’t for the “bad” composites, we could say that the
    probability of failure is at most 2^(−k).
    */

    /**
     * using the Miller-Rabin algorithm uses s random tries to determine whether
     * n is definitely a composite number or possibly is prime.
     * implements algorithm from Chap 31 of Cormen et al. Introduction to Algorithms.
     * @param n
     * @param s number of random tries to use
     * @return 
     */
    public static boolean probablyPrime(long n, int s) {
        
        ThreadLocalRandom rand = ThreadLocalRandom.current();
        long a;
        int j;
        
        for (j = 0; j < s; ++j) {
            a = rand.nextLong(1, n-1);
            if (witness(a, n, rand)) {
                return false;
            }
        }
        return true;
    }

    /*** implements algorithm from Chap 31 of Cormen et al. Introduction to Algorithms.
     * @param n
     * @return 
     */
    static boolean witness(long a, long n, ThreadLocalRandom rand) {
                
        // a^(n-1) = (a^u)^(2^t)
        
        /*
        let n -1 = u*2^t 
            where u is odd and t >= 1
        x_0 = modularExponentiation(a, u, n)
        for (i = 1 to t
            x_i = (x_(i-1))^2 mod n
            if (x_i == 1 and x_(i-1) != 1 and x_(i-1) != (n-1)) 
                return true
        if (x_i != 1
            return true
        return false
        */
        
        
        /* let n-1 = 2^t * u   where t>=1 and u is odd
        
        also, t < n
        
        2^t = (n-1)/u
        log_2( 2^t ) = t = log_2( (n-1)/u )
        
        u = (n-1)/(2^t)
        */
        // randomly choose t as 1 <= t <= (n-1), but use a maximum value of integer
        //    if n is a long larger than 31 bits
        int t = 1;
        long u = 1;
        
        do {
            if (n > Integer.MAX_VALUE) {
                t = rand.nextInt(1, Integer.MAX_VALUE - 1);
            } else {
                t = rand.nextInt(1, (int)n - 1);
            }
            u = (n - 1)/((long)Math.pow(2, t));
        } while ((u & 1) == 0);
        
        long[] x = new long[t + 1];
        int i;
        long xim;
        x[0] = MiscMath0.modularExponentiation(a, u, n);
        
        for (i = 1; i < t; ++i) {
            xim = x[i - 1];
            x[i] = (xim*xim) % n;
            if (x[i] == 1 && xim != 1 && xim != (n-1)) {
                return true;
            }
        }
        if (x[t - 1] != 1){
            return true;
        }
        return false;
    }

    private static long multiply(TLongSet factors) {
        TLongIterator iter = factors.iterator();
        long m = 1;
        while (iter.hasNext()) {
            m *= iter.next();
        }
        return m;
    }
    
    private static long max(TLongSet factors) {
        TLongIterator iter = factors.iterator();
        long m = Long.MIN_VALUE;
        long r;
        while (iter.hasNext()) {
            r = iter.next();
            if (r > m) {
                m = r;
            }
        }
        return m;
    }
}
