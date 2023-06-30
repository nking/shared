package algorithms.misc;

import algorithms.VeryLongBitString;
import gnu.trove.iterator.TLongIterator;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TLongArrayList;
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
     * find the prime factors of the integer factorization of n using
     * the Pollard-Rho algorithm and repeating on factors until they're
     * prime.
     * 
     @param n
     @return
     * @throws NoSuchAlgorithmException 
     */
    public static TLongSet findPrimeFactors(final long n) 
        throws NoSuchAlgorithmException {
        
        ThreadLocalRandom rand = ThreadLocalRandom.current();
        
        return findPrimeFactors(n, rand);
    }
    
    /**
     * find prime factors of the integer factorization of n using
     * the Pollard-Rho algorithm and repeating on factors until they're
     * prime.  it should return at least 1 prime.
     * 
     @param n
     @param rand
     @return returns at least one prime factor of n
     * @throws NoSuchAlgorithmException 
     */
    public static TLongSet findPrimeFactors(final long n, ThreadLocalRandom rand) 
        throws NoSuchAlgorithmException {
        
        int s = 10;
        
        TLongSet factors = pollardRhoFactorization(n, rand);
        TLongList factorsList = new TLongArrayList(factors);
        
        TLongSet primeFactors = new TLongHashSet();
        
        TLongSet factors2;
        
        long factor;
        
        long factor2;
        
        TLongIterator iter;
        
        int i = 0;
        
        while (i < factorsList.size()) {
            factor = factorsList.get(i);
            ++i;
            if (probablyPrime(factor, s)) {
                primeFactors.add(factor);
                continue;
            }
                        
            factors2 = pollardRhoFactorization(factor, rand);
            iter = factors2.iterator();
            while (iter.hasNext()) {
                factor2 = iter.next();
                if (probablyPrime(factor2, s)) {
                    primeFactors.add(factor2);
                    continue;
                }
                if (!factors.contains(factor2)) {
                    factorsList.add(factor2);
                }
            }
        }
                
        return primeFactors;
    }
    
    /**
     * integer factorization into primes following Pollard-Rho algorithm in Cormen, Leiserson, Rivest, and Stein 
     * "Introduction to Algorithms".  Usually can find at least one small integer
     * that divides the number n.  The runtime is usually O(n^1/4).
     * "the algorithms is only a heuristic, neither its running time nor its success 
     * is guaranteed, although the procedure is highly effective in practice."
     * 
     * NOTE: for factoring large numbers, may want to implement:
     * "Factoring integers with the number field sieve"
     J. P. Buhler, H. W. Lenstra, Jr., Carl Pomerance
     * http://www.math.leidenuniv.nl/~hwl/PUBLICATIONS/1993e/art.pdf
     * 
     * NOTE: the method returns numbers that are factors, but not prime also.
     * One can use probablyPrime to test for primality and use
     * pollardRhoFactorization() on the non-prime factor.
     * 
     @param n
     @return 
     * @throws java.security.NoSuchAlgorithmException 
     */
    public static TLongSet pollardRhoFactorization(final long n) throws NoSuchAlgorithmException {
                    
        ThreadLocalRandom rand = ThreadLocalRandom.current();
        
        return pollardRhoFactorization(n, rand);
    }
    
    /**
     * integer factorization into primes following Pollard-Rho algorithm in Cormen, Leiserson, Rivest, and Stein 
     * "Introduction to Algorithms".  Usually can find at least one small integer
     * that divides the number n.  The runtime is usually O(n^1/4).
     * "the algorithm is only a heuristic, neither its running time nor its success
     * is guaranteed, although the procedure is highly effective in practice."
     * 
     * NOTE: for factoring large numbers, may want to implement:
     * "Factoring integers with the number field sieve"
      J. P. Buhler, H. W. Lenstra, Jr., Carl Pomerance
     * http://www.math.leidenuniv.nl/~hwl/PUBLICATIONS/1993e/art.pdf
     * 
     * NOTE: the method returns numbers that are factors, but not prime also.
     * One can use probablyPrime to test for primality and use
     * pollardRhoFactorization() on the non-prime factor.

     @param n number to decompose into prime factors
     @param rand rand number generator
     @return set of primes that are factors of n
     */
    public static TLongSet pollardRhoFactorization(final long n, ThreadLocalRandom rand) {
                            
        long i = 1;
        
        long x1 = rand.nextLong(2, n - 1);
        
        long x2;
        
        long y = x1;
        
        long k = 2;
        
        long xCycle = -1;
        
        // first prime expected to be found within n^(1/4) steps:
        long maxI = (long)Math.ceil(Math.pow(n, 0.25));
        
        TLongSet factors = new TLongHashSet(); 
        
        long d;
                
        //System.out.printf("i=%d xi=%d\n", i, x1);
                
        while (true) {
            
            ++i;
            
            x2 = Math.floorMod(x1*x1 - 1, n);
            
            if (x2 == xCycle) {
                break;
            }
            
            d = NumberTheory.extendedEuclid(y - x2, n)[0];
            d = Math.abs(d);
            
            //System.out.printf("i=%d x2=%d, d=(%d,%d) xc=%d\n", 
            //    i, x2, d, NumberTheory.euclid(y-x2, n), xCycle);
                   
            if ((d != 1) && (d != n) && !factors.contains(d)) {
                factors.add(d);
                if (xCycle == -1) {
                    xCycle = x1;
                }
            }
            
            if (i == k) {
                y = x2;
                k *= 2;
            }
            
            x1 = x2;
                        
            if (i > maxI*maxI) {
                break;
            }
        };
        
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
     * The probability that the Miller-Rabin algorithm errs in the result
     * possibly prime is at most 2^(-s) where is any odd integer greater than 2.
     * The return false indicates n is a composite and the result is definite.
     * e.g. 2^(-10) = 1E-3; 2^(-100) ~ 1E-30.
     * The probability that n is prime, given that MILLER-RABIN has returned PRIME
     * is Pr{A|B} which is the alternate form of Bayes’s theorem (equation (C.18)).
     * <pre>
     * Pr{A|B} ≈ 1/(1 + ((2^(-s)) * (ln(n) - 1))).
     * Pr{A|B} does not exceed 1/2 until s .gt. log_2(log(n)-1).
     * so choose s ≥ log_2(log(n)-1) = math.log(beta/1.443)/math.log(2)
     * where beta is the bitlength of n.
     * 
     * reference: Chap 31 of Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms.
     * </pre>
     @param n number to test for primality, must be odd and .gt. 2.
     @param s number of randomly chosen base numbers to try
     @return 
     */
    public static boolean probablyPrime(final long n, final int s) {
        
        if (n < 3) {
            throw new IllegalArgumentException("n must be > 2");
        }
        
        if ((n & 1) == 0) {
            // even numbers are composites of 2
            return false;
        }
        
        if (n < 4) {
            return true;
        }
        
        ThreadLocalRandom rand = ThreadLocalRandom.current();
        
        long a;
        int j;
        
        for (j = 0; j < s; ++j) {
            // 1 <= a <= (n-1)
            a = rand.nextLong(1, n-1);
            if (witness(a, n, rand)) {
                // a is definitely composite
                return false;
            }
        }
        return true;
    }

    /*** implements algorithm from Chap 31 of Cormen, Leiserson, Rivest, and Stein Introduction to Algorithms.
     * tests whether number n is composite (hence not prime) using the
     * number a or possibly prime.
     * witness is a more effective extension of the test 
     * <pre>
     * a^(n-1) ≢ 1 (mod n)
     * </pre>
     @param a number in the range [1, n-1] inclusive, that is a random number 
     * which may prove that n is a composite number, and hence not prime
     @param n the number being tested for primality.  must be odd and .gt. 2.
     @param rand
     @return true when n is composite, else false when n is possibly prime.
     */
    static boolean witness(long a, long n, ThreadLocalRandom rand) {
        
        if ((n & 1) == 0 || n < 3) {
            throw new IllegalArgumentException("n must be odd and > 2");
        }
               
        //n - 1 = (2^t)*u where t >=1 and u is an odd integer
        // n-1/u is an integer 
        // ((n-1)/u) = 2^t
        // t = math.log( (n-1)/u )/math.log(2)

        long u = 1;
        int t = 1;
        int div;
                        
        for (long i = 3; i < n; i+=2) {
            if (Math.floorMod(n - 1, i) == 0) {
                u = i;
                div = (int)((n-1)/u);
                if (MiscMath0.isAPowerOf2(div)) {
                    t = (int)(Math.log( (n-1)/u )/Math.log(2));
                    assert(t > 0);
                    break;
                }
            }
        }
        
        if (u == 1 && t == 1) {
            // keep u=1
            t = (int)(Math.log(n-1)/Math.log(2));
        }
        //System.out.printf("(1<<t)*u=%d, (n-1)=%d, n=%d t=%d u=%d\n", (1<<t)*u, n-1, n, t, u);
        
        assert( (1<<t)*u == (n-1));
                
        long[] x = new long[t + 1];
        int i;
        
        // X = <a^u, a^(2u), a^(2^(2u)),...a^(2^(tu))> (all computations are performed modulo n).
        
        // compute x0 = a^u mod n
        x[0] = NumberTheory.modularExponentiation(a, u, n);
        
        for (i = 1; i <= t; ++i) {
            x[i] = Math.floorMod(x[i - 1]*x[i - 1], n);
            if (x[i] == 1 && x[i - 1] != 1 && x[i - 1] != (n-1)) {
                // x[i-1] is a nontrivial square root of 1, modulo n.
                // n is composite.
                return true;
            }
        }
        if (x[t] != 1){
            //x_t ≢ (a^(n-1)) (mod n) != 1
            return true;
        }
        return false;
    }
    
    /**
     *
     @param bitlength
     @return
     */
    public static long naivePrimeGenerator(int bitlength) {
        ThreadLocalRandom rand = ThreadLocalRandom.current();
        return naivePrimeGenerator(bitlength, rand);
    }
    
    /**
     * generate a prime, naively.
     * <pre>
     * Reference:
     * Joye, Paillier, and Vaudenay "Efficient Generation of Prime Numbers"
     * </pre>
     @param bitLength
     @param rand
     @return 
     */
    public static long naivePrimeGenerator(int bitLength, ThreadLocalRandom rand) {
        if (bitLength < 2) {
            throw new IllegalStateException("bitLength must be > 1");
        }
        // generate a random n-bit odd number and test it for primality and return
        //    when true.
        //    the expected number of iterations is bitLength * log_2(2)/2 = 0.347*bitLength
        // it should find a 32 bit number in 11 iterations.
        int s = 10; // can increase to 100 for very high accuracy
        int maxIter = (int) Math.ceil(0.347 * bitLength);
        int nIter = 0;
        long q;
        do {
            q = generateRandomOdd(bitLength, rand);
            nIter++;
        } while ( !probablyPrime(q, s));
        
        System.out.printf("nIter=%d, max expected=%d prime=%d\n", nIter, maxIter, q);
        
        return q;
    }
    
    static long generateRandomOdd(int bitlength, ThreadLocalRandom rand) {
        // for bit length n=3:  [(1<<(n-1)), (1<<n)-1 ] 4 thru 7, '100' thru '111'
        long q = rand.nextLong((1<<(bitlength - 1)), (1 << bitlength) - 1) ;
        System.out.printf("bitlength=%d rand=%d\n", bitlength, q);
        q <<= 1L;
        q |= 1;
        return q;
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
    
    /**
     * estimate the number of primes less than or equal to n.
     * uses Theorem 31.37 (Prime number theorem) of Cormen, Leiserson, Rivest, and Stein
     * Introduction to Algorithms.
     * The method is off by less than 6% at n = 1E9.
     @param n
     @return 
     */
    public static int numberOfPrimes(int n) {
        return (int)Math.floor((double)n/Math.log(n));
    }
    
     
    /**
     * assuming that the platform word size is either 32 bit or 64 bit, return the
     * largest prime less than the word size
     @return 
     */
    public static long getLargestPrimeForPlatformWordSize() {
        // see http://en.wikipedia.org/wiki/Mersenne_prime
        // for 32 bit  2147483647 which is a Mersenne prime
        // for 62 bit  2305843009213693951
        //     64 bits 9223372036854775783
        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;
        
        if (is32Bit) {
            return 2147483647l;
        } else {
            return 9223372036854775783l;
        }
    }

    /**
     * return a bit vector with set bits for the primes between 2 and n.
     * The "Sieve of Eratosthenes" is used.
     * <pre>
     *     reference:  https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
     * </pre>
     * @param n
     * @return
     */
    public static int[] allPrimesLessThanN(int n) {
        if (n < 2) {
            throw new UnsupportedOperationException("n must be > 1");
        } else if (n == 2) {
            return new int[]{};
        } else if (n == 3) {
            return new int[]{3};
        }
        VeryLongBitString b = new VeryLongBitString(n);
        int i;
        for (i = 2; i < n; ++i) {
            b.setBit(i);
        }
        int j;
        for (i = 2; i < (int)Math.sqrt(n); ++i) {
            if (b.isSet(i)) {
                for (j = i*i; j < n && j <= (Integer.MAX_VALUE - i); j+=i) {
                    // remove all multiples of i
                    b.clearBit(j);
                }
            }
        }
        b.clearBit(2);

        return b.getSetBits();
    }
    
}
