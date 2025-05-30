package algorithms.misc;

/**
 * implemented from pseudocode from Cormen, Leiserson, Rivest, and Stein 
 * "Introduction to Algorithms", Chap 31
 *
 * first implemented in project
     http://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)
   then moved here to share with other projects.
   
 * @author nichole
 */
public class NumberTheory {

    /**
     * return the greatest common denominator of the 2 integers.
     *
     * runtime complexity is less than the number of base10 digits in a.
     * r.t.c. O(nlogn)
     * 
     @param a
     @param b
     @return
     */
    public static int euclid(int a, int b) {
        int t;
        while (b != 0) {
            t = b;
            b = Math.floorMod(a, b);//a % b;  a - floorDiv(a, b) * b
            a = t;
        }
        return Math.max(a, -a);
    }

    /**
     * return the greatest common denominator of the 2 integers.
     *
     * runtime complexity is O(nlogn)
     * 
     @param a
     @param b
     @return
     */
    public static long euclid(long a, long b) {
        long swap;
        while (b != 0) {
            swap = b;
            b = Math.floorMod(a, b);
            a = swap;
        }
        return Math.max(a, -a);
    }
    
    /**
     * solves for x in the equation a * x ≡ b (mod n) (which is actually (a*x) % n = b)
     * where d is the gcd of number n and d|b (a.k.a. d divides b).
     * finds the smallest gcd for which a*x + b*y = d.
     * The equation may have zero, one, or more than one such solution.
     * performs O(lg n + gcd(a, n)) arithmetic operations.
     * <pre>
     * Section 31.4, Cormen, Leiserson, Rivest, and Stein Introduction to Computer Algorithms.
     * </pre>
     @param a positive number greater than 0
     @param b positive number greater than 0
     @param n
     @return 
     */
    public static long[] gcdModularLinearEqnSolver(final long a, final long b, final long n) {
        
        /*
        https://en.wikipedia.org/wiki/B%C3%A9zout%27s_identity
        Bézout's identity (also called Bézout's lemma) is a theorem in 
        elementary number theory: let a and b be nonzero integers and let d be 
        their greatest common divisor. Then there exist integers x and y such 
        that 
            ax+by=d.  
        In addition, the greatest common divisor d is the 
        smallest positive integer that can be written as ax + by every integer 
        of the form ax + by is a multiple of the greatest common divisor d.
        The integers x and y are called Bézout coefficients for (a, b); they 
        are not unique. A pair of Bézout coefficients can be computed by the 
        extended Euclidean algorithm.
        */
        
        long min = Long.MAX_VALUE;
        
        long[] dXY = extendedEuclid(a, n);
        long d = dXY[0];
        if (d == 0L || ((b % d) != 0)) {
            return new long[]{min};
        }
        long m = dXY[1] * (b / d);
        // use the floor modulo operator instead of the default truncated which is '%'
        long x0 = Math.floorMod(m, n);
        long[] s = new long[(int)d];
        
        for (int i = 0; i < d; ++i) {
            s[i] = (x0 + i * (n / d)) % n;
        }
        return s;
    }
    
    /**
     * calculate a^b mod n.
     * <pre>
     * Chap 31, MODULAR-EXPONENTIATION(a, b, n) from Cormen, Leiserson, Rivest, and Stein Introduction
     * to Algorithms (a.k.a. CLRS).
     * </pre>
     @param a non-negative integer
     @param b non-negative integer
     @param n positive integer
     @return 
     */
    public static int modularExponentiation(int a, int b, int n) {
        
        if (a == 0) {
            return 0;
        }
        if (b == 0 && a < n) {
            return a;
        }
        int c = 0;
        int d = 1;
        int k = MiscMath0.numberOfBits(b);
        int mask = 1 << (k - 1);
        for (int i = k - 1; i >= 0; --i) {
            c *= 2;
            d = Math.floorMod(d*d, n);//(d*d) % n;
            if ((b & mask) != 0) {
                c++;
                d = Math.floorMod(d*a, n);//(d*a) % n;
            } 
            mask >>= 1;
        }
        return d;
    }
    
    /**
     * calculate a^b mod n.
     * <pre>
     * Chap 31, MODULAR-EXPONENTIATION(a, b, n) from Cormen, Leiserson, Rivest, and Stein Introduction
     * to Algorithms (a.k.a. CLRS).
     * </pre>
     @param a non-negative integer
     @param b non-negative integer
     @param n positive integer
     @return 
     */
    public static long modularExponentiation(long a, long b, long n) {
        
        if (a == 0) {
            return 0;
        }
        if (b == 0 && a < n) {
            return a;
        }
        int c = 0;
        long d = 1;
        int k = MiscMath0.numberOfBits(b);
        long mask = 1 << (k - 1);
        for (int i = k - 1; i >= 0; --i) {
            c *= 2;
            d = Math.floorMod(d*d, n);//(d*d) % n;
            if ((b & mask) != 0) {
                c++;
                d = Math.floorMod(d*a, n);//(d*a) % n;
            } 
            mask >>= 1;
        }
        return d;
    }
    
    /*
    r.t. complexity of multiplying 2 n bit numbers is
         O(n log n log log n).
    */

    /**
     * 
     * extended euclid returns d, x, and y where
     * d = gcd(a, b) = a*x + b*y where x and y may be zero or negative.
     * x and y are useful for forming multiplicative inverses.
     * 
     * if a .gt. b .geq. 0, runtime complexity is O(log_2(b)).
     * 
     @param a
     @param b
     @return returns d, x, and y where
     * d = gcd(a, b) = a*x + b*y where x and y may be zero or negative.
     * x and y are useful for forming multiplicative inverses.
     * 
     */
    public static long[] extendedEuclid(long a, long b) {
        if (b == 0) {
            return new long[]{a, 1, 0};
        }
        
        long[] dxyP = extendedEuclid(b, a % b);
        
        long t = (long)Math.floor((double)a/(double)b);
        long r = dxyP[1] - t*dxyP[2];
        
        return new long[] {dxyP[0], dxyP[2], r};
    }
    
    /**
     * the least common multiple, lowest common multiple, or smallest common 
     * multiple of two integers a and b, usually denoted by lcm(a, b), 
     * is the smallest positive integer that is divisible by both a and b.
     *
     * e.g. for a=12, b=18, result= 12/2 * 18/2 = 6*6 = 36
     * <pre>
     * reference:
     * https://en.m.wikipedia.org/wiki/Least_common_multiple
     * </pre>
     @param a
     @param b
     @return 
     */
    public static long leastCommonMultiple(long a, long b) {
        if (a == 0 && b == 0) {
            throw new IllegalArgumentException("a and b cannot both be0");
        }
        
        long m = Math.abs(a * b);
        long gcd = euclid(a, b);
        
        return m/gcd;
    }

    /**
     * gererate a Pythagorean triple of numbers, that is a, b and c where
     * (a*a) + (b*b) = (c*c) using Euclid's formula.
     * @param m
     * @param n
     * @return
     */
    public static long[] generatePythagoreanTriples(int n, int m) {
        if (m <= 0 || n <= 0 || m < n) {
            throw new IllegalArgumentException("m, n must be > 0, and m > n.");
        }
        return generatePythagoreanTriples(n, m, 1);
    }
    public static long[] generatePythagoreanTriples(int n, int m, int k) {
        if (m <= 0 || n <= 0 || m < n || k<=0) {
            throw new IllegalArgumentException("m, n, k must be > 0, and m > n.");
        }
        long m2 = m*m;
        long n2 = n*n;
        return new long[]{k * (m2 - n2), k * (2 * m * n), k * (m2 + n2)};
    }
}
