package algorithms.misc;

/**
 * implemented from pseudocode from Cormen et al. 
 * "Introduction to Algorithms", Chap 31
 *
 * first implemented in project
     http://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)
   then moved here to share with other projects.
   
 * @author nichole
 */
public class GreatestCommonDenominator {

    /**
     * return the greatest common denominator of the 2 integers.
     *
     * runtime complexity is (?)
     * 
     * @param a
     * @param b
     * @return
     */
    public static int euclid(int a, int b) {
        if (b == 0) {
            return a;
        }
        count++;
        return euclid(b, a % b);
    }
    public static int count = 0;

    /**
     * return the greatest common denominator of the 2 integers.
     *
     * runtime complexity is (?)
     * 
     * @param a
     * @param b
     * @return
     */
    public static long euclid(long a, long b) {
        //System.out.format("a=%d b=%d\n", a, b); System.out.flush();
        if (b == 0) {
            //System.out.format("   euclid=%d\n", a); System.out.flush();
            return a;
        }
        count++;
        return euclid(b, a % b);
    }
    
    /**
     * solves the equation a * x = b mod n to
     * find the smallest gcd for which a*x + b*y = d where d is a
     * gcd of number n.
     * @param a
     * @param b
     * @param n
     * @return 
     */
    public static long gcdModularLinearEqnSolver(long a, long b, long n) {
        
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
        
        long[] d_xp_yp = extendedEuclid(a, n);
        if ((d_xp_yp[0] != 0) || d_xp_yp[2] != 0) {
            long d = d_xp_yp[0];
            long x0 = d_xp_yp[1] * (b/d) % n;
            for (int i = 0; i < d; ++i) {
                
                long x1 = (x0 + i*(n/d)) % n;
                
                //System.out.println(" " + d);
                
                if (d > 0 && d < min) {
                    min = d;
                }
            }
        }
        return min;
    }
    
    public static int modularExponentiation(int a, int b, int n) {
        int c = 0;
        int d = 1;
        int nbits = MiscMath0.numberOfBits(b);
        for (int i = nbits - 1; i >= 0; --i) {
            c *= 2;
            d = (d*d) % n;
            if ((b & (1 << i)) != 0) {
                c++;
                d = (d*a) % n;
            }
        }
        return d;
    }
    
    /*
    r.t. complexity of multiplying 2 n bit numbers is
         O(n log n log log n).
    */

    /**
     * 
     * extended euclid 
     * 
     * @param a
     * @param b
     * @return 
     */
    public static long[] extendedEuclid(long a, long b) {
        if (b == 0) {
            return new long[]{a, 1, 0};
        }
        
        long[] dxy_p = extendedEuclid(b, a % b);
        
        long t = (long)Math.floor((double)a/(double)b);
        
        //System.out.format("a=%d b=%d (a/b)=%d d=%d x=%d y=%d\n", 
        //    a, b, t, dxy_p[0], dxy_p[2], dxy_p[1] - t*dxy_p[2]);
        
        long[] dxy = new long[] {
            dxy_p[0], dxy_p[2], (dxy_p[1] - t*dxy_p[2])
        };
        return dxy;
    }
    
}
