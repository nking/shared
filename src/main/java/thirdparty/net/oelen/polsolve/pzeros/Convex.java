package thirdparty.net.oelen.polsolve.pzeros;


/************************************************************************
*                             CLASS CONVEX                              *
*************************************************************************
* Compute  the upper convex hull of the set (i,a(i)), i.e., the set of  *
* vertices (i_k,a(i_k)), k=1,2,...,m, such that the points (i,a(i)) lie *
* below the straight lines passing through two consecutive vertices.    *
* The abscissae of the vertices of the convex hull equal the indices of *
* the TRUE  components of the logical output vector H.                  *
* The used method requires O(nlog n) comparisons and is based on a      *
* divide-and-conquer technique. Once the upper convex hull of two       *
* contiguous sets  (say, {(1,a(1)),(2,a(2)),...,(k,a(k))} and           *
* {(k,a(k)), (k+1,a(k+1)),...,(q,a(q))}) have been computed, then       *
* the upper convex hull of their union is provided by the subroutine    *
* CMERGE. The program starts with sets made up by two consecutive       *
* points, which trivially constitute a convex hull, then obtains sets   *
* of 3,5,9... points,  up to  arrive at the entire set.                 *
* The program uses the subroutine  CMERGE; the subroutine CMERGE uses   *
* the subroutines LEFT, RIGHT and CTEST. The latter tests the convexity *
* of the angle formed by the points (i,a(i)), (j,a(j)), (k,a(k)) in the *
* vertex (j,a(j)) up to within a given tolerance TOLER, where i<j<k.    *
*************************************************************************/

strictfp class Convex {
    
    private final static double TOLER = 0.4;	// slope tolerace
    
    private static int MIN(int x, int y) {
        return (x < y) ? x : y;
    }
    
    
    
    

    /************************************************************************
    *                             SUBROUTINE LEFT                           *
    ************************************************************************* 
    * Given as input the integer I and the vector H of logical, compute the *
    * the maximum integer IL such that IL<I and H(IL) is TRUE.              *
    ************************************************************************* 
    * Input variables:                                                      *
    *     N   : length of the vector H                                      *
    *     H   : vector of logical                                           *
    *     I   : integer                                                     *
    ************************************************************************* 
    * Output variable:                                                      *
    *     IL  : maximum integer such that IL<I, H(IL)=.TRUE.                *
    *************************************************************************/
    private static int left(int i, int lo, boolean[] h) {
        if (i == lo) {
            return lo;
        }
        for (i--; i > lo; i--) {
            if (h[i]) {
                break;
            }
        }
        return i;
    }
    
    
    
    

    /************************************************************************
    *                             SUBROUTINE RIGHT                          *
    *************************************************************************
    ************************************************************************* 
    * Given as input the integer I and the vector H of logical, compute the *
    * the minimum integer IR such that IR>I and H(IL) is TRUE.              *
    *************************************************************************
    ************************************************************************* 
    * Input variables:                                                      *
    *     N   : length of the vector H                                      *
    *     H   : vector of logical                                           *
    *     I   : integer                                                     *
    ************************************************************************* 
    * Output variable:                                                      *
    *     IR  : minimum integer such that IR>I, H(IR)=.TRUE.                *
    *************************************************************************/
    
    private static int right(int i, int up, boolean[] h) {
        if (i == up) {
            return up;
        }
        for (i++; i < up; i++) {
            if (h[i]) {
                break;
            }
        }
        return i;
    }
    
    
    
    

    /************************************************************************
    *                             FUNCTION CTEST                            *
    ************************************************************************* 
    * Test the convexity of the angle formed by (IL,A(IL)), (I,A(I)),       *
    * (IR,A(IR)) at the vertex (I,A(I)), up to within the tolerance         *
    * TOLER. If convexity holds then the function is set to .TRUE.,         *
    * otherwise CTEST=.FALSE. The parameter TOLER is set to 0.4 by default. *
    ************************************************************************* 
    * Input variables:                                                      *
    *     N       : length of the vector A                                  *
    *     A       : vector of double                                        *
    *     IL,I,IR : integers such that IL<I<IR                              *
    *************************************************************************
    * Output:                                                               *
    *     .TRUE. if the angle formed by (IL,A(IL)), (I,A(I)), (IR,A(IR)) at *
    *            the vertex (I,A(I)), is convex up to within the tolerance  *
    *            TOLER, i.e., if                                            *
    *            (A(I)-A(IL))*(IR-I)-(A(IR)-A(I))*(I-IL)>TOLER.             *
    *     .FALSE.,  otherwise.                                              *
    *************************************************************************/
    
    private static boolean ctest(int il, int i, int ir, double[] a) {
        double s1, s2;

        s1 = (a[i] - a[il]) * (ir - i);
        s2 = (a[ir] - a[i]) * (i - il);
        return (s1 - s2 > TOLER);
    }
    
    
    
    

    /************************************************************************
    *                             SUBROUTINE CMERGE                         *
    *************************************************************************
    * Given the upper convex hulls of two consecutive sets of pairs         *
    * (j,A(j)), compute the upper convex hull of their union                *
    *************************************************************************
    * Input variables:                                                      *
    *     N    : length of the vector A                                     *
    *     A    : vector defining the points (j,A(j))                        *
    *     I    : abscissa of the common vertex of the two sets              *
    *     M    : the number of elements of each set is M+1                  *
    *************************************************************************
    * Input/Output variable:                                                *
    *     H    : vector defining the vertices of the convex hull, i.e.,     *
    *            H(j) is .TRUE. if (j,A(j)) is a vertex of the convex hull  *
    *            This vector is used also as output.                        *
    *************************************************************************/
    
    private static void cmerge(int lo, int i, int up, double a[], boolean[] h) {
        int il, ir, ill, irr;
        boolean tstl, tstr;

        ill = lo;
        irr = up;
        il = left(i, lo, h);
        ir = right(i, up, h);
        if (ctest(il, i, ir, a)) {
            return;
        }
        h[i] = false;
        do {
            if (il == lo) {
                tstl = true;
            } 
            else {
                ill = left(il, lo, h);
                tstl = ctest(ill, il, ir, a);
            }

            if (ir == up) {
                tstr = true;
            } 
            else {
                irr = right(ir, up, h);
                tstr = ctest(il, ir, irr, a);
            }

            if (!tstl) {
                h[il] = false;
                il = ill;
            }
            if (!tstr) {
                h[ir] = false;
                ir = irr;
            }
        } while (!(tstl && tstr));
    }

    
    
    
    
    /**************************************************************
    *                       SUBROUTINE CNVEX                      *
    * *************************************************************
    * Compute the convex hull of the data set a[]. The result     *
    * is in the bool vector h[]. The algorithm successively       *
    * merges adjacent convex hulls of sizes 2, 4, 8, ...          *
    ***************************************************************/
    
    static void cnvex(int n, double a[], boolean[] h) {
        int m, c;

        for (m = 0; m <= n; m++) {
            h[m] = true;
        }

        for (m = 1; m < n; m <<= 1) {
            for (c = m; c < n; c += 2 * m) {
                cmerge(c - m, c, MIN(n, c + m), a, h);
            }
        }
    }
}
