package algorithms.compGeometry;

import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.compGeometry.convexHull.GrahamScan.CHL;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.misc.MiscMath0;
import algorithms.util.PairInt;

import java.util.Set;

/**
 * finds the pair of points most distant from one another a.k.a. the
 * furthest pair problem.
 * one can imagine a rough estimate by fitting an ellipse around the given points
 * and determining the major axis of the points.
 * even more precisely, one can fit a convex hull around the points and
 * determine which pair of points has the largest distance between them.
 * Note that the convex hull returns points ordered in a counter-clockwise manner
 * so is suitable for ordered traversal already.
 * The runtime complexity for this convex hull, which uses a linear time sort,
 * is O(n).
 *
 * One can speculate that a "rotating calipers" approach could reduce the
 * further search by half, that is, only traversing half of the convex hull points
 * (the savings by roughly half arises by terminating each leading points' search
 * for distance pair as soon as the distances start to decrease),
 * but random input testing shows that the possibly very uneven distribution
 * of points on the convex hull means that all points should be visited.
 * 
 * In this algorithm, since integers are given as input, the convex hull
 * uses an internal sorting by polar angles rounded to integers between 0
 * and 359, inclusive, making the convex hull algorithm runtime complexity
 * O(N) where N is the number of x, y points given to this class.
 * If higher resolution than 1 degree is needed for a specific hull, then
 * this algorithm for finding the largest distance pair of points is not
 * the right algorithm for that problem.
 * 
 * The largest distance between a pair algorithm then proceeds to use the number of
 * points on the convex hull, n_H.  The runtime complexity of this part of
 * the algorithm is O(n_h^2).
 * 
 * We can speculate about the size of the hull from information in Cormen
 * et al. Introduction to Algorithms, Exercise 33-5
 * for sparse-hulled distributions of a unit-radius disk, a convex polygon with k sides,
 * and a 2-D normal distribution respectively as n^(1/3), log_2(n), sqrt(log_2(n)).
 * 
 * Then, assuming that N > n_h^2, this algorithm's total runtime complexity is
 * O(N).

first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 * 
 * @author nichole
 */
public class LargestDistancePair {

    /**
     * find the largest distance between the given pairs of points.
     * The runtime complexity is <em>max(O(N), O(n_h^2))</em>.
     *
     * @param x input array of x coordinates.  note that this method modifies the order of x
     *          upon CCW sorting of x and y, so copy the arrays in if need to kep the
     *          original order.
     * @param y x input array of x coordinates.  note that this method modifies the order of x
     *          upon CCW sorting of x and y, so copy the arrays in if need to kep the
     *          original order.
     * @return an array of the pair of points furthest from one another.  the
     * format in the return array is [xa, ya, xb, yb].
     * If there are more than one points with same maximum distance,
     * only one is returned, and that one is the latest point in roughly the 1st half of the
     * convex hull and it last pairing in terms of CCW ordering of the hull points.
     * @throws algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException if the unique polar angles w.r.t. the smallest y point in the arrays
     *                                                                            are fewer than 3 in number.
     */
    public static long[] find(long[] x, long[] y) throws GrahamScanTooFewPointsException {

        // x and y
        CHL ch = GrahamScan.computeHull(x, y, true);

        int n = ch.getXH().length;

        double meanX = MiscMath0.getAvgAndStDev(ch.getXH(), n)[0];
        double meanY = MiscMath0.getAvgAndStDev(ch.getYH(), n)[0];

        long maxDist = Long.MIN_VALUE;
        int iMaxDist = 0;
        long distSq;
        long xd;
        long yd;

        int iRef = getPointFurthestFromMean(meanX, meanY, ch.getXH(), ch.getYH());
        long xRef = ch.getXH()[iRef];
        long yRef = ch.getYH()[iRef];

        //TODO: simplify.  this still has the debugging code in it.

        // find furthest point from the reference point.  that will be the scan range;
        int i = iRef + 1;
        int iter = 0;
        // scan to the 2nd to last point because the last point in the hull is the same as the first point
        while (iter < (n - 1)) {
            if (i > (n - 2)) {
                i = 0;
            }
            xd = ch.getXH()[i] - xRef;
            yd = ch.getYH()[i] - yRef;
            distSq = xd * xd + yd * yd;

            if (distSq >= maxDist) {
                maxDist = distSq;
                iMaxDist = i;
            } else {
                // because the hull points are sorted in CCW order, the max dist
                // will increase or stay the same and then decrease.
                break;
            }
            ++iter;
            ++i;
        }

        // using a larger scan range than iter is also necessary
        int nScan = n - 1;//iter;

        //System.out.printf("nScan = %d, n-1=%d\n", iter, n - 1);

        //System.out.printf("i0=%d, i1=%d, distSq=%d\n", iRef, iMaxDist, maxDist);

        int j;

        long maxDistJ;
        int jMaxDist = 0;

        // this O(n_hull^2) section is necessary

        // try pairs i0: [2, nScan) to see if an i1 pairing has dist > maxDist
        i = iRef + 1;
        iter = 0;
        int iterJ;
        while (iter <= nScan) {
            if (i > (n - 2)) {
                i = 0;
            }
            xRef = ch.getXH()[i];
            yRef = ch.getYH()[i];
            //System.out.printf("%d)\n", i);

            iterJ = 0;
            j = i + 1;

            maxDistJ = Long.MIN_VALUE;
            while (iterJ <= nScan) {
                if (j > (n - 2)) {
                    j = 0;
                }
                xd = ch.getXH()[j] - xRef;
                yd = ch.getYH()[j] - yRef;
                distSq = xd * xd + yd * yd;
                //System.out.printf("   %10d  %20d", j, distSq);
                if (distSq >= maxDistJ) {
                    maxDistJ = distSq;
                    jMaxDist = j;
                    //System.out.printf("*\n");
                } else {
                    //System.out.printf("\n");
                    break;
                }
                ++iterJ;
                ++j;
            } // end loop over j
            if (maxDistJ > maxDist) {
                iRef = i;
                iMaxDist = jMaxDist;
                maxDist = maxDistJ;
                //System.out.printf("    j: i0=%d, i1=%d, distSq=%d\n", iRef, iMaxDist, maxDist);
            }
            ++iter;
            ++i;
        } //end loop over i
        
        /*System.out.printf("MAXDISTSQ=%d  (x[%d]=%d, y[%d]=%d), (x[%d]=%d, y[%d]=%d)\n", 
            maxDist, 
            iRef, ch.getXH()[iRef], 
            iRef, ch.getYH()[iRef],
            iMaxDist, ch.getXH()[iMaxDist],
            iMaxDist, ch.getYH()[iMaxDist]);*/

        long[] result = new long[]{ch.getXH()[iRef], ch.getYH()[iRef], ch.getXH()[iMaxDist], ch.getYH()[iMaxDist]};

        assert (maxDist > Long.MIN_VALUE);
        assert (maxDist == ((result[0] - result[2]) * (result[0] - result[2]) +
                (result[1] - result[3]) * (result[1] - result[3])));

        //TODO: edit to return the entire hull and the furthest pair indexes

        return result;
    }

    protected static int getPointFurthestFromMean(double meanX, double meanY,
                                                  long[] xh, long[] yh) {

        int i;
        double xd;
        double yd;
        double distSq;
        double maxDist = Double.NEGATIVE_INFINITY;
        int iMaxDist = -1;
        for (i = 0; i < xh.length; ++i) {
            xd = xh[i] - meanX;
            yd = yh[i] - meanY;
            distSq = xd * xd + yd * yd;
            if (distSq > maxDist) {
                maxDist = distSq;
                iMaxDist = i;
            }
        }
        return iMaxDist;
    }

    public static long[] find(Set<PairInt> points) throws GrahamScanTooFewPointsException {
        long[] x = new long[points.size()];
        long[] y = new long[points.size()];
        int i = 0;
        for (PairInt p : points) {
            x[i] = p.getX();
            y[i] = p.getY();
            ++i;
        }
        return find(x, y);
    }

}
