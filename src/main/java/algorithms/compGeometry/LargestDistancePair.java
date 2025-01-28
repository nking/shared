package algorithms.compGeometry;

import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.compGeometry.convexHull.GrahamScan.CHL;
import algorithms.compGeometry.convexHull.GrahamScanPairInt;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.matrix.MatrixUtil;
import algorithms.util.PairInt;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.TDoubleSet;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 * finds the pair of points most distant from one another a.k.a. the
 * furthest pair problem.
 *
 * Note that there is a maximum manhattan distance between pairs in Distance.maxManhattan(points)
 * that is O(n*8) ~ O(n) for dimensions = 2.
 *
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
     * The runtime complexity is <em>max(O(N), O(n_h^2))</em> if useLinear is true,
     * else is O(N*log_2(N)) if useLinear is false.
     *
     * @param x input array of x coordinates.  note that this method modifies the order of x
     *          upon CCW sorting of x and y, so copy the arrays in if need to kep the
     *          original order.
     * @param y x input array of x coordinates.  note that this method modifies the order of x
     *          upon CCW sorting of x and y, so copy the arrays in if need to kep the
     *          original order.
     * @param useLinear if true, uses a O(n) runtime algorithm to compute the convex hull.  if false, uses
     *                  a O(n*log_2(n)) algorithm to compute the convex hull.
     *                  The linear algorithm is usually correct.   random input testing shows that
     *                  for 3 out of 100 tests, the largest separation pair is wrong, shorter than the true value
     *                  by less than 10 percent.
     *                  TODO: test more to improve the possible error stats empirically, and consider the geometry for
     *                  theoretical limits to error.
     * @return the convex hull and the indexes for the pairs of points on the hull which are furthest from each other.
     *
     * @throws algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException if the unique polar angles w.r.t. the smallest y point in the arrays
     *                                                                            are fewer than 3 in number.
     */
    public static PairAndHull find(long[] x, long[] y, boolean useLinear) throws GrahamScanTooFewPointsException {

        // x and y
        CHL ch = GrahamScan.computeHull(x, y, useLinear);

        int n = ch.getXH().length;

        long maxDistSq = Long.MIN_VALUE;
        int maxI = -1;
        int maxJ = -1;
        long maxDistSqJ;
        int maxJ2 = -1;

        long distSq;
        long xd;
        long yd;

        long xRef;
        long yRef;

        int nScan = n - 1;

        int j;
        int i;
        for (i = 0; i < nScan; ++i) {
            xRef = ch.getXH()[i];
            yRef = ch.getYH()[i];
            maxDistSqJ = Long.MIN_VALUE;

            for (j = i+1; j < nScan; ++j) {
                xd = ch.getXH()[j] - xRef;
                yd = ch.getYH()[j] - yRef;
                distSq = xd * xd + yd * yd;
                if (distSq >= maxDistSqJ) {
                    maxDistSqJ = distSq;
                    maxJ2 = j;
                } else {
                    break;
                }
            } // end loop over j
            if (maxDistSqJ > maxDistSq) {
                maxI = i;
                maxJ = maxJ2;
                maxDistSq = maxDistSqJ;
            }
        } //end loop over i

        PairAndHull ph = new PairAndHull(maxI, maxJ, maxDistSq, ch);

        return ph;
    }

    public static class PairAndHull {
        public final int i0;
        public final int i1;
        public final double[] xHull;
        public final double[] yHull;
        public final double distSq;
        public PairAndHull(int i0, int i1, long distSq, CHL hull) {
            this.i0 = i0;
            this.i1 = i1;
            this.distSq = distSq;
            this.xHull = MatrixUtil.copyLongToDouble(hull.getXH());
            this.yHull = MatrixUtil.copyLongToDouble(hull.getYH());
        }
        public PairAndHull(int i0, int i1, double distSq, GrahamScan.CHD hull) {
            this.i0 = i0;
            this.i1 = i1;
            this.distSq = distSq;
            int n = hull.getXH().length;
            this.xHull = Arrays.copyOf(hull.getXH(), n);
            this.yHull = Arrays.copyOf(hull.getYH(), n);
        }
        public double[] getXY0() {
            return new double[]{xHull[i0], yHull[i0]};
        }
        public double[] getXY1() {
            return new double[]{xHull[i1], yHull[i1]};
        }

        /**
         *
         * @param coordStringFormat format for String.format to use in printing each coordinate of the hull.  e.g. "%.2f"
         * @return
         */
        public String toString(String coordStringFormat) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("furthest pair indexes = %d, %d\nhull:\n", i0, i1));
            String cFmt = String.format("(%s, %s)\n", coordStringFormat, coordStringFormat);
            for (int i = 0; i < xHull.length; ++i) {
                sb.append(String.format(cFmt, xHull[i], yHull[i]));
            }
            return sb.toString();
        }
    }

    /**
     * find the pair of points with the largest separation.  The runtime complexity is O(N*log_2(N)).
     * Note this method can be modified for a linear runtime option upon need.
     * @param points
     * @return
     * @throws GrahamScanTooFewPointsException
     */
    public static <T extends PairInt> T[] find(T[] points) throws GrahamScanTooFewPointsException {

        GrahamScanPairInt<T> scan = new GrahamScanPairInt<T>();
        scan.computeHull(points);

        List<T> hull = scan.getHull();

        int n = hull.size();

        long maxDistSq = Long.MIN_VALUE;
        int maxI = -1;
        int maxJ = -1;
        long maxDistSqJ;
        int maxJ2 = -1;

        long distSq;
        long xd;
        long yd;

        long xRef;
        long yRef;

        int nScan = n - 1;

        int j;
        int i;
        for (i = 0; i < nScan; ++i) {
            xRef = hull.get(i).getX();
            yRef = hull.get(i).getY();
            maxDistSqJ = Long.MIN_VALUE;

            for (j = i+1; j < nScan; ++j) {
                xd = hull.get(j).getX() - xRef;
                yd = hull.get(j).getY() - yRef;
                distSq = xd * xd + yd * yd;
                if (distSq >= maxDistSqJ) {
                    maxDistSqJ = distSq;
                    maxJ2 = j;
                } else {
                    break;
                }
            } // end loop over j
            if (maxDistSqJ > maxDistSq) {
                maxI = i;
                maxJ = maxJ2;
                maxDistSq = maxDistSqJ;
            }
        } //end loop over i

        // or specifiy type of class in method arguments
        T[] pairs = (T[]) Array.newInstance(PairInt.class, 2);
        pairs[0] = hull.get(maxI);
        pairs[1] = hull.get(maxJ);
        return pairs;
    }

}
