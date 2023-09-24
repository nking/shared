package algorithms.compGeometry.convexHull;

import algorithms.sort.MiscSorter;
import algorithms.sort.MultiArrayMergeSort;
import algorithms.sort.PolarAngleMergeSort;
import algorithms.util.Stack;

import java.util.Arrays;

/**
 * adapted from
 * https://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/compGeometry/convexHull/
 * under MIT License (MIT), Nichole King 2013
 *
  <pre>
  Solves the Convex Hull problem w/ a stack S of candidate points.

  Given a set of Q points returns the vertices of the ConvexHull(Q) in clockwise
  order.   a convex hull is the smallest convex polygon that will include all points in Q.

  Graham's Scan runs in O(n lg n).
    (in contrast to Jarvis's March which runs in O(nh) where h is the number of
    vertices in the convex hull.)
  Will adjust this after estimates...

  Both use a technique called 'rotational sweep' to process vertices in the order
  of the polar angles they form with a reference vertex.

  constructed from pseudo-code in Cormen et al. "Introduction to Algorithms
 </pre>

first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 * @author nichole
 */
public class GrahamScan {

    public static class CH {
        private float[] xH;
        private float[] yH;
        public CH(float[] x, float[] y) {
            this.xH = x;
            this.yH = y;
        }

        /**
         * @return the xH
         */
        public float[] getXH() {
            return xH;
        }

        /**
         * @return the yH
         */
        public float[] getYH() {
            return yH;
        }

        @Override
        public String toString() {
            if (xH == null) {
                return "[]";
            }
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < xH.length; ++i) {
                if (i > 0) {
                    sb.append(", ");
                }
                sb.append("(").append(xH[i]).append(",").append(yH[i]).append(")");
            }
            return sb.toString();
        }
    }

    public static class CHD {
        private double[] xH;
        private double[] yH;
        public CHD(double[] x, double[] y) {
            this.xH = x;
            this.yH = y;
        }

        /**
         * @return the xH
         */
        public double[] getXH() {
            return xH;
        }

        /**
         * @return the yH
         */
        public double[] getYH() {
            return yH;
        }

        @Override
        public String toString() {
            if (xH == null) {
                return "[]";
            }
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < xH.length; ++i) {
                if (i > 0) {
                    sb.append(", ");
                }
                sb.append("(").append(xH[i]).append(",").append(yH[i]).append(")");
            }
            return sb.toString();
        }
    }

    public static class CHL {
        private long[] xH;
        private long[] yH;
        public CHL(long[] x, long[] y) {
            this.xH = x;
            this.yH = y;
        }

        /**
         * @return the xH
         */
        public long[] getXH() {
            return xH;
        }

        /**
         * @return the yH
         */
        public long[] getYH() {
            return yH;
        }

        @Override
        public String toString() {
            if (xH == null) {
                return "[]";
            }
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < xH.length; ++i) {
                if (i > 0) {
                    sb.append(", ");
                }
                sb.append("(").append(xH[i]).append(",").append(yH[i]).append(")");
            }
            return sb.toString();
        }
    }

	public GrahamScan() {
	}

    /**
     * find the convex hull of the given (x, y) points.  Note that the resulting
     * hull points have the same first point as last point.
     * The runtime complexity is O(N*log_2(N)) if useLinearSort is false, wlse is O(N) if useLinearSort is true.
     * @param x
     * @param y
     * @param useLinearSort if true, a linear runtime complexity algorithm is used (CountingSort or BucketSort
     *                      depending upon the range of data and number of points.  If the range of data
     *                      is 0 to 359 and x.length >= 360, then counting sort is used with a resolution of 1 degree.
     * @throws GrahamScanTooFewPointsException
     * @return the convex hull of (x, y)
     */
    public CHD computeHull(float[] x, float[] y, boolean useLinearSort) throws GrahamScanTooFewPointsException {

        if (x == null) {
            throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
            throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
            throw new IllegalArgumentException("x must have the same number of items as y");
        }
        if (x.length < 3) {
            throw new IllegalArgumentException("x must have at least 3 items");
        }

        double[] x2 = new double[x.length];
        double[] y2 = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            x2[i] = x[i];
            y2[i] = y[i];
        }

        CHD chd = computeHull(x2, y2, useLinearSort);

        return chd;
    }

    /**
     * find the convex hull of the given (x, y) points.  Note that the resulting
     * hull points have the same first point as last point.
     * The runtime complexity is O(N*log_2(N)) if useLinearSort is false, wlse is O(N) if useLinearSort is true.
     * @param x
     * @param y
     * @param useLinearSort if true, a linear runtime complexity algorithm is used (CountingSort or BucketSort
     *                      depending upon the range of data and number of points.  If the range of data
     *                      is 0 to 359 and x.length >= 360, then counting sort is used with a resolution of 1 degree.
     * @throws GrahamScanTooFewPointsException
     * @return the convex hull of (x, y)
     */
    public CHD computeHull(double[] x, double[] y, boolean useLinearSort) throws GrahamScanTooFewPointsException {

        if (x == null) {
            throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
            throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
            throw new IllegalArgumentException("x must have the same number of items as y");
        }
        if (x.length < 3) {
            throw new IllegalArgumentException("x must have at least 3 items");
        }

        // make copy of points to avoid modifying input
        x = Arrays.copyOf(x, x.length);
        y = Arrays.copyOf(y, y.length);

        /*
         * Q is a stack of candidate points which have been pushed once onto the stack
         * and removed if they are not vertices of the stack.
         *
         * when complete, the stack S contains the vertices of the hull in counterclockwise order.
         *
         * Q > 3
         *
         * 1 -- let p0 be the point in Q w/ min y-coordinate, or leftmost point of a tie
         * 2 -- let <p1, p2, ... pm> be the remaining points in Q.
         *      sorted by polar angle in counter clockwise order around p0.
         *      ** if more than one point has the same angle, remove all but the one that is furthest from p0. **
         * 3 -- push p0 onto S
         * 4 -- push p1 onto S
         * 5 -- push p2 onto S
         * 6 -- for i=3 to m
         * 7 --     do while the angle formed by points NEXT-TO-TOP(S), TOP(S), and p_i makes a nonleft turn
         * 8 --         pop(S)
         * 9 --     push(S)
         * 10 -return S
         */

        // (1) let p0 be the point in Q w/ minimum yCoord,
        //     or the leftmost point if more than one w/ same minimum yCoord.
        int p0Index = findIndexOfMinY(x, y);

        if (p0Index != 0) {
            // move the point at index [iP0] to index [0]
            double swapX = x[p0Index];
            double swapY = y[p0Index];
            x[p0Index] = x[0];
            y[p0Index] = y[0];
            x[0] = swapX;
            y[0] = swapY;
            p0Index = 0;
        }

        // (2) let <p1, p2, ..., pm> be the remaining points in Q, sorted
        //     by polar angle in counterclockwise order around p0
        //     (if more than one pt has same angle, keep only the furthest from p0)

        // this step uses angles rounded to degrees between 0 and 360.
        // the runtime complexity is O( max(x.length, 360) )
        int nPointsUsable;
        boolean reduceToUniqueAngles = true;
        if (useLinearSort) {
            nPointsUsable = MiscSorter.sortCCWBy1stPointLinear(x, y, reduceToUniqueAngles);
        } else {
            double[] outputPolarAngle = new double[x.length];
            double tolerance = 1E-4;
            nPointsUsable = PolarAngleMergeSort.sortCCWBy1stPoint(x, y, outputPolarAngle, reduceToUniqueAngles, tolerance);
        }

        if (nPointsUsable < 3) {
            throw new GrahamScanTooFewPointsException("polar angle sorting has reduced the number of points to less than 3");
        }

        Stack<Integer> s = new Stack<Integer>();
        s.push(p0Index);
        s.push(1);
        s.push(2);

        int topS;
        int nextToTopS;
        int si;

        // for i = 3 to m
        //    while angle between next-to-top(S), top(S) and p_i makes a nonleft turn
        //        do pop(S)
        //    push(pi, S)
        for (int i = 3; i < nPointsUsable; i++) {

            topS = s.peek();
            nextToTopS = s.peekPopNext();
            si = i;

            double direction = ((x[topS] - x[nextToTopS])*(y[si] - y[nextToTopS])) - ((y[topS]
                    - y[nextToTopS])*(x[si] - x[nextToTopS]));

            while (direction <= 0) {

                s.pop();

                if (s.size() < 2) {
                    break;
                }

                topS = s.peek();
                nextToTopS = s.peekPopNext();
                si = i;
                direction = ((x[topS] - x[nextToTopS])*(y[si] - y[nextToTopS])) - ((y[topS]
                        - y[nextToTopS])*(x[si] - x[nextToTopS]));
            }

            s.push(i);
        }

        int n = s.size() + 1;

        double[] xHull = new double[n];
        double[] yHull = new double[n];
        for (int i = 0; i < (n - 1); ++i) {
            si = s.pop();
            xHull[i] = x[si];
            yHull[i] = y[si];
        }
        xHull[n-1] = xHull[0];
        yHull[n-1] = yHull[0];

        return new CHD(xHull, yHull);
    }

    /**
     * find the convex hull of the given (x, y) points.  Note that the resulting
     * hull points have the same first point as last point.
     * The runtime complexity is O(N*log_2(N)) if useLinearSort is false, wlse is O(N) if useLinearSort is true.
     * @param x
     * @param y
     * @param useLinearSort if true, a linear runtime complexity algorithm is used (CountingSort or BucketSort
     *                      depending upon the range of data and number of points.  If the range of data
     *                      is 0 to 359 and x.length >= 360, then counting sort is used with a resolution of 1 degree.
     * @throws GrahamScanTooFewPointsException
     * @return the convex hull of (x, y)
     */
    public static CHL computeHull(long[] x, long[] y, boolean useLinearSort) throws GrahamScanTooFewPointsException {

        if (x == null) {
            throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
            throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
            throw new IllegalArgumentException("x must have the same number of items as y");
        }
        if (x.length < 3) {
            throw new IllegalArgumentException("x must have at least 3 items");
        }

        // make copy of points to avoid modifying input
        x = Arrays.copyOf(x, x.length);
        y = Arrays.copyOf(y, y.length);

        /*
         * Q is a stack of candidate points which have been pushed once onto the stack
         * and removed if they are not vertices of the stack.
         *
         * when complete, the stack S contains the vertices of the hull in counterclockwise order.
         *
         * Q > 3
         *
         * 1 -- let p0 be the point in Q w/ min y-coordinate, or leftmost point of a tie
         * 2 -- let <p1, p2, ... pm> be the remaining points in Q.
         *      sorted by polar angle in counter clockwise order around p0.
         *      ** if more than one point has the same angle, remove all but the one that is furthest from p0. **
         * 3 -- push p0 onto S
         * 4 -- push p1 onto S
         * 5 -- push p2 onto S
         * 6 -- for i=3 to m
         * 7 --     do while the angle formed by points NEXT-TO-TOP(S), TOP(S), and p_i makes a nonleft turn
         * 8 --         pop(S)
         * 9 --     push(S)
         * 10 -return S
         */

        // (1) let p0 be the point in Q w/ minimum yCoord,
        //     or the leftmost point if more than one w/ same minimum yCoord.
        int p0Index = findIndexOfMinY(x, y);

        if (p0Index != 0) {
            // move the point at index [iP0] to index [0]
            long swapX = x[p0Index];
            long swapY = y[p0Index];
            x[p0Index] = x[0];
            y[p0Index] = y[0];
            x[0] = swapX;
            y[0] = swapY;
            p0Index = 0;
        }

        // (2) let <p1, p2, ..., pm> be the remaining points in Q, sorted
        //     by polar angle in counterclockwise order around p0
        //     (if more than one pt has same angle, keep only the furthest from p0)

        // this step uses angles rounded to degrees between 0 and 360.
        // the runtime complexity is O( max(x.length, 360) )
        int nPointsUsable;
        boolean reduceToUniqueAngles = true;
        if (useLinearSort) {
            nPointsUsable = MiscSorter.sortCCWBy1stPointLinear(x, y, reduceToUniqueAngles);
        } else {
            double[] outputPolarAngle = new double[x.length];
            double tolerance = 1E-4;
            nPointsUsable = PolarAngleMergeSort.sortCCWBy1stPoint(x, y, outputPolarAngle, reduceToUniqueAngles, tolerance);
        }

        if (nPointsUsable < 3) {
            throw new GrahamScanTooFewPointsException("polar angle sorting has reduced the number of points to less than 3");
        }

        Stack<Integer> s = new Stack<Integer>();
        s.push(p0Index);
        s.push(1);
        s.push(2);

        int topS;
        int nextToTopS;
        int si;

        // for i = 3 to m
        //    while angle between next-to-top(S), top(S) and p_i makes a nonleft turn
        //        do pop(S)
        //    push(pi, S)
        for (int i = 3; i < nPointsUsable; i++) {

            topS = s.peek();
            nextToTopS = s.peekPopNext();
            si = i;

            double direction = ((x[topS] - x[nextToTopS])*(y[si] - y[nextToTopS])) - ((y[topS]
                    - y[nextToTopS])*(x[si] - x[nextToTopS]));

            while (direction <= 0) {

                s.pop();

                if (s.size() < 2) {
                    break;
                }

                topS = s.peek();
                nextToTopS = s.peekPopNext();
                si = i;
                direction = ((x[topS] - x[nextToTopS])*(y[si] - y[nextToTopS])) - ((y[topS]
                        - y[nextToTopS])*(x[si] - x[nextToTopS]));
            }

            s.push(i);
        }

        int n = s.size() + 1;

        long[] xHull = new long[n];
        long[] yHull = new long[n];
        for (int i = 0; i < (n - 1); ++i) {
            si = s.pop();
            xHull[i] = x[si];
            yHull[i] = y[si];
        }
        xHull[n-1] = xHull[0];
        yHull[n-1] = yHull[0];

        return new CHL(xHull, yHull);
    }

    protected static int findIndexOfMinY(double[] x, double[] y) {
        if (x.length < 1) {
            throw new IllegalArgumentException("x and y must be longer than 1");
        }
        int iMin = 0;
        double minY = y[0];
        double minX = x[0];
        int i;
        for (i = 1; i < x.length; ++i) {
            if ((y[i] < minY) || ( (y[i] == minY) && (x[i] < minX))) {
                minY = y[i];
                minX = x[i];
                iMin = i;
            }
        }
        return iMin;
    }

    protected static int findIndexOfMinY(long[] x, long[] y) {
        if (x.length < 1) {
            throw new IllegalArgumentException("x and y must be longer than 1");
        }
        int iMin = -1;
        long minY = Integer.MAX_VALUE;
        long minX = -1;
        int i;
        for (i = 0; i < x.length; ++i) {
            if ((y[i] < minY) || (y[i] == minY && x[i] < minX)){
                minY = y[i];
                minX = x[i];
                iMin = i;
            }
        }
        return iMin;
    }

}
