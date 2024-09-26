package algorithms.sort;

import algorithms.util.AngleUtil;

import java.util.Arrays;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/compGeometry/convexHull/
 * under MIT License (MIT), Nichole King 2013
 * 
 * The algorithm uses a merge sort which has a worse case runtime is O(N * log_2(N)),
 * but also includes an additional set of operations to remove inner points with the
 * same polar angle from P0.
 *
first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 */
public class PolarAngleMergeSort {

    public static int sort(float xP0, float yP0, float[] x, float[] y) {

        if (x == null) {
        	throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
        	throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
        	throw new IllegalArgumentException("number of items in x must be the same as in y");
        }
        
        if (x.length == 1) {
            return 1;
        }

        // for angles which are same, a delete operation is needed after all processing
        //    and ability to ignore the point to be deleted.
        double[] polarAngle = new double[x.length];
        
        for (int i = 1; i < x.length; i++) {
            polarAngle[i] = AngleUtil.polarAngleCCW((double)(x[i] - xP0), 
                (double)(y[i] - yP0));
        }
        
        sort(xP0, yP0, x, y, 1, x.length - 1, polarAngle);
        
        int nUsable = PolarAngleQuickSort.reduceToUniquePolarAngles(xP0, yP0, x, y, polarAngle);
        
        return nUsable;
    }
    
    protected static void sort(float xP0, float yP0, float[] x, float[] y, 
        int indexLo, int indexHi, double[] polarAngle) {

        if (indexLo < indexHi) {

            int indexMid = (indexLo + indexHi) >> 1;

            sort(xP0, yP0, x, y, indexLo, indexMid, polarAngle);
            sort(xP0, yP0, x, y, indexMid + 1, indexHi, polarAngle);
            merge(xP0, yP0, x, y, indexLo, indexMid, indexHi, polarAngle);
        }
    }
    
    private static void merge(float xP0, float yP0, float[] x, float[] y, 
        int indexLo, int indexMid, int indexHi, final double[] polarAngle) {

        int nLeft = indexMid - indexLo + 1;
        int nRight = indexHi - indexMid;

        float[] xLeft = Arrays.copyOfRange(x, indexLo, indexMid + 2);       // add 1 for sentinel
        float[] yLeft = Arrays.copyOfRange(y, indexLo, indexMid + 2);
        double[] angleLeft = Arrays.copyOfRange(polarAngle, indexLo, indexMid + 2);

        float[] xRight = Arrays.copyOfRange(x, indexMid + 1, indexHi + 2);  // add 1 for sentinel
        float[] yRight = Arrays.copyOfRange(y, indexMid + 1, indexHi + 2);
        double[] angleRight = Arrays.copyOfRange(polarAngle, indexMid + 1, indexHi + 2);

        xLeft[nLeft] = Float.MAX_VALUE;
        yLeft[nLeft] = Float.MAX_VALUE;
        angleLeft[nLeft] = Double.MAX_VALUE;
        xRight[nRight] = Float.MAX_VALUE;
        yRight[nRight] = Float.MAX_VALUE;
        angleRight[nRight] = Double.MAX_VALUE;

        int i = 0;
        int j = 0;

        for (int k = indexLo; k <= indexHi; k++) {

            if (angleLeft[i] <= angleRight[j]) {

                y[k] = yLeft[i];
                x[k] = xLeft[i];
                polarAngle[k] = angleLeft[i];
                i++;
            } else {

                y[k] = yRight[j];
                x[k] = xRight[j];
                polarAngle[k] = angleRight[j];
                j++;
            }
        }
    }

    static double relativeLengthOfLine(double x1, double y1, double x2, double y2) {
        double dx2 = (x2 - x1);
        dx2 *= dx2;
        double dy2 = (y2 - y1);
        dy2 *= dy2;
        //double d = Math.sqrt(dx2 + dy2);
        return dx2 + dy2;
    }


    /**
     * sort the points x, y by their polar angles with respect to the first point (x[0],y[0])
     * in counterclockwise order, ascending.
     * It uses a non-recursive merge sort for a runtime complexity of O(n*log_2(n)) and space complexity of O(n).
     * @param x
     * @param y
     * @param outPolarAngle the output array that should already be instantiated to same size
     *                      as x.  it will be filled with the polar angles of x,y from x[0],y[0]
     * @param reduceToUniqueAngles if true, for any points having the same polar angle, only the point furthest
     *                             from the first point in (x,y) will be kept. if false, all points are kept.
     * @param tolerance if reduceToUniqueAngles is true, tolerance is used to determing if two angles are equivalent.
     *                  The angles are in radians and so is the tolerance.
     * for reference, 0.1 degrees = 0.001745 radians.
     * @return the number of usable points in the sorted arrays. if reduceToUniqueAngles is true,
     * then for points having same angle,
     * only the furthest from x[0], y[0] is kept, and so there may be points at the ends of arrays
     * x,y, and outPolarAngle that are unusable.
     */
    public static int sortCCWBy1stPoint(double[] x, double[] y, double[] outPolarAngle, 
        boolean reduceToUniqueAngles, double tolerance) {

        if (x.length != outPolarAngle.length || x.length != y.length) {
            throw new IllegalArgumentException("x, y, and outPolarAngle must be same lengths");
        }

        double x0 = x[0];
        double y0 = y[0];

        double[] x1 = Arrays.copyOfRange(x, 1, x.length);
        double[] y1 = Arrays.copyOfRange(y, 1, y.length);
        double[] p1 = new double[x1.length];

        int n1 = x1.length;
        int i;
        for (i = 0; i < n1; i++) {
            // these are in radians
            p1[i] = Math.atan2((y1[i] - y0), (x1[i] - x0));
            if (p1[i] < 0) {
                p1[i] += 2. * Math.PI;
            }
            // to convert to degrees: multiply by (180/math.pi)
        }

        int[] sortedIndexes1 = MiscSorter.mergeSortIncreasing(p1);
        for (i = 0; i < n1; ++i) {
            x[i+1] = x1[sortedIndexes1[i]];
            y[i+1] = y1[sortedIndexes1[i]];
            outPolarAngle[i+1] = p1[i];
        }

        int nUsable;
        // need tolerance in radians.  for reference, 0.1 degrees = 0.001745 radians
        if (reduceToUniqueAngles) {
            nUsable = MiscSorter.reduceToUniquePolarAngles(x, y, outPolarAngle, tolerance);
        } else {
            nUsable = x.length;
        }

        return nUsable;
    }


    /**
     * sort the points x, y by their polar angles with respect to the first point (x[0],y[0])
     * in counterclockwise order, ascending.
     * It uses a non-recursive merge sort for a runtime complexity of O(n*log_2(n)) and space complexity of O(n).
     * @param x
     * @param y
     * @param outPolarAngle the output array that should already be instantiated to same size
     *                      as x.  it will be filled with the polar angles of x,y from x[0],y[0]
     * @param reduceToUniqueAngles if true, for any points having the same polar angle, only the point furthest
     *                             from the first point in (x,y) will be kept. if false, all points are kept.
     * @param tolerance if reduceToUniqueAngles is true, tolerance is used to determing if two angles are equivalent.
     *                  The angles are in radians and so is the tolerance.
     * for reference, 0.1 degrees = 0.001745 radians.
     * @return the number of usable points in the sorted arrays. if reduceToUniqueAngles is true,
     * then for points having same angle,
     * only the furthest from x[0], y[0] is kept, and so there may be points at the ends of arrays
     * x,y, and outPolarAngle that are unusable.
     */
    public static int sortCCWBy1stPoint(long[] x, long[] y, double[] outPolarAngle, boolean reduceToUniqueAngles,
                                        double tolerance) {

        if (x.length != outPolarAngle.length || x.length != y.length) {
            throw new IllegalArgumentException("x, y, and outPolarAngle must be same lengths");
        }

        long x0 = x[0];
        long y0 = y[0];

        long[] x1 = Arrays.copyOfRange(x, 1, x.length);
        long[] y1 = Arrays.copyOfRange(y, 1, y.length);
        double[] p1 = new double[x1.length];

        int n1 = x1.length;
        int i;
        double theta;
        for (i = 0; i < n1; i++) {
            // these are in radians
            theta = Math.atan2((y1[i] - y0), (x1[i] - x0));
            if (theta < 0) {
                theta += 2. * Math.PI;
            }
            p1[i] = theta;
            // to convert to degrees: multiply by (180/math.pi)
        }

        int[] sortedIndexes1 = MiscSorter.mergeSortIncreasing(p1);
        for (i = 0; i < n1; ++i) {
            x[i+1] = x1[sortedIndexes1[i]];
            y[i+1] = y1[sortedIndexes1[i]];
            outPolarAngle[i+1] = p1[i];
        }

        int nUsable;
        // need tolerance in radians.  for reference, 0.1 degrees = 0.001745 radians
        if (reduceToUniqueAngles) {
            nUsable = MiscSorter.reduceToUniquePolarAngles(x, y, outPolarAngle, tolerance);
        } else {
            nUsable = x.length;
        }

        return nUsable;
    }

}
