package algorithms.misc;

import algorithms.sort.CountingSort;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.List;
import java.util.ArrayList;
import java.util.Collection;
import gnu.trove.set.TIntSet;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.Arrays;
import algorithms.sort.MiscSorter;
import gnu.trove.list.array.TIntArrayList;
import java.math.RoundingMode;

/**
    miscellaneous math methods. some could probably be improved.

    first implemented in projects
     https://github.com/nking/two-point-correlation
     w/ Copyright (c) 2013-2015 Nichole King
     http://nking.github.io/two-point-correlation/
     using The MIT License (MIT)
     and
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 * @author nichole
 */
public class MiscMath0 {

    public static float[] calculateOuterRoundedMinAndMax(float[] a) {

        // find the powers of 10 for the data min and max
        float xmin = Float.MAX_VALUE;
        float xmax = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > xmax) {
                xmax = a[i];
            }
            if (a[i] < xmin) {
                xmin = a[i];
            }
        }

        xmax = MiscMath0.roundUpByLargestPower(xmax);

        xmin = MiscMath0.roundDownByLargestPower(xmin);

        // xmax > 1 and xmin is between 0 and 1, round xmin down
        if ((xmax > 1) && (xmin > 0) && (xmin < 1.0)) {
            xmin = 0;
        }

        return new float[]{xmin, xmax};
    }

    /**
     Round number up to the next digit in the largest power.

        MiscMath0.roundUpByLargestPower(31.1f) == 40.0f;

        MiscMath0.roundUpByLargestPower(0.11f) == 1.0f;

        MiscMath0.roundUpByLargestPower(-0.011f) == -0.02f;

        MiscMath0.roundUpByLargestPower(310.1f) == 400.0f;

        MiscMath0.roundUpByLargestPower(-3.1) == -4.0f;

        MiscMath0.roundUpByLargestPower(3.1) == 4.0f;

    @param f
    @return
    */
    public static float roundUpByLargestPower(float f) {

        int pow = findPowerOf10(f);

        double pow10 = Math.pow(10, pow);

        int d = (int)(f/pow10);
        int m;
        if (f == pow10) {
            m = d;
        } else if (f > 0) {
            // residual ?
            if (f > (d*pow10)) {
                m = d + 1;
            } else {
                m = d;
            }
        } else {
            // decimals
            m = d - 1;
        }

        float r = (float) (m * pow10);

        return r;
    }

    /**
     * Round number down to in the largest power.
     * For example,
     *     roundDownByLargestPower(3.1) returns 3.0
     *     roundDownByLargestPower(-3.1) returns -4.0
     *     roundDownByLargestPower(31.1) returns 30.0
     *
     * @param f
     * @return
     */
    public static float roundDownByLargestPower(float f) {

        if (f == 0) {
            return 0;
        }

        int pow = findPowerOf10(f);
        double pow10 = Math.pow(10, pow);

        double r;
        if (f > 0) {
            double m = f % pow10;
            r = f - m;
        } else {
            int m = (int)(f/pow10);
            r = pow10 * (m - 1);
        }

        return (float)r;
    }
    
    public static int findPowerOf10(float a) {

        if (a == 0) {
            return 0;
        }
        if (a < 0.f) {
            a *= -1.0f;
        }
        double b = Math.log10(a);
        if (b > 0) {
            return (int)b;
        } else {
            if (b >= 1) {
                return (int)Math.round(b);
            } else if (b > -1) {
                // fractions between -1 and +1
                 return findPowerOf10_2(a);
            } else {
                return (int)Math.round(b);
            }
        }
    }

    public static int findPowerOf10_2(float a) {

        if (a == 0) {
            return 0;
        }

        int power = 0;

        if (a <= 1.0f) {
            while (a < 1.0) {
                a *=  10.0;
                power--;
            }
        } else {
            // precision errors in multiplication here are trouble for non base2 numbers such as powers of 10
            while (a >= 1.0) {
                a /= 10.0;
                power++;
            }
            power--;
        }

        return power;
    }

    /**
     * find max but ignore values such as FLOAT.MAX_VALUE, infinity, and NAN
     * @param a
     * @return
     */
    public static int findYMaxIndex(float[] a) {
        if (a == null || a.length == 0) {
            return -1;
        }
        float max = Float.NEGATIVE_INFINITY;
        int index = 0;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) 
                && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }
    
    /**
     * find max but ignore values such as infinity, and NAN
     * @param a
     * @return
     */
    public static int findYMaxIndex(double[] a) {
        if (a == null || a.length == 0) {
            return -1;
        }
        double max = Double.NEGATIVE_INFINITY;
        int index = 0;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && !Double.isInfinite(a[i]) && !Double.isNaN(a[i]) 
                ) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }

    public static int findYMaxIndex(int[] a) {
        if (a == null || a.length == 0) {
            return -1;
        }
        int max = Integer.MIN_VALUE;
        int index = 0;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }

    public static int findYMaxIndex(TIntList a) {
        if (a == null || a.size() == 0) {
            return -1;
        }
        int max = Integer.MIN_VALUE;
        int index = 0;
        for (int i = 0; i < a.size(); i++) {
            if (a.get(i) > max) {
                max = a.get(i);
                index = i;
            }
        }
        return index;
    }

    public static int findYMaxIndex(long[] a) {
    
        if (a == null || a.length == 0) {
            return -1;
        }
        long max = Long.MIN_VALUE;
        int index = 0;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }
 
    public static float findMin(float[] a) {
        float min = Float.POSITIVE_INFINITY;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] < min) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i])) {
                min = a[i];
            }
        }
        return min;
    }
    
    public static double findMin(double[] a) {
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] < min) && !Double.isInfinite(a[i]) && !Double.isNaN(a[i])) {
                min = a[i];
            }
        }
        return min;
    }

    public static int findMin(int[] a) {
        int min = Integer.MAX_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] < min) {
                min = a[i];
            }
        }
        return min;
    }

    /**
     * find max 
     * @param a
     * @return
     */
    public static int findMax(int[] a) {
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
            }
        }
        return max;
    }
    
    public static int findMax(int[][] a) {
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                if (a[i][j] > max) {
                    max = a[i][j];
                }
            }
        }
        return max;
    }
    
    public static int findMin(int[][] a) {
        int min = Integer.MAX_VALUE;
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                if (a[i][j] < min) {
                    min = a[i][j];
                }
            }
        }
        return min;
    }
    
    public static float findMax(float[][] a) {
        float max = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                if (a[i][j] > max) {
                    max = a[i][j];
                }
            }
        }
        return max;
    }
    
    public static float findMin(float[][] a) {
        float min = Float.POSITIVE_INFINITY;
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                if (a[i][j] < min) {
                    min = a[i][j];
                }
            }
        }
        return min;
    }

    /**
     * find max 
     * @param a
     * @return
     */
    public static float findMax(float[] a) {
        float max = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
            }
        }
        return max;
    }
    
    public static double findMax(double[] a) {
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
            }
        }
        return max;
    }
    
    public static double findMin(double[][] img) {
        
        double min = Double.POSITIVE_INFINITY;
        
        for (int i = 0; i < img.length; ++i) {
            for (int j = 0; j < img[i].length; ++j) {
                double v = img[i][j];
                if (v < min) {
                    min = v;
                }
            }
        }
        
        return min;
    }
    
    public static double findMax(double[][] img) {
        
        double max = Double.NEGATIVE_INFINITY;
        
        for (int i = 0; i < img.length; ++i) {
            for (int j = 0; j < img[i].length; ++j) {
                double v = img[i][j];
                if (v > max) {
                    max = v;
                }
            }
        }
        
        return max;
    }

    /**
     * for the given histogram, returns the indexes of the primary peak and
     * any peaks which are larger than frac*maxPeak above their neighboring
     * values.
     * @param h
     * @param frac the fraction of the y value of the maximum peak that is used
     * as a critical limit that any other peaks must have in excess of their
     * neighboring points.  For example, 0.1.
     * @return 
     */
    public static List<Integer> findStrongPeakIndexes(HistogramHolder h, float frac) {
        
        float[] x = h.getXHist();
        int[] y = h.getYHist();
        
        if (x == null || y == null || x.length == 0 || y.length == 0) {
            return null;
        }
        
        if (x.length == 1) {
            List<Integer> minMaxIndexes = new ArrayList<Integer>();
            if (y[0] > 0) {
                minMaxIndexes.add(Integer.valueOf(0));
            }
            return minMaxIndexes;
        }
        
        int yPeakIdx = MiscMath0.findYMaxIndex(y);
        
        if (yPeakIdx == -1) {
            return null;
        }
        
        int yMaxPeak = y[yPeakIdx];
        
        /*
        storing the minima and maxima in the same array list.
        the minima have -1*index within k
        and the maxima keep their positive values of the index within k.
        */
        List<Integer> minMaxIndexes = new ArrayList<Integer>();
        
        float lastY = y[0];
        boolean incr = true;
        for (int ii = 1; ii < y.length; ii++) {
            if ((y[ii] < lastY) && incr) {
                minMaxIndexes.add(Integer.valueOf(ii - 1));
                incr = false;
            } else if ((y[ii] > lastY) && !incr) {
                minMaxIndexes.add(Integer.valueOf(-1*(ii - 1)));
                incr = true;
            }
            lastY = y[ii];
        }
        
        if (incr) {
            // add the last point
             minMaxIndexes.add(Integer.valueOf(y.length - 1));
        }
        
        // for the histograms of euclidean point combination differences,
        // this should usually be singly peaked
        if (minMaxIndexes.size() == 1) {
            return minMaxIndexes;
        }
        
        float limit = frac * yMaxPeak;
        
        // find peaks where y[ii] is > limit above adjacent local minima
        
        List<Integer> peaks = new ArrayList<Integer>();

        for (int ii = 0; ii < minMaxIndexes.size(); ii++) {

            int idx = minMaxIndexes.get(ii).intValue();

            if (idx > -1) {
                // this is a maxima
                
                boolean found = false;
                
                // compare to preceding minimum
                for (int iii = (ii - 1); iii > -1; iii--) {
                    int idx2 = minMaxIndexes.get(iii).intValue();
                    if (idx2 < 0) {
                        float compare = y[-1*idx2];
                        float diff = y[idx] - compare;
                        if (diff >= limit) {
                            peaks.add(Integer.valueOf(idx));
                            found = true;
                        }
                        break;
                    }
                }
                if (found) {
                    continue;
                }

                //compare to proceeding minimum
                for (int iii = (ii + 1); iii < minMaxIndexes.size(); iii++) {
                    int idx2 = minMaxIndexes.get(iii).intValue();
                    if (idx2 < 0) {
                        float compare = y[-1*idx2];
                        float diff = y[idx] - compare;
                        if (diff >= limit) {
                            peaks.add(Integer.valueOf(idx));
                        }
                        break;
                    }
                }
            }
        }
       
        return peaks;
    }

    public static List<Integer> findStrongPeakIndexesDescSort(
        HistogramHolder hist, float fracMax) {

        List<Integer> indexes = findStrongPeakIndexes(hist, fracMax);
        
        if (indexes == null || indexes.size() < 2) {
            return indexes;
        }
        
        int[] idxs = new int[indexes.size()];
        int[] c = new int[idxs.length];
        
        int maxC = Integer.MIN_VALUE;
        for (int i = 0; i < indexes.size(); ++i) {
            idxs[i] = indexes.get(i).intValue();
            c[i] = hist.getYHist()[idxs[i]];
            if (c[i] > maxC) {
                maxC = c[i];
            }
        }
        
        CountingSort.sortByDecr(c, idxs);
        
        indexes.clear();
        
        for (int i = 0; i < idxs.length; ++i) {
            indexes.add(Integer.valueOf(idxs[i]));
        }
        
        return indexes;
    }

    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return 
     */
    public static double[] getAvgAndStDev(int[] x) {
        
        return getAvgAndStDev(x, x.length);
    }

    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return 
     */
    public static double[] getAvgAndStDev(int[] x, int length) {
        
        long sumX = 0;
        for (int i = 0; i < length; i++) {
            sumX += x[i];
        }
        
        double avgX = (double)sumX/(double)length;
        
        sumX = 0;
        for (int i = 0; i < length; i++) {
            double diffX = x[i] - avgX;
            sumX += (diffX * diffX);
        }
        double stdDevX = Math.sqrt(sumX/(length - 1.0f));
        
        return new double[]{avgX, stdDevX};
    }
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return 
     */
    public static double[] getAvgAndStDev(long[] x, int length) {
        
        long sumX = 0;
        for (int i = 0; i < length; i++) {
            sumX += x[i];
        }
        
        double avgX = (double)sumX/(double)length;
        
        double diffX;
        sumX = 0;
        for (int i = 0; i < length; i++) {
            diffX = x[i] - avgX;
            sumX += (diffX * diffX);
        }
        double stdDevX = Math.sqrt(sumX/(length - 1.0f));
        
        return new double[]{avgX, stdDevX};
    }
    
    /**
     * calculate the quartiles Q1, Q2, and Q3 of x as
     * <pre>
     * x = ascending sort(x)
     * n = x.length;
     * Q1 = median of x[0] through x[n/2] which is x[n/4];
     * Q2 = median of x = x[n/2];
     * Q3 = median of x[n/2] through x[n-1] = x[(int)0.75*n];
     * </pre>
     * The IQR is Q3 - Q1;
     * @param x
     * @return an array holding Q1, Q2, Q3, min and max
     */
    public static double[] calculateQuartiles(double[] x) {
        
        x = Arrays.copyOf(x, x.length);
        Arrays.sort(x);
        
        int n = x.length;
        
        double[] q = new double[]{
            x[n/4], x[n/2], x[(int)(0.75*n)], x[0], x[n-1]
        };
        
        return q;
    }
    
    /**
    uses Tukey fences to find inlier indexes in x with factor k=1.5
    <pre>
    calculated as
       q = calculateQuartiles(x);
       iqr = q[2]-q[0];
       inliers are in the range [q[0] - k*iqr, q[2] + k*iqr], inclusive.
     * @param x
     * @return 
    */
    public static int[] findInliersUsingTukeyFences(double[] x) {
        return findInliersUsingTukeyFences(x, 1.5);
    }
    
    /**
    uses Tukey fences to find inlier indexes in x with factor k=3
    <pre>
    calculated as
       q = calculateQuartiles(x);
       iqr = q[2]-q[0];
       inliers are in the range [q[0] - k*iqr, q[2] + k*iqr], inclusive.
     * @param x
     * @return 
    */
    public static int[] findFarInliersUsingTukeyFences(double[] x) {
        return findInliersUsingTukeyFences(x, 3);
    }
    
    /**
    uses Tukey fences to find inlier indexes in x.
    <pre>
    calculated as
       q = calculateQuartiles(x);
       iqr = q[2]-q[0];
       inliers are in the range [q[0] - k*iqr, q[2] + k*iqr], inclusive.
     * @param x
     * @param k
     * @return 
    */
    public static int[] findInliersUsingTukeyFences(double[] x, double k) {
        // from https://en.m.wikipedia.org/wiki/Outlier
        double[] q = calculateQuartiles(x);
        int n = x.length;
        int[] inlierIndexes = new int[n];
        int i;
        int count = 0;
        double iqr = q[2] - q[0];
        double r0 = q[0] - k*iqr;
        double r1 = q[2] + k*iqr;
        for (i = 0; i < n; ++i) {
            if (x[i] >= r0 && x[i] <= r1) {
                inlierIndexes[count] = i;
                count++;
            }
        }
        return Arrays.copyOf(inlierIndexes, count);
    }
    
    /**
     * calculate the median of absolute deviation (MAD) as
     * <pre>
     *  absDev_i = abs(x_i - median(x))
        mad = median( absDev)
     * </pre>
     * and returns the median, followed by the median of absolute deviation.
     * to use the MAD as one uses standard deviation in determining outliers,
     * use stDev = k*MAD where k is a constant scale factor.  
     * For gaussian k~1.4826.
       Then and outlier is outside the range [median - 3*MAD, median + 3*MAD].
     * @param x
     * @return an array holding the median of absolute deviation of x, 
     * the median, the min, and the max.
     */
    public static double[] calculateMedianOfAbsoluteDeviation(double[] x) {
        
        x = Arrays.copyOf(x, x.length);
        Arrays.sort(x);
        
        int n = x.length;
        
        double median = x[n/2];
        
        double[] d = new double[n];
        
        int i;
        for (i = 0; i < n; ++i) {
            d[i] = Math.abs(x[i] - median);
        }
        
        Arrays.sort(d);
        
        return new double[]{d[n/2], median, x[0], x[n-1]};
    }
    
    /**
     * calculate the median and the interquartile range
     * @param x
     * @return 
     */
    public static double[] calcMedianAndIQR(double[] x) {
        x = Arrays.copyOf(x, x.length);
        MiscSorter.mergeSortIncreasing(x);
        double[] r = new double[2];
        int n = x.length;
        r[0] = x[n/2];
        r[1] = (x[3*n/4] - x[n/2]);
        return r;
    }
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return 
     */
    public static double[] getAvgAndStDev(double[] x) {
        
        double sumX = 0;
        for (int i = 0; i < x.length; i++) {
            sumX += x[i];
        }
        
        double avgX = (double)sumX/(double)x.length;
        
        sumX = 0;
        for (int i = 0; i < x.length; i++) {
            double diffX = x[i] - avgX;
            sumX += (diffX * diffX);
        }
        double stdDevX = Math.sqrt(sumX/(x.length - 1.0));
        
        return new double[]{avgX, stdDevX};
    }    

    public static PairIntArray get20NeighborOffsets() {
        
        PairIntArray r2Offsets = new PairIntArray();
        
                              r2Offsets.add(-1, 2); r2Offsets.add(0, 2); r2Offsets.add(1, 2);
        r2Offsets.add(-2, 1); r2Offsets.add(-1, 1); r2Offsets.add(0, 1); r2Offsets.add(1, 1); r2Offsets.add(2, 1);
        r2Offsets.add(-2, 0); r2Offsets.add(-1, 0);                       r2Offsets.add(1, 0); r2Offsets.add(2, 0);
        r2Offsets.add(-2, -1); r2Offsets.add(-1, -1); r2Offsets.add(0, -1); r2Offsets.add(1, -1); r2Offsets.add(2, -1);
                              r2Offsets.add(-1, -2); r2Offsets.add(0, -2); r2Offsets.add(1, -2);
        
        return r2Offsets;
    }

     /**
     *
     * @param a
     * @return
     */
    public static int[] findMinMaxValues(int[][] a) {
        
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                int v = a[i][j];
                if (v < min) {
                    min = v;
                }
                if (v > max) {
                    max = v;
                }
            }
        }
        
        return new int[]{min, max};
    }

    /**
     * find the minima and maxima of x and y and return them as
     * int[]{xMin, xMax, yMin, yMax}
     * @param points
     * @return minMaxXY int[]{xMin, xMax, yMin, yMax}
     */
    public static int[] findMinMaxXY(Collection<PairInt> points) {
        
        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;
        int yMin = Integer.MAX_VALUE;
        int yMax = Integer.MIN_VALUE;
        
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            if (x < xMin) {
                xMin = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y > yMax) {
                yMax = y;
            }
        }
        return new int[]{xMin, xMax, yMin, yMax};
    }
    
    /**
     * find the minima and maxima of x and y and return them as
     * int[]{xMin, xMax, yMin, yMax}
     * @return minMaxXY int[]{xMin, xMax, yMin, yMax}
     */
    public static int[] findMinMaxXY(TIntSet pixelIdxs, int imgWidth) {
        
        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;
        int yMin = Integer.MAX_VALUE;
        int yMax = Integer.MIN_VALUE;
        
        TIntIterator iter = pixelIdxs.iterator();
        
        while (iter.hasNext()) {
            
            int pixIdx = iter.next();
            int y = pixIdx/imgWidth;
            int x = pixIdx - (y * imgWidth);
            if (x < xMin) {
                xMin = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y > yMax) {
                yMax = y;
            }
        }
        return new int[]{xMin, xMax, yMin, yMax};
    }

    /**
     * find the minima and maxima of x and y and return them as
     * int[]{xMin, xMax, yMin, yMax}
     * @param points
     * @return minMaxXY int[]{xMin, xMax, yMin, yMax}
     */
    public static int[] findMinMaxXY(PairIntArray points) {
        
        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;
        int yMin = Integer.MAX_VALUE;
        int yMax = Integer.MIN_VALUE;
        
        for (int i = 0; i < points.getN(); ++i) {
            int x = points.getX(i);
            int y = points.getY(i);
            if (x < xMin) {
                xMin = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y > yMax) {
                yMax = y;
            }
        }
        return new int[]{xMin, xMax, yMin, yMax};
    }
    
    /**
     * test for whether n is a power of 2.
     * @param n a non-negative number.
     * @return true if is a power of , else false.
     */
    public static boolean isAPowerOf2(int n) {
        // bitmask test e.g. 128 & 0x7f = 0
        return ((n == 0) || ((n & (n - 1)) == 0));
    }
    
    /**
     * convert integer a to the given base. 
     * java.lang.Integer already has radix toString operations,
     * so this is just here for convenience.
     * @param a
     * @param base
     * @return 
     */
    public static String convertToBase(int a, int base) {
        StringBuilder sb = new StringBuilder();
        while (a > 0) {
            sb.append(Integer.toString(a % base));
            a /= base;
        }
        sb.reverse();
        return sb.toString();
    }
    
    /**
     * get fractional part of a.  
     * e.g. for a = 385.55, fractional part = 0.55.
     * for a = -385.55, fractional part = -0.55.
     * @param a 
     * @return 
     */
    public static double getFractionalPart(double a) {
        if (a < 0) {
            return -1*(-a % 1);
        }
        return a % 1;
    }
    
    /**
     * determine the number of bits, that is, the msb position + 1.
     * Note that a value of 0 returns a bit length of 1.
     * @param v
     * @return 
     */
    public static int numberOfBits(int v) {
        
        if (v < 0) {
            v *= -1;
        } else if (v == 0) {
            return 1;
        }
        return 32 - Integer.numberOfLeadingZeros(v);
    }
    
    /**
     * determine the number of bits, that is the msb position + 1.
     * Note that a value of 0 returns a bit length of 1.
     * @param v
     * @return 
     */
    public static int numberOfBits(long v) {
        if (v < 0) {
            v *= -1;
        } else if (v == 0) {
            return 1;
        }
        return 64 -  Long.numberOfLeadingZeros(v);
    }

    /**
     * determine the number of set bits.  this method uses the
     * hamming weight which uses binary magic numbers.
     <pre>
     reference:
     https://en.wikipedia.org/wiki/Hamming_weight
     method popcount64c
     </pre>
     * @param x the bitstring with the set bits to count
     * @return the number of bits set in x
     */
    public static int numberOfSetBits(long x) {
        //put count of each 2 bits into those 2 bits, where the mask is 62 bits of repeated '10's
        x -= (x >> 1) & 0x5555555555555555L;
        //put count of each 4 bits into those 4 bits
        x = (x & 0x3333333333333333L) + ((x >> 2) & 0x3333333333333333L);
        //put count of each 8 bits into those 8 bits
        x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fL;
        //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
        x = ((x * 0x0101010101010101L) >> 56);
        return (int) (x & 0x7f);
    }
    
    /**
     * 
     * @param v
     * @return 
     */
    public static int bitReverse(int v, int nBits) {

        int r = v;
                
        int rev = 0;
        for (int i = 0; i < nBits; i++) {
            rev = (rev << 1) | (r & 1);
            r >>= 1;
        }
        
        return rev;
        
    }
    
       /**
     * method from scipy's numpy package to return numbers spaced evenly 
     * on a log scale.
     <pre>
     numpy is licensed under BSD 3-clause "New" or "Revised" License
      
     https://github.com/numpy/numpy/blob/master/LICENSE.txt
     </pre>
    
    @param start is the starting value of the sequence.
    @param stop float is the final value of the sequence, 
       unless `endpoint` is False.  
       In that case, ``num + 1`` values are spaced over the
        interval in log-space, of which all but the last 
        (a sequence of length ``num``) are returned.
    @param num : Number of samples to generate.  Default is 50.
    @param endpoint : If true, `stop` is the last sample. 
        Otherwise, it is not included.
        Default is True.
    */
    public static TDoubleArrayList logspace(float start, float stop, int num,
        boolean endpoint) {
                
        TFloatArrayList y = linspace(start, stop, num, endpoint);

        TDoubleArrayList logY = new TDoubleArrayList(y.size());

        for (int i = 0; i < y.size(); ++i) {
            float v = y.get(i);
            double l = Math.pow(10, v);
            logY.add(l);
        }
        
        return logY;
    }
    
    /**
     * adapted method from scipy's numpy package to return evenly spaced numbers 
     * over a specified interval
     <pre>
     numpy is licensed under BSD 3-clause "New" or "Revised" License
      
     https://github.com/numpy/numpy/blob/master/LICENSE.txt
     </pre>
    
    @param start is the starting value of the sequence.
    @param stop float is the final value of the sequence, 
       unless `endpoint` is False.  
       In that case, ``num + 1`` values are spaced over the
        interval in log-space, of which all but the last 
        (a sequence of length ``num``) are returned.
    @param num : Number of samples to generate.  Default is 50.
    @param endpoint : If true, `stop` is the last sample. 
        Otherwise, it is not included.
        Default is True.
    */
    public static TFloatArrayList linspace(
        float start, float stop, int num,
        boolean endpoint) {
        
        if (num < 0) {
            throw new IllegalArgumentException("num must be .gte. 0");
        }
        
        int n2 = endpoint ? (num - 1) : num;
        
        float dt = (stop - start)/n2;
        
        TFloatArrayList y = new TFloatArrayList(num);
        for (int i = 0; i <= num; ++i) {
            float v = start + (dt * i);
            y.add(v);
        }
        
        return y;
    }
    
    /**
     * calculates the mean of the data per dimension and returns it as a point 
     * in all dimensions.
     * @param data nDimensional data points in format [ point_0 in all dimensions,
     *   point_1 in all dimensions, ... point_{n-1} in all dimensions[
     * @param nDimensions the number of dimensions of a point
     * @return an array holding the mean per dimension, e.g. point_0_0,
     * point_0_1, point_0_2...point_0_{nDimensions-1}]
     */
    public static double[] mean(double[] data, int nDimensions) {
        
        if ((data.length % nDimensions) != 0) {
            throw new IllegalArgumentException("data.length must be a multiple of nDimensions");
        }
        int nData = data.length/nDimensions;
        
        double[] c = new double[nDimensions];
        int i, j, d;
        for (i = 0; i < nData; ++i) {
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                c[d] += data[j];
            }
        }
        for (d = 0; d < nDimensions; ++d) {
            c[d] /= (double) nData;
        }
        return c;
    }
    /**
     * calculates the mean of the data
     * @param data data points
     * @return the sample mean
     */
    public static double mean(double[] data) {

        int n = data.length;
        int i;
        double sum = 0;
        for (i = 0; i < n; ++i) {
            sum += data[i];
        }
        sum /= (double)n;
        return sum;
    }

    /**
     * calculates the mean per dimension and standard deviation per dimension
     * of the data and returns them in a double array of size [2][nDimensions] 
     * @param data nDimensional data points in format [ point_0 in all dimensions,
     *   point_1 in all dimensions, ... point_{n-1} in all dimensions[
     * @param nDimensions the number of dimensions of a point
     * @return the standard deviation of the data in each dimension.  return format
     * is a double array of size [2][nDimensions] 
     */
    public static double[][] standardDeviation(double[] data, int nDimensions) {
        
        if ((data.length % nDimensions) != 0) {
            throw new IllegalArgumentException("data.length must be a multiple of nDimensions");
        }
        int nData = data.length/nDimensions;
        
        double[] c = mean(data, nDimensions);
        double[][] out = new double[2][];
        out[0] = c;
        out[1] = new double[nDimensions];
        
        int i, j, d;
        double diff;
        for (i = 0; i < nData; ++i) {
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                diff = data[j] - out[0][d];
                out[1][d] += (diff*diff);
            }
        }
        for (d = 0; d < nDimensions; ++d) {
            out[1][d] = Math.sqrt(out[1][d]/(nData - 1.0f)); 
        }
        
        return out;
    }
    
    public static double[] cumulativeSum(double[] a) {
        
        double[] s = Arrays.copyOf(a, a.length);
        
        for (int i = 1; i < a.length; ++i) {
            s[i] += s[i - 1];
        }
        
        return s;
    }
    
    public static int[] cumulativeSum(int[] a) {
        
        int[] s = Arrays.copyOf(a, a.length);
        
        for (int i = 1; i < a.length; ++i) {
            s[i] += s[i - 1];
        }
        
        return s;
    }
    
    public static double[][] cumulativeSumMatlabPort(double[][] a) {
        if (a.length > 1) {
            return cumulativeSumAlongColumns(a);
        } else if (a[0].length > 1) {
            return cumulativeSumAlongRows(a);
        }
        throw new IllegalArgumentException("a must have at least one dimension of size > 1");
    }
    
    public static double[][] cumulativeSumAlongColumns(double[][] a) {
        
        double[][] s = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; ++i) {
            s[i] = Arrays.copyOf(a[i], a[i].length);
        }
        
        for (int row = 1; row < a.length; ++row) {
            for (int col = 0; col < a[0].length; ++col) {
               s[row][col] += s[row - 1][col];
            }
        }
        
        return s;
    }
    
    public static double[][] cumulativeSumAlongRows(double[][] a) {
        
        double[][] s = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; ++i) {
            s[i] = Arrays.copyOf(a[i], a[i].length);
        }
        
        for (int row = 0; row < a.length; ++row) {
            for (int col = 1; col < a[0].length; ++col) {
               s[row][col] += s[row][col - 1];
            }
        }
        
        return s;
    }

    public static TIntIntMap makeFrequencyMap(int[] a) {
        TIntIntMap f = new TIntIntHashMap();
        int c;
        for (int key : a) {
            if (f.containsKey(key)) {
                c = f.get(key) + 1;
            } else {
                c = 1;
            }
            f.put(key, c);
        }
        return f;
    }
    
     /**
     * compute n!/k!(n-k)!.  Note that if n or k are larger than 12,
     * computeNDivKTimesNMinusKBigIntegerExact is used and in that case,
     * if the result is larger than Long.MAX_VALUE an exception is thrown.
     *
     * (Aho & Ullman "Foundations of Computer Science")
     * @param n
     * @param k
     * @return
     * @throws ArithmeticException thrown when result is out of range of type long
     */
    public static long computeNDivKTimesNMinusK(int n, int k) {

        if (n == k) {
            return 1;
        }
        
        if (k > 12 || n > 12) {
            BigInteger result = computeNDivKTimesNMinusKBigInteger(n, k);
            if (result.bitLength() > 63) {
                throw new ArithmeticException("the result will not fit in a long");
            }
            return result.longValue();
        }

        int i;
        double result = 1;
        for (i = n; i > (n-k); i--) {
            result *= i;
            result /= (i - n + k);
        }
        
        return Math.round(result);
    }
    
     /**
     * compute n!/(n-k)!... needed for large numbers
     *
     * @param n
     * @param k
     * @return
     */
    public static long computeNDivNMinusK(int n, int k) {

        if (n == k) {
            return 1;
        }

        long result = 1;
        for (int i = n; i > (n-k); i--) {
            result *= i;
        }
        return result;
    }
    
    /**
     * compute n!/k!(n-k)!
     * @param n
     * @param k
     * @return 
     * @throws ArithmeticException thrown when result is out of range of type long
     */
    public static BigInteger computeNDivKTimesNMinusKBigInteger(int n, int k) {
        
        if (n == k) {
            return BigInteger.ONE;
        }
        
        MathContext ctx = MathContext.DECIMAL128;//MathContext.UNLIMITED;
        
        /*
        following Aho & Ullman "Foundations of Computer Science".
        
        runtime O(k)
        */
        
        BigDecimal result = BigDecimal.ONE;
        int i;
        BigDecimal m;
        
        for (i = n; i > (n - k); i--) {
            m = new BigDecimal(Integer.toString(i), ctx);
            result = result.multiply(m, ctx);
            m = new BigDecimal(Integer.toString(i - n + k), ctx);
            result = result.divide(m, ctx);
        }
        result = result.setScale(0, RoundingMode.HALF_UP);
        return result.toBigInteger();
    }
    
    /**
     * compute n!
     *
     * @param n
     * @return
     */
    public static long factorial(int n) {

        if (n < 3) {
            return n;
        }
        
        if (n > 20) {
            throw new IllegalArgumentException("use factorialBigInteger instead");
        }

        long result = 1;
        for (int i = 2; i <= n; i++) {
            result *= i;
        }
        return result;
    }
    
    /**
     * compute n!
     *
     * @param n
     * @return
     */
    public static BigInteger factorialBigInteger(int n) {

        BigInteger result = BigInteger.ONE;

        for (int i = 2; i <= n; i++) {
            
            byte[] bytes = MiscMath0.writeToBigEndianBytes(i);
            
            BigInteger v = new BigInteger(bytes);
            
            result = result.multiply(v);
        }
        
        return result;
    }
    
    /**
     * write value to a byte array in big endian, that is LSB in highest order bit
     * (MSB is in lowest memory address).
     * these are signed values stored as twos complement and can be input
     * to BigInteger's constructor.
     * @param value
     * @return 
     */
    public static byte[] writeToBigEndianBytes(long value) {
    
        long nBits = numberOfBits(value);
        
        int nBytes = (int) Math.ceil((float)nBits/(float)4);
        
        //System.out.println("nBits=" + nBits + " value=" + value + " nBytes=" + nBytes);
        
        byte[] bytes = new byte[nBytes];

        for (int i = 0; i < nBytes; i++) {
            long shift = i * 8;
            long a = (value >> shift);
            byte b = (byte)a;
            bytes[nBytes - i - 1] = b;
        }

        return bytes;
    }

    public static double[] getMinMax(double[] d) {
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for (double a : d) {
            if (a < min) {
                min = a;
            }
            if (a > max) {
                max = a;
            }
        }
        return new double[]{min, max};
    }
    
    public static float[] getMinMax(float[] d) {
        float min = Float.POSITIVE_INFINITY;
        float max = Float.NEGATIVE_INFINITY;
        for (float a : d) {
            if (a < min) {
                min = a;
            }
            if (a > max) {
                max = a;
            }
        }
        return new float[]{min, max};
    }

    /**
     * third standardized moment of skewness as sample skewness:
     * <pre>
     *  https://en.m.wikipedia.org/wiki/Skewness
     *     g1 = kappa_3/((kappa_2)^(3/2))
     *        = (1/n) * summation from i=1 to n ( (x_i - mean)^3)
     *           / ([ (1/n) * summation from i=1 to n ( (x_i - mean)^2) ]^(3/2) )
     *
     *           where kappa_2 is the 2nd cumulant = variance.
     *          kappa_3 = central moment which is skewness
     * </pre>
     * @param data
     * @return
     */
    public static double calcSampleSkewness(double[] data) {
        double mean = MiscMath0.mean(data);
        int n = data.length;
        double sum2 = 0;
        double sum3 = 0;
        int i;
        double diff;
        double d2;
        for (i = 0; i < n; ++i) {
            diff = data[i] - mean;
            d2 = diff * diff;
            sum2 += d2;
            sum3 += (diff * d2);
        }
        sum3 /= (double)n;
        sum2 /= (double)n;
        sum2 = Math.pow(sum2, 1.5);

        return sum3/sum2;
    }

    public int sign(int v) {
        return v >>> 31;
    }
    
    static int[] MultiplyDeBruijnBitPosition = new int[]{
        0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
        8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
    };
    /**
     * determine the number of bits without branching and using only an int
     * @param v
     * @return 
     */
    public static int numberOfBitsWOB(int v) {
        if (v == 0) {
            return 1;
        }
        //from http://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
        // then edited for negative numbers and signed int
        int sign = v >>> 31;
        // + sign=0 -->  0
        // - sign=1 --> -1
        v = v + sign * (-2) * v;
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        int idx = (v * 0x07C4ACDD) >> 27;
        sign = idx >>> 31;
        //System.out.println("    v=" + v + " sign=" + sign);
        idx += sign*32;
        int r = MultiplyDeBruijnBitPosition[idx];
        return r + 1;
    }
   
    /**
     * http://en.wikipedia.org/wiki/Mersenne_prime
     * @param powerOf2
     * @return 
     */
    public static long getMarsennePrime(int powerOf2) {
        return (1L << powerOf2) - 1;
    }
    
    public static double acosh(double z) {
        //https://mathworld.wolfram.com/InverseHyperbolicCosine.html
        return Math.log(z + Math.sqrt(z + 1) * Math.sqrt(z - 1));
    }

    public static double areaOfTriangle(double[] x, double[] y) {
        if (x.length != 3 || y.length != 3) {
            throw new IllegalArgumentException("x and y must be lengths 3");
        }
        // choose the first x,y pair to be the origin (0, 0)
        double oX = -x[0];
        double oY = -y[0];

        // translate all points by adding xO and oY, then the area is
        double x1 = x[1] + oX;
        double x2 = x[2] + oX;
        double y1 = y[1] + oY;
        double y2 = y[2] + oY;
        // 0.5 * det|x y| = 0.5 * (x1*y2 - x2*y1)
        return 0.5 * ((x1 * y2) - (y1 * x2));
    }
    
    public static double asinh(double z) {
        //https://mathworld.wolfram.com/InverseHyperbolicSine.html
        return Math.log(z + Math.sqrt(z*z + 1));
    }
    
    public static void reverse(int[] a) {
        
        int n = a.length;
        
        if (n < 2) {
            return;
        }
                
        int end = n >> 1;
        int swap, i2;
        // 0 1 2 3 4
        for (int i = 0; i < end; i++) {
            i2 = n - i - 1;
            swap = a[i];
            a[i] = a[i2];
            a[i2] = swap;
        }
    }

    /**
     * return the Euler-Mascheroni constant.
     * It is the first derivative of the gamma function w.r.t. n at n=1.
     * <pre>
     *     Chap 19 of "Statistical Distributions" by Evans et al.
     * </pre>
     * @return
     */
    public static double eulerMascheroniConstant() {
        return 0.5772156649;
    }

    public static float[] convertDoubleToFloat(double[] a) {
        if (a == null) {
            return null;
        }
        float[] b = new float[a.length];
        for (int i = 0; i < a.length; ++i) {
            b[i] = (float)a[i];
        }
        return b;
    }

    public static double[] convertFloatToDouble(float[] a) {
        if (a == null) {
            return null;
        }
        double[] b = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            b[i] = a[i];
        }
        return b;
    }
    public static double[] convertIntToDouble(int[] a) {
        if (a == null) {
            return null;
        }
        double[] b = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            b[i] = a[i];
        }
        return b;
    }
}
